#!/bin/Rscript

library(ggplot2)
library(parallel)

############################################################################################################
working_dir <- "path_to_your_working_dir"
all_RNA_seq_file='./SRR_list.txt'
output_dir <- './output/ABS_event'
cores <- 6
setwd(working_dir)
dir.create(output_dir)

cl<-makeCluster(cores)

## According to the total_filt & hfilt circRNAs lists, identify ABS events
# total_filt circRNAs
# import datasets
# total_filt_circ<-read.csv('./output/circQTLs/log_anno_array_filt_anno.csv',header = T)
# hfilt_circRNAs
hfilt_circ<-read.csv('./output/circQTLs/log_anno_array_hfilt_3soft_anno.csv',header = T)
total_filt_circ<-hfilt_circ
total_filt_info<-total_filt_circ[,c((length(total_filt_circ)-2):length(total_filt_circ))]
total_filt_count<-total_filt_circ[,-c((length(total_filt_circ)-2):length(total_filt_circ))]>0
total_filt_readscount<-total_filt_circ[,-c((length(total_filt_circ)-2):length(total_filt_circ))]
total_filt_info$ABS<-'Not ABS'
total_filt_info$ABS_class<-NA
# ABS identification function
ABS_iden<-function(x,total) {

  A5BS<-total$chrom==unlist(x[1])&total$start==unlist(x[2])
  A3BS<-total$chrom==unlist(x[1])&total$end==unlist(x[3])
  # ABS<-total$chrom==unlist(x[1])&total$start==unlist(x[2])&total$end==unlist(x[3])
  n<-length(total_filt_info[,1])
 
  if (sum(A5BS)>=2 & sum(A3BS)>=2) {
      abs_5b<-total[A5BS,];abs_3b<-total[A3BS,]
      ABS<-'A5BS & A3BS'
      abs_chr<-unlist(x[1]);abs_start<-min(abs_5b$start);abs_end<-max(abs_5b$end);
      abs_3b_start<-min(abs_3b$start);abs_3b_end<-max(abs_3b$end);
      abs<-c(abs_chr,abs_start,abs_end,abs_3b_start,abs_3b_end)
  } else if (sum(A5BS)>=2 & sum(A3BS)==1) {
    abs_each<-total[A5BS,];ABS<-'A5BS'
    A3BS<-rep(FALSE,n)
    abs_chr<-unlist(x[1]);abs_start<-min(abs_each$start);abs_end<-max(abs_each$end);
    abs<-c(abs_chr,abs_start,abs_end)
  } else if (sum(A3BS)>=2 & sum(A5BS)==1) {
    abs_each<-total[A3BS,];ABS<-'A3BS'
    A5BS<-rep(FALSE,n);
    abs_chr<-unlist(x[1]);abs_start<-min(abs_each$start);abs_end<-max(abs_each$end);
    abs<-c(abs_chr,abs_start,abs_end)
  } else {
    abs_each<-total[ABS,];ABS<-'not ABS'
    # exclude itself
    A5BS<-rep(FALSE,n);A3BS<-rep(FALSE,n);
    abs<-NA
  }
  
  return(list(abs,A5BS,A3BS,ABS))
  
}

# ABS identify in for-loop
for (i in 1:length(total_filt_info[,1])) {
  
  x = list(total_filt_info[i,])[[1]]
  res<-ABS_iden(x,total_filt_info)
  abs_info<-res[[1]];
  A5BS<-res[[2]]*(1:length(total_filt_info[,1]));
  A3BS<-res[[3]]*(1:length(total_filt_info[,1]));
  ABS<-res[[4]]
  
  total_filt_info$ABS[i]<-ABS
  if (sum(A5BS)>=2&sum(A3BS)>=2){
    total_filt_info$ABS_class[i]<-paste(abs_info[1],':',abs_info[2],'-',abs_info[3],';',
                                          abs_info[1],':',abs_info[4],'-',abs_info[5],sep='')
    
  } else  {
    total_filt_info$ABS_class[i]<-paste(abs_info[1],':',abs_info[2],'-',abs_info[3],sep='')
    
  }
  # total_filt_info$ABS[unlist(a53bs)]<-'A5BS & A3BS';
  # if (i == 1) {
    # write.table('./output/ABS_event/ABS_events.txt',sep='\t',quote = F,row.names = F,append = F)
  # } else {
    # write.table('./output/ABS_event/ABS_events.txt',sep='\t',quote = F,row.names = F,append = T)
  # }
}
total_filt_abs<-cbind(total_filt_info,total_filt_count)
total_filt_abs_reads<-cbind(total_filt_info,total_filt_readscount)
# write.table(total_filt_abs,'./output/ABS_event/ABS_events.txt',sep='\t',quote = F,row.names = F)
# write.table(total_filt_abs_reads,'./output/ABS_event/ABS_events_readscount.txt',sep='\t',quote = F,row.names = F)
# # hfilt
write.table(total_filt_abs,'./output/ABS_event/ABS_events_hfilt.txt',sep='\t',quote = F,row.names = F)
write.table(total_filt_abs_reads,'./output/ABS_event/ABS_events_readscount_hfilt.txt',sep='\t',quote = F,row.names = F)

## Analysis of ABS events
# total circRNA
total_filt_info<-read.table('./output/ABS_event/ABS_events.txt',sep='\t',header = T)
# total circRNA in matrixeqtl
total_filt_info<-read.table('./output/ABS_event/ABS_events_hfilt.txt',sep='\t',header = T)
rownames(total_filt_info)<-paste(total_filt_info$chrom,':',total_filt_info$start,'-',total_filt_info$end,sep='')

# circRNA associated with circQTL in introns
# read the snps-circRNAs pairs
hvaild_circqtl<-read.table('./output/circQTLs/hvaild_circQTLs.txt',header = T)
# import circQTLs
intron_snps<-read.table('./output/circQTLs/high_vaild_significant_log_quant.txt',header = T)
# intron
intron_snps<-intron_snps[intron_snps$ConsequenceType=='intron_variant',]

# gwas
inter_temp_frep<-read.table('./output/coloc/inter_gwas_loci_circQTL.txt',sep='\t',header = T)
gwas_circqtl<-inter_temp_frep$VariantID
snps<-gwas_circqtl
# head(inter_temp_frep)
# dim(inter_temp_frep)
snps<-intersect(intron_snps$VariantID,hvaild_circqtl$snps)
idx<-unlist(lapply(snps, function(x){which(hvaild_circqtl$snps==x)}))
circ_intron<-unique(hvaild_circqtl$gene[idx])
idx_intron<-unlist(lapply(circ_intron, function(x){which(rownames(total_filt_info)==x)}))
total_filt_info<-total_filt_info[unique(idx_intron),]
# splice-sites (too little circQTLs)
# ss_snps<-read.table('./output/circQTLs/high_vaild_significant_log_quant.txt',header = T)
# ss_snps<-ss_snps[ss_snps$ConsequenceType=='intron_variant',]
# snps<-intersect(intron_snps$VariantID,hvaild_circqtl$snps)
# idx<-unlist(lapply(snps, function(x){which(hvaild_circqtl$snps==x)}))
# circ_intron<-hvaild_circqtl$gene[idx]
# idx_intron<-unlist(lapply(circ_intron, function(x){which(rownames(total_filt_info)==x)}))
# total_filt_info<-total_filt_info[unique(idx_intron),]

############################################################################################################
#####################################
## Figure S11
#####################################
## Fig S11A-C
# plot Number of ABS events in total datasets
library(VennDiagram)
library(reshape2)
A5BS_list<-rownames(total_filt_info)[total_filt_info$ABS=='A5BS']
A3BS_list<-rownames(total_filt_info)[total_filt_info$ABS=='A3BS']
ABS_list<-rownames(total_filt_info)[total_filt_info$ABS=='A5BS & A3BS']
A5BS_list<-append(A5BS_list,ABS_list);A3BS_list<-append(A3BS_list,ABS_list)
venn.diagram(list('A3BS'=A3BS_list,'A5BS'=A5BS_list), fill=c("green","red"), 
             sub='A5BS & A3BS',sub.pos=c(0.5,0.9),sub.cex=2,sub.fontfamily = 'serif',
             alpha=c(0.5,0.5), cex=2, fontfamily='serif',
             cat.cex = 2, cat.fontface = 2,cat.dist=c(0.03,0.03),
            cat.pos=c(20,-20),inverted=T,
             filename="./Distribution of ABS_events in intron.tiff")

############################################################################################################			 
#####################################
## Figure 3
#####################################
## Fig 3A-C
# plot Number or expression of ABS events in different breeds or samples
breed_info = read.csv('./output/matrixeqtl/SRR_breeds.CSV')
breeds = breed_info$Breeds[!duplicated(breed_info$Breeds)]
total_filt_breeds<-total_filt_info[,1:5]
for (i in 1:length(breeds)) {
  # i=1
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  total_filt_breeds[,breeds[i]]<-rowSums(total_filt_info[,srr])
  
}
total_breeds_number<-total_filt_breeds[,c(4,6:length(total_filt_breeds[1,]))]
total_breeds_number<-melt(data = total_breeds_number,id.vars = 'ABS',variable.name = 'Breeds',value.name = 'Number')
total_breeds_number$Number<-total_breeds_number$Number>0

##Fig 3A-C
ABS=pal_nejm()(8)[1:4]
names(ABS)<-unique(annotation_row$ABS)[c(3,1,4,2)]
# plot
ggplot(total_breeds_number[total_breeds_number$Number==TRUE,])+geom_bar(aes(x=Breeds,fill=ABS),stat = 'count')+
  scale_fill_manual(values = ABS)+
  theme_bw()+ylab('Number of circRNA \nwith ABS events')+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=18,angle = 45,hjust=1),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

## Fig 3G Clusterheatmap of circRNA with different ABS events for breeds
total_filt_abs_reads<-read.table('./output/ABS_event/ABS_events_readscount_hfilt.txt',sep='\t',header = T)

# convert individuals to breeds
breed_info = read.csv('./SRR_breeds.csv')
# annotation col
cols<-data.frame("SRR_list"=colnames(total_filt_abs_reads)[6:length(total_filt_abs_reads)])
annotation_col<-merge(cols,breed_info,by.y ='SRR_list');SRR_lisr_filt<-annotation_col$SRR_list
rownames(annotation_col)<-annotation_col$SRR_list
annotation_col<-subset(annotation_col,select = 'Breeds')
# annotation_row
annotation_row<-data.frame("id"=paste0("circ",rownames(total_filt_abs_reads)),"ABS"=total_filt_abs_reads$ABS)
rownames(annotation_row)<-annotation_row$id
annotation_row<-subset(annotation_row,select = 'ABS')
rownames(total_filt_abs_reads)<-rownames(annotation_row)
# annotation_color
ABS=pal_nejm()(8)[1:4]
names(ABS)<-unique(annotation_row$ABS)[c(3,1,4,2)]
Breeds=c("#ff9289","#82b7ff",'#00dca6','#ff83ff','#e199ff',"#d3ba00")
names(Breeds)=breeds
annotation_color=list(
  ABS=ABS,
  Breeds=Breeds
)
annotation_colors = annotation_color
library(pheatmap)
# filt it by the Mean value of log2(CPM)>=5
idx<-which(rowMeans(total_filt_abs_reads[,cols$SRR_list])>=8)
pheatmap(total_filt_abs_reads[idx,cols$SRR_list],scale='column',annotation_row = annotation_row,annotation_col = annotation_col,fontsize=18,
         show_rownames = F,show_colnames = F,annotation_colors = annotation_color)

## Differential analysis of expression for circRNA
library(reshape2)
library(ggrepel)

diff_df<-total_filt_abs_reads[idx,cols$SRR_list]
diff_df<-data.frame(scale(diff_df))
diff_df$circ<-rownames(diff_df);annotation_row$circ<-rownames(annotation_row)
diff_df<-merge(diff_df,annotation_row,by='circ')
diff_df<-melt(diff_df,id.vars = c("circ","ABS"),variable.name = "SRR_list",value.name = "log2CPM")

## ABS vs other ABS in overall
medians<-c(round(median(diff_df_abs$log2CPM[diff_df_abs$isABS=='ABS']),3),
           round(median(diff_df_abs$log2CPM[diff_df_abs$isABS!='ABS']),3))
## Fig 3H
ggplot(diff_df_abs)+geom_boxplot(aes(x=isABS,y=log2CPM,fill=isABS))+coord_flip()+
  # geom_jitter(aes(x=ABS,y=log2CPM,color=ABS),size=1)+
  scale_fill_manual(values = isABScolor)+
  scale_color_manual(values = isABScolor)+
  # stat_compare_means(comparisons = comparisons)+
  # geom_signif(test="wilcox.test", map_signif_level = F)
  stat_compare_means(aes(x=isABS,y=log2CPM,label = ifelse(..p.format.. < 0.05, sprintf("p = %2.1e", as.numeric(..p.format..)),
                                                          sprintf("p = %s",..p.format..))),label.x.npc='center',
                                                  hjust=1.5,vjust = 0.5,size=5)+
  # stat_pvalue_manual(df_wilcox_filt,color ="red",bracket.nudge.y = 0.5,tip.length = 0,step.increase = 0.02)+
  xlab('')+ylab(expression(paste('Z-score of ', log[2],'(CPM)')))+
  annotate('text',x=c(1,2),y=c(4.5,4.2),label=paste0("Median = ",medians),size=4,hjust=0.2)+
  theme_bw()+theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15),
                   strip.text = element_text(size = 18))

## all ABS types vs other ABS in overall
# pairwise comparison
df_wilcox <- diff_df_abs %>%
  pairwise_wilcox_test(log2CPM ~ ABS) %>%
  add_y_position(step.increase = 0.02)

df_wilcox_filt<-df_wilcox[df_wilcox$p.adj<=0.05,]

# medians<-c(round(median(diff_df_abs$log2CPM[diff_df_abs$isABS=='ABS']),3),
#            round(median(diff_df_abs$log2CPM[diff_df_abs$isABS!='ABS']),3))

#####################################
## Figure S12
#####################################
## Fig S12A
ggplot(diff_df_abs)+geom_boxplot(aes(x=ABS,y=log2CPM,fill=ABS))+
  # geom_jitter(aes(x=ABS,y=log2CPM,color=ABS),size=1)+
  scale_fill_manual(values = ABS)+
  scale_color_manual(values = ABS)+
  # stat_compare_means(comparisons = comparisons)+
  # geom_signif(test="wilcox.test", map_signif_level = F)
  stat_compare_means(aes(x=ABS,y=log2CPM),label = 'p',label.x.npc='center',hjust=0.5,vjust = 0,size=5)+
  stat_pvalue_manual(df_wilcox_filt,color ="red",bracket.nudge.y = -2.5,tip.length = 0,step.increase = 0.03)+
  xlab('')+ylab(expression(paste('Z-score of ', log[2],'(CPM)')))+
  # annotate('text',x=c(1,2),y=c(4.5,4.2),label=paste0("Median = ",medians),size=4)+
  theme_bw()+theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(angle = 45,hjust = 1,size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15),
                   strip.text = element_text(size = 18))
				   
############################################################################################################
## Figure 3I
# dim(diff_df)
diff_df<-merge(diff_df,breed_info,by=c('SRR_list'))
## ABS vs not ABS
diff_df_abs<-diff_df
diff_df_abs$isABS<-"not ABS"
diff_df_abs$isABS[diff_df_abs$ABS!="not ABS"]<-"ABS"
#shapiro test
shapiro_res<-shapiro.test(diff_df_abs$log2CPM[1:5000])$p.value
shapiro_res
# pairwise comparison
df_wilcox <- diff_df_abs %>%
  group_by(Breeds) %>%
  pairwise_wilcox_test(log2CPM ~ isABS) %>%
  add_y_position(step.increase = 0.02)

df_wilcox_filt<-df_wilcox[df_wilcox$p.adj<=0.05,]
# set isABS color 
isABScolor<-c(pal_nejm()(8)[5],ABS['not ABS']);names(isABScolor)<-c("ABS","not ABS")

# isABS<-list(isABScolor)
# plotfacet_wrap(~Breeds)+
ggplot(diff_df_abs)+geom_boxplot(aes(x=isABS,y=log2CPM,fill=isABS))+facet_wrap(~Breeds)+
  # geom_jitter(aes(x=ABS,y=log2CPM,color=ABS),size=1)+
  scale_fill_manual(values = isABScolor)+
  scale_color_manual(values = isABScolor)+
  # stat_compare_means(comparisons = comparisons)+
  # geom_signif(test="wilcox.test", map_signif_level = F)
  stat_compare_means(aes(x=isABS,y=log2CPM,label = ifelse(..p.format.. < 0.05, sprintf("p = %2.1e", as.numeric(..p.format..)),
                                                          sprintf("p = %s",..p.format..))),label.x.npc='center',hjust=0.5,vjust = 0.5,size=5)+
  stat_pvalue_manual(df_wilcox_filt,color ="red",bracket.nudge.y = 0.5,tip.length = 0,step.increase = 0.02)+
  xlab('')+ylab(expression(paste('Z-score of ', log[2],'(CPM)')))+
  theme_bw()+theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(angle = 45,hjust = 1,size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15),
                   strip.text = element_text(size = 18))

############################################################################################################
## Figure S12B
## all ABS comparison pairs
# normal distribution test
#shapiro test
shapiro_res<-shapiro.test(diff_df$log2CPM)$p.value
# comparisons<-list(c("A5BS","not ABS"),c("A3BS","not ABS"))
library(rstatix)
df_wilcox <- diff_df %>%
  group_by(Breeds) %>%
  pairwise_wilcox_test(log2CPM ~ ABS) %>%
  add_y_position(step.increase = 0.02)

df_wilcox_filt<-df_wilcox[df_wilcox$p.adj<=0.05,]
# compare_means()
library(ggpubr)
ggplot(diff_df)+geom_boxplot(aes(x=ABS,y=log2CPM,fill=ABS))+facet_wrap(~Breeds)+
  # geom_jitter(aes(x=ABS,y=log2CPM,color=ABS),size=1)+
  scale_fill_manual(values = ABS)+
  scale_color_manual(values = ABS)+
  # stat_compare_means(comparisons = comparisons)+
  # geom_signif(test="wilcox.test", map_signif_level = F)
  stat_compare_means(aes(x=ABS,y=log2CPM),label='p.format',label.x.npc='center',hjust=0.5,vjust = 0.5,size=5)+
  stat_pvalue_manual(df_wilcox_filt,color ="red",step.group.by="Breeds",tip.length = 0,step.increase = 0.02)+
  xlab('')+ylab(expression(paste('Z-score of ', log[2],'(CPM)')))+
  theme_bw()+theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(angle = 45,hjust = 1,size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15),
                   strip.text = element_text(size = 18))

# breed_info<-breed_info[unlist(lapply(breed_info$SRR_list,function(x){x%in%all_gtf})),]
breeds = breed_info$Breeds[!duplicated(breed_info$Breeds)]
total_filt_abs_reads_convert=data.frame(total_filt_abs_reads[,1:5])
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  total_filt_abs_reads_breed = subset(total_filt_abs_reads,select = srr)
  total_filt_abs_reads_convert[,as.character(breeds[i])]=apply(total_filt_abs_reads_breed,1,mean)
}

############################################################################################################
## Fig 3D-F
# plot the cumulative fraction of ABS circRNAs in different expression levels
total_filt_abs_reads<-read.table('./output/ABS_event/ABS_events_readscount.txt',sep='\t',header = T)
# hfilt
total_filt_abs_reads<-read.table('./output/ABS_event/ABS_events_readscount_hfilt.txt',sep='\t',header = T)
rownames(total_filt_abs_reads)<-paste(total_filt_abs_reads$chrom,':',total_filt_abs_reads$start,'-',total_filt_abs_reads$end,sep='')
#intron
total_filt_abs_reads<-total_filt_abs_reads[unique(idx_intron),]
# gwas
total_filt_abs_reads<-total_filt_abs_reads[unique(idx_intron),]

total_filt_abs_total<-cbind(total_filt_abs_reads[,1:5],rowMeans(total_filt_abs_reads[,6:length(total_filt_abs_reads[1,])]))
colnames(total_filt_abs_total)[6]<-'MeanExpre'
types<-total_filt_abs_total$ABS[!duplicated(total_filt_abs_total$ABS)]
for (i in 1:length(types)) {
  # i=1
  tmp<-total_filt_abs_total[total_filt_abs_total$ABS==types[i],]
  tmp<-tmp[order(tmp$MeanExpre),]
  idx<-rownames(tmp)
  total_filt_abs_total[idx,'Frac']<-seq(0+1/length(tmp[,1]),1,1/length(tmp[,1]))
}
# ggplot(total_filt_abs_total)+geom_histogram(aes(x=))
ggplot(total_filt_abs_total)+geom_line(aes(x=MeanExpre,y=Frac,color=ABS),size=1.5)+
  scale_color_manual(values = ABS)+
  theme_bw()+xlab('Mean expression of circRNAs \nassociated with circQTLs (log2(CPM))')+ylab('Cumulative frequency of circRNAs')+
  theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15))

############################################################################################################
## Fig S11D-F
## plot the ABS_event
library(Gviz)
library(lattice)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(biomaRt)

txdb<-TxDb.Btaurus.UCSC.bosTau9.refGene
ensembl = useMart("ensembl",dataset="btaurus_gene_ensembl")

# construct the object in Track of Gviz
# There are some problems in here!!
abs_vis<-function(total_filt_info,abs_info) {
  library(stringr)
  # abs_info<-'chr27:6220180-6221115'
  len<-str_count(strsplit(abs_info,':',fixed = T)[[1]][1])
  trans_loc<-c(substr(strsplit(abs_info,':',fixed = T)[[1]][1],4,len),
               strsplit(strsplit(abs_info,':',fixed = T)[[1]][2],
                        '-',fixed = T)[[1]][1],
               strsplit(strsplit(abs_info,':',fixed = T)[[1]][2],
                        '-',fixed = T)[[1]][2])
  print(trans_loc)
  # trans_loc<-c('1',"21171970","21181886")
  tran_info<-getBM(attributes = c('chromosome_name','transcript_start','transcript_end','transcript_biotype',
                                  'entrezgene_id','external_gene_name','ensembl_gene_id'),filters = c('chromosome_name','start','end'),
                   values = list(chrom=trans_loc[1],start=trans_loc[2],
                                 end=trans_loc[3]),mart = ensembl)
  # tran_info<-getBM(attributes = c('chromosome_name',
  #                                 'external_gene_name'),filters = c('chromosome_name','start','end'),
  #                  values = list(chrom=trans_loc[1],start=trans_loc[2],
  #                                end=trans_loc[3]),mart = ensembl)
  
  # exon_info_filt<-0
  for (j in 1:length(na.omit(unique(tran_info$entrezgene_id)))) {
    # j=1
    exon_info<-exons(txdb,filter=list(gene_id=unique(tran_info$entrezgene_id)[j]))
    if (nrow(data.frame(exon_info))==0) {
      exon_info_filt<-exon_info
      next
    }
    exon_info<-data.frame(exon_info)
    exon_info_filt<-exon_info[exon_info$start>=as.numeric(trans_loc[2])&
                                exon_info$end<=as.numeric(trans_loc[3]),]
    if (nrow(data.frame(exon_info_filt))==0) {
      break
    }
    exon_info_filt$gene<-unique(tran_info$entrezgene_id)[j];exon_info_filt$transcript<-unique(tran_info$entrezgene_id)[j];
    exon_info_filt$symbol<-unique(tran_info$external_gene_name)[j];
    
  }
  if (nrow(data.frame(exon_info_filt))==0) {
    symbol<-unique(tran_info$external_gene_name)[j]
    return(symbol)
    # pass
  } else {
  symbol<-unique(exon_info_filt$symbol)
  # exon_info_filt<-subset(exon_info_filt,select=-c(exon_))
  # geneModels<-total_filt_info[total_filt_info$ABS_class==abs_class[i],1:3]
  # colnames(geneModels)[1]<-c('chromosome')
  # geneModels$strand<-'+'
  gen <- 'bosTau9'
  chr <- substr(as.character(unique(exon_info_filt$seqnames)),4,4)
  
  gtrack <- GenomeAxisTrack()
  # Exons track
  grtrack <- GeneRegionTrack(exon_info_filt,genome = gen,chromosome = chr, name = "Gene Model",
                             transcriptAnnotation = "symbol", showExonId=F,showId=TRUE,geneSymbol=TRUE)
  
  
  alTrack <-AlignmentsTrack("./Alignment/ERR3555854_quant.sorted.bam", chromosome=chr,from=as.numeric(trans_loc[2]), 
                            to = as.numeric(trans_loc[3]),isPaired = T,options(ucscChromosomeNames=F))
  
  idx<-unlist(lapply(total_filt_info$ABS_class,function(x){is.na(str_match(x,abs_info))}))
  introns <- GRanges(chr,IRanges(start = na.omit(total_filt_info$start[idx]), 
                                 end = na.omit(total_filt_info$end[idx])))
  png(paste0('./ABS_event/output_each/',stringr::str_replace(abs_info,':','-'),'.png'),width = 2400, height = 2000, units = "px",res = 300)
  plotTracks(list(gtrack,alTrack,grtrack),
             chromosome =chr,type=c('coverage',"sashimi"),from = as.numeric(trans_loc[2])-100,to=as.numeric(trans_loc[3])+100,
             col.sashimi='red',lwd.sashimiMax=4,
             sashimiFilter = introns,sashimiFilterTolerance=1e6L,size=c(10,5,5,10),fontsize=15,fontfamily="serif",
             col='black',fontcolor='black',col.border.title='black',fill.confint='black',fontcolor.legend='black',
             fontcolor.item='black',fontcolor.group='black',col.axis=c('black'),col.border.title="white",col.symbol='black',
             col.title="black" ,fontsize.title=15,fontsize.symbol=15,fontsize.group=15,
             collapse=F,innerMargin=0,margin = 0)
  dev.off()
  return(symbol)}
}
# plot the combined tracks of ABS events
abs_class<-na.omit(replace(unique(total_filt_info$ABS_class[total_filt_info$ABS!='A5BS & A3BS']),2,NA))
for (i in 1:length(abs_class)) {
  abs_info<-abs_class[i]
  abs_vis(total_filt_info,abs_info)
}

############################################################################################################
# Comparison of the complexity of ABS events
library(stringr)
total_complex_info<-total_filt_info[!is.na(total_filt_info$ABS_class),1:5]
total_complex_data<-data.frame('ABS_class'=abs_class)
total_complex_info$count<-1
for (i in 1:length(abs_class)) {
  abs_info<-abs_class[i]
  a5bs<-sum(total_complex_info$count[total_complex_info$ABS_class==abs_class[i]&
                                       total_complex_info$ABS=='A5BS'])
  a3bs<-sum(total_complex_info$count[total_complex_info$ABS_class==abs_class[i]&
                                        total_complex_info$ABS=='A3BS'])
  
  idx_5<-!is.na(str_match(total_complex_info$ABS_class,paste0(abs_info,';')))
  a53bs_5<-sum(total_complex_info$count[idx_5&
                                        total_complex_info$ABS=='A5BS & A3BS'])
 
  idx_3<-!is.na(str_match(total_complex_info$ABS_class,paste0(';',abs_info)))
  
  a53bs_3<-sum(total_complex_info$count[idx_3&
                                          total_complex_info$ABS=='A5BS & A3BS'])
  idx<-!is.na(str_match(total_complex_info$ABS_class,paste0(abs_info)))
  a53bs<-sum(total_complex_info$count[idx&
                                        total_complex_info$ABS=='A5BS & A3BS'])
  total_complex_data$A5BS[i]<-a5bs+a53bs_5;total_complex_data$A3BS[i]<-a3bs+a53bs_3;
  total_complex_data$ABS[i]<-a53bs
}

total_complex_data<-melt(total_complex_data,id.vars = 'ABS_class',variable.name = 'ABS_type',
                         value.name = 'Complexity')

total_complex_data$Complexity<-as.numeric(total_complex_data$Complexity)
total_complex_data<-total_complex_data[total_complex_data$Complexity>=2,]
# plot the distribution of number of ABS events with different complexity
# Most complex ABS events
most_complex<-total_complex_data
most_complex_a5bs<-most_complex[most_complex$ABS_type=='A5BS'&most_complex$Complexity==12,]
most_complex_a3bs<-most_complex[most_complex$ABS_type=='A3BS'&most_complex$Complexity==13,]
mc_a5bs<-most_complex_a5bs$ABS_class;mc_a3bs<-most_complex_a3bs$ABS_class
mc_abs<-most_complex[most_complex$ABS_type=='ABS'&most_complex$Complexity==8,]$ABS_class
a5bs_symbol<-list()
for (i in 1:length(mc_a5bs)) {
  a5bs_symbol[[i]]<-abs_vis(total_filt_info,abs_info=mc_a5bs[i])
  
}
a3bs_symbol<-list()
for (i in 1:length(mc_a3bs)) {
  a3bs_symbol[[i]]<-abs_vis(total_filt_info,abs_info=mc_a3bs[i])
  
}
abs_symbol<-abs_vis(total_filt_info,abs_info=mc_abs)
# abs_symbol<-
a5bs_symbol<-abs_vis(total_filt_info,abs_info=mc_a5bs)
a5bs_symbol<-'TUT4'
a3bs_symbol<-c('ADGRA3','TTN')

############################################################################################################
#####################################
## Figure S13
#####################################
library(ggbio)

## Figure S13A
ggplot(total_complex_data)+geom_histogram(aes(x=Complexity,fill=ABS_type),
                                          position = 'dodge',stat = 'count')+
  theme_bw()+ylab('Count')+
  theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15))+
  scale_fill_manual(values = ABS)+
  annotate('text',x=7.0,y=100,label=a5bs_symbol,size=5,color=ABS[2])+
  annotate('segment', x=7.7, xend=7.0, y=0, yend=80,color=ABS[2],
           arrow=arrow(angle = 0,length = unit(0,'cm')))+
  annotate('segment', x=8.0, xend=8.2, y=0, yend=180,color=ABS[3],
          arrow=arrow(angle = 0,length = unit(0,'cm')))+
  annotate('text',x=8.2,y=180,label=a3bs_symbol[1],size=5,color=ABS[3])+
  annotate('segment', x=8.0, xend=7.5, y=0, yend=220,color=ABS[3],
           arrow=arrow(angle = 0,length = unit(0,'cm')))+
  annotate('text',x=7.5,y=220,label=a3bs_symbol[2],size=5,color=ABS[3])
  # annotate('segment', x=8.0, xend=6.7, y=0, yend=280,color='green',
  #          arrow=arrow(angle = 0,length = unit(0,'cm')))+
  # annotate('text',x=8.0,y=280,label=a3bs_symbol[3],size=5,color='green')

# Function of plot ABS pattern
abs_pattern_plot<-function(data,type,xlabel,title,symbol,line.color, curve.color,text.color, title.color,
                           line.size,curve.size, text.size, title.size){
  # Conversion the locations (including chrom, start and end) to visulization plots
  # data: data.frame class, it must include chrom, start and end columns.
  # title: auto or self-designed
  # color: includes line.color, text.color, title.color
  # size: line.size, text.size, title.size
  
  library(ggplot2)
  library(stringr)
  
  fancy_scientific <- function(l) {
   
    l <- l/1e6
    # return this as an expression
    parse(text=l)
  }
  
  data_info<-subset(data,select=c(chrom,start,end))
  chrom<-unique(data$chrom);min_start<-min(data$start);max_end<-max(data$end)
  if (is.null(title)) {
    title<-paste0(chrom,':',min_start,'-',max_end)
  }
  
  # circ_gene<-GRanges(rownames(data))
  circ_gene<-GRanges(data)
  circ_gene$symbol<-symbol
  base <-autoplot(txdb, circ_gene,
             coord = "genome",
             layout = "linear",color=border.color,
             ratio = 0.0025,
             mode = "reduce",geom =
               "alignment",fill=rect.color,gap.geom = "chevron",label.size=12,
             resize.extra = 10, space.skip = 1,names.expr = symbol,label = T,facetByRow=T,genome.axis=F)+
    xlim(min(data$start)-1e4,max(data$end)) + ylim(0,2)+ theme_clear()+
    theme(
      axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),axis.line.y = element_blank(),
    axis.text.x = element_blank(),plot.caption = element_blank(),
    plot.title = element_text(family = 'serif',hjust = 0.5,vjust=-1,color = title.color,size = title.size),
    plot.subtitle = element_text(family = 'serif',hjust = 0.5,vjust=-3,size = title.size)
    )
 
  start<-data$start;end<-data$end
  p<-base+
    geom_curve(aes(x=start, xend=end, y=1, yend=1),
               lineend='round',curvature=0.3,
               color=curve.color,size=curve.size)
  
  ## add transcript by each
  plist<-list()
  for (i in 1:length(circ_gene$ABS)){
    # i=1
    circ_gene_each<-circ_gene[i,]
    symbol_each<-circ_gene$pair[i]
    if (i!=length(circ_gene$ABS)){
      plist[[i]]<-autoplot(txdb, circ_gene_each,coord = "genome",
                           layout = "linear",color=border_each.color,
                           ratio = 0.0025,
                           mode = "reduce",geom = "alignment",fill=rect_each.color,label.size=5,
                           resize.extra = 10, space.skip = 1,label = F,facetByRow=T,genome.axis=F)+
        xlim(min(data$start)-1e4,end(circ_gene_each)) + ylim(0.75,1.25)+
        annotate('text',x =  min(data$start)-1e4, y = 1.0,label=symbol_each)+
        # ylab(symbol_each)
        theme_clear()+
        # xlab(xlabel)+
        theme(
          axis.title.x = element_blank(),axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),axis.line.y = element_blank(),
          axis.text.x = element_blank(),plot.caption = element_blank(),
          
        )
    } else {
      fancy_scientific <- function(l) {
        l <- l/1e6
        # return this as an expression
        parse(text=l)
      }
      
      plist[[i]]<-autoplot(txdb, circ_gene_each,coord = "genome",
                           layout = "linear",color=border_each.color,
                           ratio = 0.0025,
                           mode = "reduce",geom = "alignment",fill=rect_each.color,label.size=5,
                           resize.extra = 10, space.skip = 1,label = F,facetByRow=T,genome.axis=F)+
        xlim(min(data$start)-1e4,end(circ_gene_each)) + ylim(0.75,1.25)+
        annotate('text',x =  min(data$start)-1e4, y = 1.0,label=symbol_each)+
        xlab(xlabel)+
        theme_clear()+
        theme(
          axis.title.x = element_text(size = 25),
          axis.title.y = element_blank(),axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),axis.line.y = element_blank(),
          axis.text.x = element_text(size = 18),plot.caption = element_blank(),
        )+scale_x_continuous(labels=fancy_scientific) 
    }
    
  }
  tracks(p,plist[[1]],plist[[2]],plist[[3]],plist[[4]],plist[[5]],plist[[6]],plist[[7]],plist[[8]],heights = c(3,rep(1,7),1),
         xlab = 'Position on chr3 (bp)',label.text.cex=18,
         theme = theme(
           # axis.text.x = element_text(size=18),
           axis.title.x = element_text(size = 25),
         ))

}

## Pattern plots of ABS events
type<-'A5BS'
symbol<-a3bs_symbol[1]
abs_info<-mc_a3bs[1]
symbol<-a5bs_symbol[1]
abs_info<-mc_a5bs[1]

idx<-unlist(lapply(total_complex_info$ABS_class,function(x){grepl(abs_info,x)}))
data<-total_complex_info[idx,]
data$ABS_class<-paste0(data$ABS_class,':*')
data$pair<-paste0(symbol,"circRNA_",seq(1,length(data$chrom)))

rect.color="#a00000";border.color="#a00000";rect_each.color<-'grey';border_each.color<-'black'
line.color='black'; curve.color='#a00000';text.color='black'; title.color='black';
line.size=2; curve.size=1.0; text.size=10; title.size=30;xlabel="Position on chr3 (Mb)";
## Figure S13B
abs_pattern_plot(data,type,xlabel="Position on chr3 (Mb)",title = NULL,symbol,line.color,curve.color,text.color, title.color,
                 line.size,curve.size, text.size, title.size)
