##!/bin/Rscript

library(ggplot2)
library(stringr)
library(ggsci)
library(ggrepel)
library(reshape2)
############################################################################################################
working_dir <- "path_to_your_working_dir"
circRNAprofiler_dir<-paste(working_dir,'output/RBPpred/predCirc',sep='/')
all_RNA_seq_file='./SRR_list_CIRI.txt'
output_dir <- paste(working_dir,'./output/circQTL_effect',sep='/')
setwd(working_dir)
dir.create(output_dir)
cores<-3
cl<-makeCluster(cores)

## import snps or circQTL within circRNA data
# import  all circRNA internal seq info
load(paste0(circRNAprofiler_dir,'/targetsFTS_circ.Rdata'))

labels<-c("Total snps","circQTLs","GWAS-related circQTLs")
# import location info of different type SNPs
index<-2
snps_in_circRNA<-paste0(labels[index],"_in_circRNA.txt")
all_snps_in_circ<-read.table(snps_in_circRNA,header = T,sep="\t")
all_snps_in_circ<-melt(all_snps_in_circ,id.vars = "id",variable.name = "VariantID",value.name = "isin")
all_snps_in_circ_filt<-all_snps_in_circ[!is.na(all_snps_in_circ$isin),]
head(all_snps_in_circ_filt)
dim(all_snps_in_circ_filt)
## import all loci info
all_snps<-read.table('./output/combined/combined.vcf',header = T)
head(all_snps)
# import the correct ref and alt Alleles 
all_snps_ref_alt<-read.table(paste0(circRNAprofiler_dir,'/combined_data_ref_alt.vcf'),sep="\t",header=F,comment.char = "#")
colnames(all_snps_ref_alt)[1:5]<-c("Chrom","ARS_UCD1_2_Pos","VariantID","Ref","Alt")
all_snps_ref_alt<-all_snps_ref_alt[,1:5]
all_snps_ref_alt$ARS_UCD1_2_Pos<-paste0(all_snps_ref_alt$Chrom,":",all_snps_ref_alt$ARS_UCD1_2_Pos)
head(all_snps_ref_alt)
# circQTLs
circQTLs<-read.table('./output/circQTLs/high_vaild_significant_log_quant.txt',header = T)
head(circQTLs)

# circQTL_loc<-merge(all_snps)
dim(circQTLs)
# GWAS-loci related circQTLs
gwas_circqtls<-read.table('./output/coloc/intersect.txt',sep='\t',header = T)
# gwas_circqtls<-read.table('./output/coloc/eqtl.txt',sep='\t',header = T)
head(gwas_circqtls)
dim(gwas_circqtls)
head(all_snps)
dim(all_snps)
gwas_circQTL<-merge(gwas_circqtls,all_snps,by='VariantID')
head(gwas_circQTL)
dim(gwas_circQTL)
snps<-list(all_snps,circQTLs,gwas_circQTL)

## convert circRNA internal seq to mutation type
all_snps_in_circRNA_genotype<-merge(all_snps_in_circ_filt,all_snps_ref_alt,by="VariantID")
head(all_snps_in_circRNA_genotype)
all_snps_in_circRNA_genotype<-all_snps_in_circRNA_genotype[all_snps_in_circRNA_genotype$isin!="NA,NA",]
## convert function
alt_genotype_circseq<-targetsFTS_circ$circ
base<-c("A","U","C","G")
names(base)<-c("U","A","G","C")
convert_seq<-function(x){
  id<-x[1];mu_loc<-x[2];ref<-x[3];alt<-x[4]
  if (ref == "T"){
    ref<-"U"
  }
  if (alt == "T"){
    alt<-"U"
  }
  
  ref_seq<-targetsFTS_circ$circ$seq[targetsFTS_circ$circ$id==id]
  mu_loc<-strsplit(mu_loc,",",fixed = T)[[1]]
  for (i in 1:length(mu_loc)){
    loc<-as.numeric(mu_loc[i])
    raw_geno<-str_sub(ref_seq,loc,loc)
    if (length(raw_geno)==0){
      print("Not find!")
      next
    }
    if (!is.na(raw_geno) & raw_geno!=''){
      if (raw_geno==ref){
        # str_sub(ref_seq,mu_loc,mu_loc)
        ref_seq<-strsplit(ref_seq,'')[[1]]
        ref_seq[loc]<-alt
        alt_seq<-paste0(ref_seq,collapse = '')
        print('ref-->alt')
        return(alt_seq)
      } else if(raw_geno==names(base)[base==ref]){
        ref_seq<-strsplit(ref_seq,'')[[1]]
        ref_seq[loc]<-names(base)[base==ref]
        alt_seq<-paste0(ref_seq,collapse = '')
        print('ref-->alt')
        return(alt_seq)
      } else if(raw_geno==alt){
        ref_seq<-strsplit(ref_seq,'')[[1]]
        ref_seq[loc]<-ref
        alt_seq<-paste0(ref_seq,collapse = '')
        print('alt-->ref')
        return(alt_seq)
      } else if(raw_geno==names(base)[base==alt]){
        ref_seq<-strsplit(ref_seq,'')[[1]]
        ref_seq[loc]<-names(base)[base==ref]
        alt_seq<-paste0(ref_seq,collapse = '')
        print('alt-->ref')
        return(alt_seq)
      } else{
        return(ref_seq)
      }
    } else{
      return(ref_seq)
    }
    
  }
  
}

alt_genotype_circseq
all_snps_in_circRNA_genotype
res<-apply(all_snps_in_circRNA_genotype[,c("id","isin","Ref","Alt")],1,function(x){return(convert_seq(x))})
# length(unlist(res))
convert_id<-all_snps_in_circRNA_genotype$id
alt_genotype_circseq$alt_seq<-alt_genotype_circseq$seq
alt_genotype_circseq$alt_seq[unlist(lapply(convert_id,function(x){which(alt_genotype_circseq$id==x)}))]<-unlist(res)
alt_genotype_circseq$isconvert<-'ref'
alt_genotype_circseq$isconvert[which(alt_genotype_circseq$seq!=alt_genotype_circseq$alt_seq)]<-'alt'

## export data
write.table(alt_genotype_circseq,'./output/circQTL_effect/alt_genotype_circseq.txt',quote = F,row.names = F,sep='\t')

## RBP motif search
# read.table('./alt_genotype_circseq.txt',header = T)
annotatedBSJs<-read.table('./output/RBPpred/predCirc/annotatedBSJs.txt',header = T,sep='\t')

targets_alt_FTS_circ$circ$seq<-alt_genotype_circseq$alt_seq
load("./output/RBPpred/predCirc/motifsRFTS_circ.Rdata")
mergedMotifsRFTS_circ <- mergeMotifs(motifsRFTS_circ)
motifs_alt_FTS_circ <-
  getMotifs(targets_alt_FTS_circ,
            width = 6,
            database = 'ATtRACT',
            species = "Btarurs",
            rbp = TRUE,
            reverse = FALSE) 

mergedMotifs_alt_FTS_circ <- mergeMotifs(motifs_alt_FTS_circ)
save(motifsFTS_circ,file = './output/RBPpred/predCirc/motifsFTS_circ_alt.Rdata')

#######################################################################################################################################
#####################################
## Figure S19
#####################################
## miRNA binding sites of ref vs alt
library(dplyr)
miRsites_circ <- read.table('./output/RBPpred/predCirc/all_miRNA_binding_sites_internal_by_targetscan-1.txt',header=T,sep='\t')
miRsites_circ_alt<-read.table('./output/RBPpred/predCirc/all_miRNA_binding_sites_alt_internal_by_targetscan.txt',header = T,sep='\t')
miRsites_circ_ref<-merge(miRsites_circ,miRsites_circ_alt[,c("a_Gene_ID","miRNA_family_ID")],by=c("a_Gene_ID","miRNA_family_ID"))
miRsites_circ_ref<-miRsites_circ_ref[!duplicated(miRsites_circ_ref),]
dim(miRsites_circ_ref)
dim(miRsites_circ_alt)
miRsites_circ_alt$genotype<-"alt"
miRsites_circ_ref$genotype<-"ref"
miRsites_circ_all<-rbind(miRsites_circ_ref,miRsites_circ_alt)

site_type<-unique(miRsites_circ_all$Site_type)

df_show<-miRsites_circ_all%>%group_by(genotype,Site_type)%>%summarize(count=n())%>%
  mutate(freq = count / sum(count))
df_show$genotype<-factor(df_show$genotype,levels=rev(unique(df_show$genotype)))
df_show<-df_show[order(df_show$freq,decreasing = F),]
df_show<-df_show%>%mutate(cumfreq=cumsum(freq))
# df_show$freq<-1
##### Figure . Distribution of miRNA site types between ref and alt genotype
# chisqure test
df_show_test<-dcast(df_show,Site_type~genotype,value.var = "count")
rownames(df_show_test)<-df_show_test$Site_type
chiq_res<-chisq.test(df_show_test[,2:3])
test_title<-paste0(chiq_res$method,":",round(chiq_res$p.value,3))
#####################################
## Figure S19A
#####################################
ggplot(df_show,aes(x=genotype,y=freq,fill=Site_type))+geom_bar(stat = 'identity',position='stack')+
  theme_bw()+ylab("Frequency")+
  scale_y_continuous(expand = c(0,0))+
  # scale_x_continuous(expand = c(0,0))+
  geom_text(aes(y=cumfreq-freq/2,label=paste0(round(freq*100,digits = 2),"%")),color='white',size=5)+
  ggtitle(test_title)+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(12)))

#####################################
## Figure S19B
#####################################
## Stack bar plot of miRNA binding sites among breeds
# import breed info
# breed info
breed_info = read.csv('G:/LHF/bovine circQTL/SNP data/output/matrixeqtl/SRR_breeds.CSV')
# breed_info<-breed_info[unlist(lapply(breed_info$SRR_list,function(x){x%in%all_gtf})),]
breeds = breed_info$Breeds[!duplicated(breed_info$Breeds)]
col=breed_info$SRR_list
## imprt circRNA expr data

circRNA_expr<-read.table('G:/LHF/bovine circQTL/SNP data/ABS_event/ABS_events_readscount.txt',header=T,sep = '\t')
circRNA_expr$circ<-paste0(circRNA_expr$chrom,':',circRNA_expr$start,'-',circRNA_expr$end)
## miRNA binding sites
split_loc_circprofile<-function(x){
  chr<-strsplit(x,':',fixed = T)[[1]][3]
  start<-strsplit(x,':',fixed = T)[[1]][4]
  end<-strsplit(x,':',fixed = T)[[1]][5]
  loc<-as.numeric(c(start,end))
  return(paste0(chr,':',min(loc),'-',max(loc)))
}
miRNA_count<-miRsites_circ_all%>%group_by(genotype,a_Gene_ID,Site_type)%>%summarize(count=n())
miRNA_count$circ<-unlist(lapply(miRNA_count$a_Gene_ID,split_loc_circprofile))
# miRNA_count$a_Gene_ID<-unlist(lapply(df_show$a_Gene_ID,function(x){strsplit(x,':',fixed = T)[[1]][1]}))

miRNA_count<-dcast(miRNA_count, genotype+circ~Site_type, fun.aggregate = sum,
                   subset = NULL, drop = TRUE, value.var = 'count')
genotypes<-unique(miRNA_count$genotype)

for (j in 1:length(genotypes)) {
  miRNA_count_temp<-miRNA_count[miRNA_count$genotype==genotypes[j],]
  idx=unlist(lapply(miRNA_count_temp$circ,function(x){which(circRNA_expr$circ==x)}))
circRNA_expr_miR<-circRNA_expr[idx,]
# circRNA_expr_miR<-rbind(circRNA_expr[idx,],circRNA_expr[idx,])
aidx=unlist(lapply(circRNA_expr_miR$circ,function(x){which(miRNA_count_temp$circ==x)}))
miRNA_count_temp<-miRNA_count_temp[aidx,]
# circRNA_expr_miR<-rbind(circRNA_expr_miR,circRNA_expr_miR)
miRNA_count_matrix<-crossprod(as.matrix(miRNA_count_temp[,3:6]),as.matrix(circRNA_expr_miR[,breed_info$SRR_list]!=0))
# merge breed by mean
miRNA_count_convert<-data.frame('Site_type'=rownames(miRNA_count_matrix))
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  miRNA_count_breed = subset(miRNA_count_matrix,select = srr)
  miRNA_count_convert[,as.character(breeds[i])]=apply(miRNA_count_breed,1,mean)/sum(apply(miRNA_count_breed,1,mean))
  
}
df_show<-melt(miRNA_count_convert,id.vars = 'Site_type',variable.name = 'breeds',value.name = 'count')
df_show$genotype<-genotypes[j]
if(j==1){
  df_show_total<-df_show
}else{
  df_show_total<-rbind(df_show_total,df_show)
}
}
ggplot()+
  geom_col(data=df_show_total,
           aes(x = breeds, y = count, group=genotype,fill=Site_type),position = "dodge",
           width=0.5, colour="grey40", size=0.4) 

ggplot()+
  geom_col(data=df_show_total %>% filter(genotype=="ref"),
           aes(x = breeds, y = count, fill=Site_type),
           width=0.5, colour="grey40", size=0.4) +
  geom_col(data=df_show_total %>% filter(genotype=="alt"),
           aes(x = breeds, y = count, fill=Site_type), 
           width=0.5, colour="grey40", size=0.4)
df_show_total$breeds_geno<-paste0(df_show_total$breeds,"_",df_show_total$genotype)
df_show_total<-df_show_total[order(df_show_total$count,decreasing = F),]
df_show_total=df_show_total%>%group_by(breeds,genotype)%>%mutate(cumfreq=cumsum(count))
## chiq test
df_show_test<-dcast(df_show_total,Site_type+breeds~genotype,value.var = "count")
rownames(df_show_test)<-df_show_test$Site_type
chiq_res<-chisq.test(df_show_test[,2:3])
test_title<-paste0(chiq_res$method,":",round(chiq_res$p.value,3))
ggplot(df_show_total,aes(x=breeds_geno,y=count,fill=Site_type,group=genotype))+
  geom_bar(stat = 'identity',position="stack",color='white')+
  # scale_y_discrete(breaks=seq(1.5,12,2),labels=breeds)+
  theme_bw()+ylab("Frequency")+xlab("")+
  # scale_y_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(y=cumfreq-count/2,label=paste0(round(count*100,digits = 2),"%")),color='white',size=5)+
  # geom_text(aes(label=paste0(round(count*100,digits = 2),"%")),color='white',size=5)+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(12)))

#####################################
## Figure S19C
#####################################
## geom_statpie
library(scatterpie)
df_show_pie<-dcast(df_show_total,genotype+breeds~Site_type,value.var = "count")
df_show_pie$x<-unlist(lapply(df_show_pie$breeds,function(x){which(breeds==x)}))
df_show_pie$y<-unlist(lapply(df_show_pie$genotype,function(x){which(genotypes==x)}))
df_show_pie$group<-as.numeric(interaction(df_show_pie$x,df_show_pie$y))
ggplot()+geom_scatterpie(data=df_show_pie,aes(x=x,y=y, group=group),pie_scale = 4.5,cols=site_type)+
  coord_equal()+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(12)))+
  geom_text(data=df_show_pie[!duplicated(df_show_pie$breeds),],aes(x=x,y=y-0.75,label=breeds),size=10)+
  geom_text(data=df_show_pie[!duplicated(df_show_pie$genotype),],aes(x=x-0.75,y=y,label=genotype),size=10)+
geom_label(data=df_show_pie[!duplicated(df_show_pie$breeds),],aes(x=x,y=y+1.6,label="ns"))

## aligned Score and mfe between ref and alt
miRbind_data<-read.table('./test_bySW-1.txt',sep='\t',header = T)
miRbind_data_alt<-read.table('./alt_test_bySW.txt',sep='\t',header = T)
# miRsites_circ_ref
miRsites_circ_ref$miRNA_name<-paste0('bta-',miRsites_circ_ref$miRNA_family_ID)
colnames(miRbind_data)<-c('miRNA_name','end2','subseqDP','a_Gene_ID','mfe','alignedScore')
colnames(miRbind_data_alt)<-c('miRNA_name','end2','subseqDP','a_Gene_ID','mfe','alignedScore')
miRbind_data_ref_filt<-merge(miRbind_data,miRsites_circ_ref[,c('a_Gene_ID','miRNA_name','Site_type','UTR_start','UTR_end')],
                         by=c('miRNA_name','a_Gene_ID'))
miRbind_data_alt_filt<-merge(miRbind_data_alt,miRsites_circ_ref[,c('a_Gene_ID','miRNA_name','Site_type','UTR_start','UTR_end')],
                             by=c('miRNA_name','a_Gene_ID'))
miRbind_data_alt_filt<-merge(miRbind_data_alt_filt,miRbind_data_ref_filt[,c('miRNA_name','a_Gene_ID')],by=c('miRNA_name','a_Gene_ID'))
miRbind_data_alt_filt<-miRbind_data_alt_filt[!duplicated(miRbind_data_alt_filt),]
miRbind_data_ref_filt$genotype<-"ref";miRbind_data_alt_filt$genotype<-"alt"
miRbind_data_combine<-rbind(miRbind_data_ref_filt,miRbind_data_alt_filt)

# dim(miRbind_data_alt_filt)
# sum(duplicated(miRbind_data_alt_filt))
miRbind_data_combine$alignedScore<-unlist(lapply(miRbind_data_combine$alignedScore,function(x){
  max(as.numeric(strsplit(x," ",fixed = T)[[1]]))
}))
df_show<-melt(miRbind_data_combine,id.vars = c(colnames(miRbind_data_combine)[1:4],'Site_type','UTR_start','UTR_end',"genotype"),
              variable.name = 'Index',value.name = 'Score')
# df_show_filt<-dcast(df_show, miRNA_name~+Index, mean,value.var = 'Score')
# df_show$Score<-as.numeric(df_show$Score)
# ggstatsplot::ggbetweenstats(df_show,aes(x=Site_type,y=Score))+facet_wrap(~Index)
library(ggpubr)
dodge <- position_dodge(width = 0.5)
#####################################
## Figure S19D
#####################################
ggplot(df_show,aes(x=Site_type,y=Score))+
  facet_wrap(~Index,scales = "free_y")+
  geom_violin(aes(fill=genotype),position = dodge)+
  
  geom_boxplot(aes(color=genotype),fill='white',width=0.2,position = dodge)+
  theme_bw()+
  
  stat_compare_means(aes(x=Site_type,y=Score,group=genotype),bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  scale_color_manual(values = c("black","black"))+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15),
        strip.text = element_text(size=25)
  )

###############################################################################################################
## RNA secondary structure
# BiocManager::install('RRNA')
# BiocManager::install("Structstrings")
library(Structstrings)
library(RRNA)
library(R4RNA)
library(stringr)

## import circRNA with different logFC of RBP
setwd(circRNAprofiler_dir)
# load motifs info of different genotypes and random seqs
load('./motifsFTS_circ.Rdata')
motifsFTS_circ_ref<-motifsFTS_circ
load('./motifsFTS_circ_alt.Rdata')
motifsFTS_circ_alt<-motifsFTS_circ
load('./motifsRFTS_circ_random.Rdata')
# load snp info of different snps
snp_loc_files<-list.files()[grep("_in_",list.files())]
snp_loc_file<-snp_loc_files[1]
snp_loc<-read.table(snp_loc_file,sep='\t',header = T)

## import energy difference data and mfe structure of ref vs alt
# ref
genotypes<-c("ref","alt")
structure_types<-c("mfe","therm","cent")
ref_dict<-'./RNAstructure/ref/'
alt_dict<-'./RNAstructure/alt/'
for (i in 1:length(genotypes)){
  # i=1
  for (j in 1:length(structure_types)){
    # j=2
    # ref_dict<-'./RNAstructure/ref/'
    # alt_dict<-'./RNAstructure/alt/'
    
    work_dict<-paste0('./RNAstructure/',genotypes[i],'/')
    if (j==2){
      circRNA_ref_dot<-read.table(paste0(work_dict,structure_types[j],"_structure.txt"),sep='\t')
      therm_socre=as.numeric(na.omit(unlist(lapply(circRNA_ref_dot$V1,function(x){str_match(x,pattern = "-\\d+.\\d+")}))))
      structure_data<-data.frame("i"=1,"j"=1,"length"=1,"value"=therm_socre)
    } else{
      structure_data<-readVienna(paste0(work_dict,structure_types[j],"_structure.txt"))
    }
    
    circRNA_ref_dot<-read.table(paste0(work_dict,structure_types[j],"_structure.txt"),sep='\t')
    circRNA_list<-unlist(lapply(circRNA_ref_dot$V1[grep(">",circRNA_ref_dot$V1)],function(x){strsplit(x,">")[[1]][2]}))
    score_ref_data<-structure_data[!duplicated(structure_data[,c("value")]),]
    score_ref_data$circRNA<-circRNA_list
    score_ref_data$genotype<-genotypes[i]
    score_ref_data$score_type<-structure_types[j]
    if (j==1){
      score_each_data<-score_ref_data
    } else{
      score_each_data<-rbind(score_each_data,score_ref_data)
    }
  }
  if (i==1){
    score_total_data<-score_each_data
  }else{
    score_total_data<-rbind(score_total_data,score_each_data)
  }
}
score_total_data$genotype<-factor(score_total_data$genotype,levels = c('alt',"ref"))
#####################################
## Figure S23
#####################################
## fig S23C comparison of different structure scores of circRNAs with different genotypes
head(score_total_data)
# install.packages('ggplot')
# devtools::install_github('IRkernel/IRkernel')
# install.packages("remotes")
# remotes::install_github("Rundmus/Useful4me-R_package")
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                              1))
    quantiles <- create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x","y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat ="ydensity", position ="identity", ..., draw_quantiles = NULL, trim = TRUE, scale ="area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
# fig2A-comparison of different 
ggplot(score_total_data,aes(x=score_type,y=log2(-value),fill=genotype))+
  geom_violin(alpha=0.6)+
  # geom_split_violin(alpha=0.6)+
  # geom_boxplot(notch = TRUE,width=0.5,position = "dodge",alpha=0.8)+
  geom_boxplot(notch = TRUE,width=0.5,position = position_dodge(width = 0.9),alpha=0.8)+
  geom_point(shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  # geom_line(aes(color=genotype,group=genotype),size=0.8,colour="#9C9C9C",linetype="dotted")+
  # ggpaired(score_total_data,x='score_type',y='value',fill='genotype')+
  stat_compare_means(aes(x=score_type,y=value,group=genotype),
                     bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  xlab("Score indices")+
  theme_bw()+theme(
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25),
    axis.text.x = element_text(angle=45,hjust=1,size=18),
    axis.text.y = element_text(size=18),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15),
    strip.text = element_text(size=25)
  )
# import the distribution of scores in position of each circRNA
distribution_3scores_files<-list.files(ref_dict,pattern = '_score.txt')
# distribution_3scores_files<-list.files(alt_dict,pattern = '_score.txt')
score_types<-c("mfe","pf","entropy")
prepare_distri_data<-function(distribution_3scores_file,circRNA_id){
  library(stringr)
  temp_distri_data<-read.table(distribution_3scores_file,sep='\t')
  sep_index<-c(0,which(temp_distri_data$V1=="&"))
  for (j in 1:length(sep_index)){
    # j=1
    if (j==length(sep_index)){
      temp<-temp_distri_data$V1[(sep_index[j]+1):(length(temp_distri_data$V1))]
      
    } else{
      temp<-temp_distri_data$V1[(sep_index[j]+1):(sep_index[j+1]-1)]
      
    }
    temp<-str_trim(temp,side = 'both')
    temp<-data.frame(t(data.frame(strsplit(temp," +"))))
    # temp<-data.frame()
    colnames(temp)<-c("position","value");
    temp$score_type<-score_types[j]
    if (j==1){
      temp_each_data<-temp
    } else{
      temp_each_data<-rbind(temp_each_data,temp)
    }
  }
  temp_each_data$circRNA<-circRNA_id
  return(temp_each_data)
}
convert_bins<-function(distri_data,bins=100){
  # library(pacman)
  # library(dlookr)
  library(zoo)
  ### bins!!!
  # bins=100
  distri_data<-as.numeric(distri_data)
  window_size=floor(length(distri_data)/bins)
  g<-rep(1:bins,each=window_size)
  x=tapply(distri_data[1:length(g)],g,mean)
  distri_data_bins<-data.frame(t(rbind(seq((1+window_size)/2,((length(g)+1)-window_size/2),window_size),
                                       x)))
  ### rollmean
  # window_size=length(distri_data)-bins+1
  # distri_data_bins<-data.frame(t(rbind(seq((1+window_size)/2,((length(distri_data)+1)-window_size/2),1),
  #                                (rollmean(distri_data,window_size)))))
  colnames(distri_data_bins)<-c("xaxis","mean")
  # distri_data_bins<-data.frame(t(rbind(c(distri_data_sum),xaxis)))
  return(distri_data_bins)
}
matrix_total<-list()
for (i in 1:length(distribution_3scores_files)){
  # i=3
  
  distribution_3scores_file<-paste0(alt_dict,distribution_3scores_files[i])
  circRNA_id<-strsplit(distribution_3scores_files[i],"_score.txt")[[1]][1]
  temp_each_data<-prepare_distri_data(distribution_3scores_file,circRNA_id)
  for (j in 1:length(score_types)){
    # j=1
    temp_each_bins<-convert_bins(temp_each_data$value[temp_each_data$score_type==score_types[j]],bins = 100)
    # dim(temp_each_bins)
    temp_each_bins$score_type<-score_types[j]
    if (j ==1){
      temp_each_bins_total<-temp_each_bins
    } else{
      temp_each_bins_total<-rbind(temp_each_bins_total,temp_each_bins)
    }
  }
  temp_each_bins_total$circRNA<-circRNA_id
  
  for (j in 1:length(score_types)){
    if (i==1){
      matrix_total[[j]]<-t(data.frame(temp_each_bins_total$mean[temp_each_bins_total$score_type==score_types[j]]))
      rownames(matrix_total[[j]])[i]<-circRNA_id
    }else{
      matrix_total[[j]]<-rbind(matrix_total[[j]],temp_each_bins_total$mean[temp_each_bins_total$score_type==score_types[j]])
      rownames(matrix_total[[j]])[i]<-circRNA_id
    }
  }

}
#####################################
## Figure S23A
#####################################
ref_matrix_total<-matrix_total
alt_matrix_total<-matrix_total
matrix_total<-list(ref_matrix_total,alt_matrix_total)
## ref vs alt matrix
library(pheatmap)
bins<-100
diff_matrix_total<-list()
for (i in 1:length(ref_matrix_total)){
  # i=2
  diff_matrix_total[[i]]<-alt_matrix_total[[i]]-ref_matrix_total[[i]]
  colnames(ref_matrix_total[[i]])<-paste0("ref",seq(1,bins));colnames(alt_matrix_total[[i]])<-paste0("alt",seq(1,bins));
  annotation_col<-rbind(data.frame('pos'=paste0("ref",seq(1,bins)),'genotype'='ref'), 
                        data.frame('pos'=paste0("alt",seq(1,bins)),'genotype'='alt'))
  rownames(annotation_col)<-annotation_col$pos;annotation_col<-annotation_col[-1]
  # annotation_col<-data.frame("genotype"=)
  combine_matrix_total<-cbind(ref_matrix_total[[i]],alt_matrix_total[[i]])
  tiff(paste0(score_types[i],"ref vs alt.tiff"),width = 2000,height = 1000,res = 150)
  pheatmap(combine_matrix_total,annotation_col = annotation_col,scale = 'column',cluster_cols = F,
           number_format= "%.2f",
           color = colorRampPalette(colors = c("blue","white","red"))(100))
  dev.off()
  
  tiff(paste0(score_types[i],"diff_matrix.tiff"),width = 1500,height = 1000,res = 150)
  pheatmap(diff_matrix_total[[i]],scale = 'column',cluster_cols = F,
           color = colorRampPalette(colors = c("blue","white","red"))(100))
  dev.off()
  
  
}
#####################################
## Figure S23B
#####################################
## PCR correlatino of different score types
library(corrplot)
for (i in 1:(length(diff_matrix_total)-1)){
    for (j in (i+1):length(diff_matrix_total)){
      print(paste0(i,"vs",j))
      # i=1
      # j=3
      res=cor(t(diff_matrix_total[[i]]),t(diff_matrix_total[[j]]),method='pearson')
      res[is.na(res)]<-0
      corr_res<-diag(res)
      df<-data.frame(x=rowMeans(diff_matrix_total[[i]]),y=rowMeans(diff_matrix_total[[j]]),corr=corr_res)
      df$corr_type[sign(df$corr)==1]<-"Positive";df$corr_type[sign(df$corr)==-1]<-"Negative"
      df$corr_type[sign(df$corr)==0]<-"No correlation"
      tiff(paste0("fig2C-correlation pcr of diff matrix between ref and alt of ",score_types[i]," vs ",score_types[j],".tiff"),
          width = 2000,height = 1500,res=300)
      p<-ggplot(df,aes(x=x,y=y))+
        geom_point(aes(size=abs(corr),color=corr_type),shape=20)+
        # geom_smooth(method=lm)+
        xlab(score_types[i])+ylab(score_types[j])+
        scale_color_manual(values = pal_aaas("default")(10)[c(1,9,2)])+
        guides(color=guide_legend(title = "Correlation type",override.aes = list(size=5)))+
        guides(size=guide_legend(title = "PCC"))+
        scale_size_continuous(range = c(0,max(abs(df$corr)*10)))+
        theme_bw()+
        theme(axis.title.x = element_text(size=25),
              axis.title.y = element_text(size=25),
              axis.text = element_text(size=18),
              legend.title = element_text(size=18),
              legend.text = element_text(size=15))
      print(p)
      dev.off()
      # cor(t(diff_matrix_total[[i]][1,]),t(diff_matrix_total[[j]][1,]),method='spearman')
      # pheatmap(res,cluster_cols = F,cluster_rows = F)
      # corrplot(res,method = 'circle',type='lower')
    }
}

for (i in 1:length(diff_matrix_total)){
  temp<-data.frame(t(data.frame(colMeans(diff_matrix_total[[i]]))))
  temp_score_ref<-data.frame(t(data.frame(colMeans(ref_matrix_total[[i]]))))
  temp_score_alt<-data.frame(t(data.frame(colMeans(alt_matrix_total[[i]]))))
  temp_score_ref$score<-score_types[i];temp_score_ref$genotype<-"ref"
  temp_score_alt$score<-score_types[i];temp_score_alt$genotype<-"alt"
  colnames(temp_score_ref)[1:bins]<-seq(1,bins,1);colnames(temp_score_alt)[1:bins]<-seq(1,bins,1)
  temp_score_each<-rbind(temp_score_ref,temp_score_alt)
  if (i==1){
    temp_score_total<-temp_score_each
  }else{
    temp_score_total<-rbind(temp_score_total,temp_score_each)
  }
  
  temp$score_type<-score_types[i]
  if (i==1){
    # i
    diff_distribution_total<-temp
    # diff_distribution_total$score_type[i]<-score_types[i]
  }else{
    
    diff_distribution_total<-rbind(diff_distribution_total,temp)
    # diff_distribution_total$score_type[i]<-score_types[i]
  }
  
}
#####################################
## Figure S23D
#####################################
## ref vs alt distribution
head(temp_score_total)
df_show<-melt(temp_score_total,id.vars = c("score","genotype"),variable.name = "Bins",value.name = "mean")
# df_show$Bins<-as.numeric(unlist(lapply(df_show$Bins,function(x){strsplit(as.character(x),"alt")[[1]][2]})))
df_show$Bins<-as.numeric(df_show$Bins)
diff<-abs(df_show$mean[df_show$score=="mfe"&df_show$genotype=="alt"]-df_show$mean[df_show$score=="mfe"&df_show$genotype=="ref"])
diff_cent<-abs(df_show$mean[df_show$score=="pf"&df_show$genotype=="alt"]-df_show$mean[df_show$score=="pf"&df_show$genotype=="ref"])
index<-intersect(which(diff>=0.3),which(diff_cent>=0.3))
split_index<-split(index,cummax(c(1,diff(index))))
ggplot(df_show,aes(x=Bins,y=mean,color=score))+geom_line(aes(linetype=genotype),size=1.5)+
  scale_y_continuous(expand = c(0,0))+
  # annotate("rect", x=index,ymin = 0, ymax = 80,alpha = .5, fill = "lightblue")+
  annotate("rect", xmin = min(unlist(split_index[2])) , xmax = max(unlist(split_index[2])), ymin = 0, ymax = 80, alpha = .5, fill = "lightblue")+
  annotate("rect", xmin = min(unlist(split_index[3])) , xmax = max(unlist(split_index[3])), ymin = 0, ymax = 80, alpha = .5, fill = "lightblue")+
  scale_color_manual(values = c(pal_npg("nrc")(10)[6:10],pal_aaas("default")(12)))+
  theme_bw()+ylab('Mean value')+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))
#####################################
## Figure S23E
#####################################
## differential matrix distribution
df_show<-melt(diff_distribution_total,id.vars = "score_type",variable.name = "Bins",value.name = "mean")
df_show$Bins<-as.numeric(unlist(lapply(df_show$Bins,function(x){strsplit(as.character(x),"X")[[1]][2]})))
ggplot(df_show,aes(x=Bins,y=mean,color=score_type))+geom_line(size=1.5)+
  scale_y_continuous(expand = c(0,0))+
  annotate("rect", xmin = min(unlist(split_index[2])) , xmax = max(unlist(split_index[2])), ymin = -2, ymax = 2.5, alpha = .5, fill = "lightblue")+
  annotate("rect", xmin = min(unlist(split_index[3])) , xmax = max(unlist(split_index[3])), ymin = -2, ymax = 2.5, alpha = .5, fill = "lightblue")+
  ylab('Mean of difference between ref and alt')+
  scale_color_manual(values = c(pal_npg("nrc")(10)[6:10],pal_aaas("default")(12)))+
  theme_bw()+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

###############################################################################################################
#####################################
## Figure 8
#####################################
### comparison of internal structure (dsRNA internal_loop<=4bp,total length(without indel)>=16bp)
# get snp loc function
get_snp_loc<-function(snp_loc,circRNA_id){
  snp_loc_narm<-na.omit(snp_loc[snp_loc$id==circRNA_id,])
  if (sum(!is.na(snp_loc[snp_loc$id==circRNA_id,]))>=2){
    snp_names<-names(snp_loc_narm)[2:length(snp_loc_narm)]
    snp_loc_narm<-as.numeric(snp_loc_narm[-1]);names(snp_loc_narm)<-snp_names
    return(snp_loc_narm)
  }else{
    # snp_loc_narm<-NA
    return(NA)
  }

  
}
# get motif loc and related annotation
get_motif_loc<-function(motifsFTS_circ,circRNA_id,count_cutoff=2){
  locations<-motifsFTS_circ$circ$locations
  motif_data<-motifsFTS_circ$circ$motifs
  motif_loc_data<-locations[locations$id==circRNA_id,];motif_loc_data<-motif_loc_data[which(motif_loc_data!="NA")]
  motif_loc_data<-melt(motif_loc_data,id.vars = "id",variable.name = "motif",value.name = "loc")
  motif_loc_data$count<-unlist(lapply(motif_loc_data$loc,function(x){length(strsplit(x,",",fixed=T)[[1]])}))
  motif_loc_data_filt<-motif_loc_data[motif_loc_data$count>=count_cutoff,]
  colnames(motif_loc_data_filt)[1]<-"circRNA_id"
  motif_loc_data_filt<-merge(motif_loc_data_filt,motif_data,by='motif')
  return(motif_loc_data_filt)
}
head(motifsFTS_circ$circ$locations)
genotypes<-c("ref","alt")
structure_types<-c("mfe","cent")
ref_dict<-'./RNAstructure/ref/'
alt_dict<-'./RNAstructure/alt/'

for (i in 1:length(genotypes)){
  # i=1
  print(paste0("i:",i))
  for (j in 1:length(structure_types)){
    print(paste0("j:",j))

    work_dict<-paste0('./RNAstructure/',genotypes[i],'/')
    if (structure_types[j]=="therm"){
      circRNA_ref_dot<-read.table(paste0(ref_dict,structure_types[j],"_structure.txt"),sep='\t')
      therm_socre=as.numeric(na.omit(unlist(lapply(circRNA_ref_dot$V1,function(x){str_match(x,pattern = "-\\d+.\\d+")}))))
      structure_data<-data.frame("i"=1,"j"=1,"length"=1,"value"=therm_socre)
    } else{
      structure_data<-readVienna(paste0(ref_dict,structure_types[j],"_structure.txt"))
    }
    idx<-which(duplicated(structure_data$value)==F)
    idx<-c(idx,dim(structure_data)[1]+1)
    circRNA_ref_dot<-read.table(paste0(ref_dict,structure_types[j],"_structure.txt"),sep='\t')
    circRNA_list<-unlist(lapply(circRNA_ref_dot$V1[grep(">",circRNA_ref_dot$V1)],function(x){strsplit(x,">")[[1]][2]}))
    for (k in 1:length(circRNA_list)){
      # k<-3
      print(paste0("k:",k))
      circRNA_id<-circRNA_list[k]
      # import pp data (paired probability)
      pp_each<-read.table(paste0(work_dict,str_replace_all(circRNA_id,":","_"),"_pp.txt"),sep=' ')
      colnames(pp_each)<-c('i',"j","pp","box")
      pp_each<-pp_each[,c("i","j","pp")]%>%group_by(i,j)%>%mutate(pp=mean(pp))
      pp_each<-pp_each[!duplicated(pp_each),]
      
      total_length<-motifsFTS_circ$circ$targets$length[motifsFTS_circ$circ$targets$id==circRNA_id]
      get_snp_loc(snp_loc,circRNA_id)
      data<-structure_data[(idx[k]):(idx[k+1]-1),]
      data$bins<-round(data$i/ceiling(total_length/bins))+1
      each_data<-matrix_total[[i]][[j]]
      each_data_value<-each_data[rownames(each_data)==str_replace_all(circRNA_id,":","_")]
      data$value<-each_data_value[data$bins]

      res=identify_subelements(data,total_length,snp_loc,rbp_loc = NULL,loop_len =5)
      duplex_res<-res[[1]];loop_res<-res[[2]]
      res_filt<-duplex_res[!is.na(duplex_res$dsRNA),];loop_res_filt<-loop_res[!is.na(loop_res$isInduplex),]
      if (length(res_filt$i)==0){
        res_filt[1,]<-c(NA,NA,0,NA,0,NA,NA,NA)
        
        res_filt$circRNA_id<-circRNA_id
        res_filt$duplex_number<-NA
        res_filt$pp<-NA
        # res_filt$paired_length<-NA
        
      } else{
        res_filt$circRNA_id<-circRNA_id
        res_filt$duplex_number<-max(unlist(lapply(res_filt$dsRNA,
                                                  function(x){as.numeric(strsplit(x,"duplex",fixed=T)[[1]][2])})))
        # colnames(pp_each)[2]<-'i'
        # if (structure_types[j]=="cent"){
        #   
        # }
        temp<-merge(res_filt,pp_each[,c('i',"pp")],by=c("i"))%>%group_by(i,j)%>%mutate(pp=mean(pp))
        
        res_filt<-temp[!duplicated(temp),]
        # res_filt$paired_length<-
      }
      if (length(loop_res_filt$isInduplex)==0){
        loop_res_filt[1,]<-c(NA,NA,NA)
        
        loop_res_filt$circRNA_id<-circRNA_id
        loop_res_filt$loop_number<-0
        
      } else{
        loop_res_filt$circRNA_id<-circRNA_id
        for (l in 1:length(unique(loop_res_filt$isInduplex))){
          # l=1
          index<-which(loop_res_filt$isInduplex==unique(loop_res_filt$isInduplex)[l])
          loop_index<-intersect(grep("inner-loop",loop_res_filt$loop_type),index)
          loop_res_filt[loop_index,"loop_number"]<-max(unlist(lapply(
          loop_res_filt$loop_type[loop_index],function(x){as.numeric(strsplit(x,"inner-loop",fixed=T)[[1]][2])})))
          bulge_index<-intersect(grep("bulge",loop_res_filt$loop_type),index)
          loop_res_filt[bulge_index,"loop_number"]<-max(unlist(lapply(
          loop_res_filt$loop_type[bulge_index],
          function(x){as.numeric(strsplit(x,"bulge",fixed=T)[[1]][2])})))
        }
        
      }
      res_filt$genotype<-genotypes[i]
      res_filt$score_type<-structure_types[j]
      loop_res_filt$genotype<-genotypes[i]
      loop_res_filt$score_type<-structure_types[j]
      if (k==1) {
        duplex_each<-res_filt
        loop_each<-loop_res_filt
      } else {
        res_filt$dsRNA<-as.character(res_filt$dsRNA)
        duplex_each<-rbind(duplex_each,res_filt)
        loop_each<-rbind(loop_each,loop_res_filt)
      }
    }
    
   
    if (j==1){
      duplex_each_data<-duplex_each
      loop_each_data<-loop_each
    } else{
      duplex_each_data<-rbind(duplex_each_data,duplex_each)
      loop_each_data<-rbind(loop_each_data,loop_each)
    }
  }
  if (i==1){
    duplex_total_data<-duplex_each_data
    loop_total_data<-loop_each_data
  }else{
    duplex_total_data<-rbind(duplex_total_data,duplex_each_data)
    loop_total_data<-rbind(loop_total_data,loop_each_data)
  }
}
# identify_loops<-function(data){
#   
# }
identify_subelements<-function(data,total_length,snp_loc,rbp_loc,loop_len=5){

  df<-data
  df_loop<-data.frame("loc"=setdiff(seq(1,total_length),c(df$i,df$j)))
  total_loc<-seq(1,total_length,1)
  res<-data.frame("location"=total_loc)
  duplex_id<-1;sep<-0;left_loc<-df$i;right_loc<-df$j;paired_loc<-c(left_loc,right_loc)
  # loop_len<-5;
  # length(df$i)
  for (i in 1:length(df$i)){
    if (i==1){
      duplex_id<-1;sep<-0;locs_index<-list(i);left_sep_list<-list();right_sep_list<-list();sep_list<-c(0)
      bulge_id<-1;inner_loop_id<-1;
      next
    }
    
    if ((df$j[i]<df$j[i-1]) & (max(sep_list)<=loop_len)) {
      left_sep<-df$i[i]-df$i[i-1]-1
      right_sep<-df$j[i-1]-df$j[i]-1
      if ((left_sep+right_sep>0)&left_sep!=right_sep){
          unpaired_left_loc<-seq(df$i[i-1],df$i[i]-1);unpaired_left_loc<-unpaired_left_loc[-1]
          unpaired_right_loc<-seq(df$j[i],df$j[i-1]-1);unpaired_right_loc<-unpaired_right_loc[-1]
          unpaired_loc<-c(unpaired_left_loc,unpaired_right_loc)
          idx<-unlist(lapply(unpaired_loc, function(x){which(df_loop$loc==x)}))
          df_loop[idx,"loop_type"]=paste0("bulge",bulge_id)
        bulge_id<-bulge_id+1
      }else if ((left_sep+right_sep>0)&left_sep==right_sep){
        unpaired_left_loc<-seq(df$i[i-1],df$i[i]-1);unpaired_left_loc<-unpaired_left_loc[-1]
        unpaired_right_loc<-seq(df$j[i],df$j[i-1]-1);unpaired_right_loc<-unpaired_right_loc[-1]
        unpaired_loc<-c(unpaired_left_loc,unpaired_right_loc)
        idx<-unlist(lapply(unpaired_loc, function(x){which(df_loop$loc==x)}))
        df_loop[idx,"loop_type"]=paste0("inner-loop",inner_loop_id)
        inner_loop_id<-inner_loop_id+1
      }
      left_sep_list<-append(left_sep_list,left_sep)
      right_sep_list<-append(right_sep_list,right_sep)
      locs_index<-append(locs_index,i)
      last<-locs_index[length(locs_index)][[1]];first<-locs_index[1][[1]]
      ds_RNA_length<-min(df$i[last]-df$i[first]+1,df$j[first]-df$j[last]+1)
      sep_list<-c(unlist(left_sep_list),unlist(right_sep_list))
      if (max(sep_list)>loop_len){
        last<-locs_index[length(locs_index)][[1]]-1;first<-locs_index[1][[1]]
      } else{
        last<-locs_index[length(locs_index)][[1]];first<-locs_index[1][[1]]
      }
      
    } else if (df$j[i]>df$j[i-1]) {
      
      ds_RNA_length<-min(df$i[last]-df$i[first]+1,df$j[first]-df$j[last]+1)
      df[first:last,"ds_RNA_length"]<-ds_RNA_length
      df[first:last,"paired_length"]<-last-first+1
      # sep_list<-c(unlist(left_sep_list),unlist(right_sep_list))
      idx<-(df_loop$loc<df$i[last]&df_loop$loc>df$i[first])|
        (df_loop$loc<df$j[first]&df_loop$loc>df$j[last])
      if ((max(sep_list)<=loop_len)&(ds_RNA_length>=16)){
        if (sum(idx)>0){
          # idx<-1
          df_loop$isInduplex[idx]<-paste0("duplex",duplex_id)
        }
        
        df[first:last,"dsRNA"]<-paste0("duplex",duplex_id)
        duplex_id<-duplex_id+1
        
      } else {
        df[first:last,"dsRNA"]<-NA
          df_loop$isInduplex[idx]<-NA
      }
     
      ## reset
      locs_index<-list();locs_index<-append(locs_index,i);left_sep_list<-list();right_sep_list<-list();sep_list<-c(0)
      bulge_id<-1;inner_loop_id<-1
      next
    }
  
  }
  res<-list(df,df_loop)
  return(res)
}

#####################################
## Figure 8A
#####################################
# comparison of duplex number
df_show<-duplex_total_data[(!duplicated(duplex_total_data[,c('value','genotype','score_type')]))&!is.na(duplex_total_data$value),]
df_show$combine<-factor(interaction(df_show$score_type,df_show$genotype),levels = c("cent.ref","cent.alt",
                                                                                    "mfe.ref","mfe.alt"))
df_show$group<-interaction(df_show$circRNA_id,df_show$score_type)
ggplot(df_show,aes(x=combine,y=ds_RNA_length,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.5)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10]))+
  # ggpaired(df_show,x='combine',y='duplex_number',fill='combine')+
  # stat_compare_means(aes(x=score_type,y=duplex_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  # geom_label()+
  annotate("text", x=1.5, y=50, label="p=0.96",size=5) +
  annotate("text", x=3.5, y=50, label="p=0.95",size=5) +
  ylab("dsRNA length (bp)")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))
  
head(df_show)

#####################################
## Figure 8C
#####################################
## comparison of pp between ref and alt
ggplot(df_show,aes(x=combine,y=paired_length,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.5)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10]))+
  # ggpaired(df_show,x='combine',y='duplex_number',fill='combine')+
  # stat_compare_means(aes(x=score_type,y=paired_length,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  # geom_label()+
  annotate("text", x=1.5, y=50, label="p=0.92",size=5) +
  annotate("text", x=3.5, y=50, label="p=1.00",size=5) +
  ylab("paired length (bp)")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

#####################################
## Figure 8D
#####################################
ggplot(df_show,aes(x=combine,y=pp,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.5)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10]))+
  # ggpaired(df_show,x='combine',y='duplex_number',fill='combine')+
  # stat_compare_means(aes(x=score_type,y=pp,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  # geom_label()+
  annotate("text", x=1.5, y=1.05, label="p=0.017",size=5) +
  annotate("text", x=3.5, y=1.05, label="p=0.026",size=5) +
  ylab("paired probability")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))
# df_show$group<-

#####################################
## Figure 8B
#####################################
ggplot(df_show,aes(x=combine,y=duplex_number,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.4)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  # stat_compare_means(aes(x=score_type,y=duplex_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  annotate("text", x=1.5, y=9.5, label="p=0.88",size=5) +
  annotate("text", x=3.5, y=9.5, label="p=1",size=5) +
  ylab("duplex number")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

###############################################################################################################
## comparison of ref vs alt within two specific regions
specific_bins<-c(split_index[[3]],split_index[[2]])
index<-unlist(lapply(specific_bins,function(x){which(df_show$bins==x)}))
df_show_specific_region<-df_show[index,]


ggplot(df_show_specific_region,aes(x=combine,y=ds_RNA_length,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.5)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10]))+
  # ggpaired(df_show,x='combine',y='duplex_number',fill='combine')+
  # stat_compare_means(aes(x=score_type,y=duplex_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  # geom_label()+
  annotate("text", x=1.5, y=50, label="p=0.85",size=5) +
  annotate("text", x=3.5, y=50, label="p=1.00",size=5) +
  ylab("dsRNA length (bp)")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

head(df_show)

## comparison of pp between ref and alt
ggplot(df_show_specific_region,aes(x=combine,y=paired_length,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.5)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10]))+
  # ggpaired(df_show,x='combine',y='duplex_number',fill='combine')+
  # stat_compare_means(aes(x=score_type,y=paired_length,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  # geom_label()+
  annotate("text", x=1.5, y=37, label="p=0.99",size=5) +
  annotate("text", x=3.5, y=37, label="p=1.00",size=5) +
  ylab("paired length (bp)")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

ggplot(df_show_specific_region,aes(x=combine,y=pp,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.5)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10]))+
  # ggpaired(df_show,x='combine',y='duplex_number',fill='combine')+
  # stat_compare_means(aes(x=score_type,y=pp,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  # geom_label()+
  annotate("text", x=1.5, y=1.05, label="p=0.29",size=5) +
  annotate("text", x=3.5, y=1.05, label="p=0.21",size=5) +
  ylab("paired probability")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))
# df_show$group<-

ggplot(df_show_specific_region,aes(x=combine,y=duplex_number,fill=genotype))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.4)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2))+
  geom_line(aes(group=group),size=0.8,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  # stat_compare_means(aes(x=score_type,y=duplex_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  annotate("text", x=1.5, y=9.5, label="p=0.85",size=5) +
  annotate("text", x=3.5, y=9.5, label="p=1.00",size=5) +
  ylab("duplex number")+xlab("Score type & genotype")+
  scale_x_discrete(labels = c("centroid ref","centroid alt","mfe ref","mfe alt"))+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

#####################################
## Figure 8E
#####################################
## comparison of loop data at whole circRNA level
# split_index
# at whole circRNA level
df_show<-loop_total_data[(!duplicated(loop_total_data[,c(
                           'loop_number','genotype','score_type')]))&
                           !is.na(loop_total_data$loc),]
df_show$combine<-factor(interaction(df_show$score_type,df_show$genotype),levels = c("cent.ref","cent.alt",
                                                                                    "mfe.ref","mfe.alt"))
df_show$loop_type<-str_sub(df_show$loop_type,1,str_count(df_show$loop_type)-1)
df_show$group<-interaction(df_show$circRNA_id,df_show$score_type)
ggplot(df_show,aes(x=combine,y=loop_number,fill=genotype))+
  facet_wrap(vars(loop_type))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.4)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2,height = 0))+
  geom_line(aes(group=group),size=1,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  # stat_compare_means(aes(x=score_type,y=loop_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  ylab("subelements number")+xlab("Score type & genotype")+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        strip.text = element_text(size=25),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))


#####################################
## Figure 8F
#####################################
## comparison of loop data at each duplex level
df_show<-loop_total_data[(!duplicated(loop_total_data[,c('isInduplex',
                                                         'loop_number','genotype','score_type')]))&
                           !is.na(loop_total_data$loc),]
df_show$combine<-factor(interaction(df_show$score_type,df_show$genotype),levels = c("cent.ref","cent.alt",
                                                                                    "mfe.ref","mfe.alt"))
df_show$loop_type<-str_sub(df_show$loop_type,1,str_count(df_show$loop_type)-1)
df_show$group<-interaction(df_show$circRNA_id,df_show$score_type)
df_show_matrix<-dcast(df_show,circRNA_id+combine+loop_type~isInduplex+loop_number,sum,value.var = "loop_number")
rownames(df_show_matrix)<-interaction(df_show_matrix$circRNA_id,df_show_matrix$combine,df_show_matrix$loop_type)
annotation_col<-df_show_matrix[,c('circRNA_id','combine','loop_type')]
df_show_matrix<-df_show_matrix[order(df_show_matrix$combine),]
df_show_matrix<-df_show_matrix[order(df_show_matrix$loop_type),]
rownames(annotation_col)<-interaction(annotation_col$circRNA_id,annotation_col$combine,annotation_col$loop_type)
annotation_col<-annotation_col[-1]
annotation_col$genotype<-unlist(lapply(annotation_col$combine,function(x){strsplit(as.character(x),".",fixed=T)[[1]][2]}))
annotation_col$score_type<-unlist(lapply(annotation_col$combine,function(x){strsplit(as.character(x),".",fixed=T)[[1]][1]}))
annotation_col<-annotation_col[-1]
pheatmap(df_show_matrix[,4:dim(df_show_matrix)[2]],cluster_rows = F,cluster_cols = T,
         show_rownames = F,border_color = "NA",annotation_row = annotation_col)
# df_show%>%group_by(circRNA_id,combine,loop_type)%>%mutate()
ggplot(df_show,aes(x=combine,y=loop_number,fill=genotype))+
  facet_wrap(vars(loop_type))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.4)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2,height = 0))+
  # geom_jitter(width = 0.4,height = 0,size=4)+
  geom_line(aes(group=group),size=1,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  
  # stat_compare_means(aes(x=score_type,y=loop_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  
  ylab("subelements number")+xlab("Score type & genotype")+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        strip.text = element_text(size=25),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

head(df_show)

####################################################################################################
### Loop number within specific region inside circRNA
# specific_bins<-c(split_index[[3]],split_index[[2]])
circRNA_len_data<-motifsFTS_circ$circ$targets[,c("id","length")]
colnames(circRNA_len_data)[1]<-"circRNA_id"
df_show<-merge(df_show,circRNA_len_data,by='circRNA_id')
df_show$bins<-round(df_show$loc/ceiling(df_show$length/bins))+1
index<-unlist(lapply(specific_bins,function(x){which(df_show$bins==x)}))
df_show_specific_region<-df_show[index,]
## specific region

ggplot(df_show_specific_region,aes(x=combine,y=loop_number,fill=genotype))+
  # facet_wrap(vars(loop_type))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.4)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2,height = 0))+
  geom_line(aes(group=group),size=1,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  # stat_compare_means(aes(x=score_type,y=loop_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  ylab("subelements number")+xlab("Score type & genotype")+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        strip.text = element_text(size=25),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

## comparison of loop data at each duplex level
df_show_matrix<-dcast(df_show_specific_region,circRNA_id+combine+loop_type~isInduplex+loop_number,sum,value.var = "loop_number")
rownames(df_show_matrix)<-interaction(df_show_matrix$circRNA_id,df_show_matrix$combine,df_show_matrix$loop_type)
annotation_col<-df_show_matrix[,c('circRNA_id','combine','loop_type')]
df_show_matrix<-df_show_matrix[order(df_show_matrix$combine),]
df_show_matrix<-df_show_matrix[order(df_show_matrix$loop_type),]
rownames(annotation_col)<-interaction(annotation_col$circRNA_id,annotation_col$combine,annotation_col$loop_type)
annotation_col<-annotation_col[-1]
annotation_col$genotype<-unlist(lapply(annotation_col$combine,function(x){strsplit(as.character(x),".",fixed=T)[[1]][2]}))
annotation_col$score_type<-unlist(lapply(annotation_col$combine,function(x){strsplit(as.character(x),".",fixed=T)[[1]][1]}))
annotation_col<-annotation_col[-1]
pheatmap(df_show_matrix[,4:dim(df_show_matrix)[2]],cluster_rows = F,cluster_cols = T,
         show_rownames = F,border_color = "NA",annotation_row = annotation_col)
# df_show%>%group_by(circRNA_id,combine,loop_type)%>%mutate()
ggplot(df_show_specific_region,aes(x=combine,y=loop_number,fill=genotype))+
  facet_wrap(vars(loop_type))+
  geom_violin(alpha=0.6)+
  geom_boxplot(alpha=0.8,width=0.4)+
  theme_bw()+
  scale_fill_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  geom_point(aes(colour=genotype,shape=score_type),size=2,position = position_jitter(0.2,height = 0))+
  # geom_jitter(width = 0.4,height = 0,size=4)+
  geom_line(aes(group=group),size=1,colour="black",linetype="dotted",alpha=0.6)+
  scale_color_manual(values = c(pal_npg("nrc")(10)[5:10],pal_aaas("default")(12)))+
  
  # stat_compare_means(aes(x=score_type,y=loop_number,group=genotype),
  # bracket.size='20',method='wilcox.test',size = 5,label = 'p.format')+
  
  ylab("subelements number")+xlab("Score type & genotype")+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        strip.text = element_text(size=25),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))
