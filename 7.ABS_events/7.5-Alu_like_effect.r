#!/bin/Rscript

library(ggplot2)
library(stringr)
library(parallel)


############################################################################################################
working_dir <- "path_to_your_working_dir"
all_RNA_seq_file='./SRR_list.txt'
output_dir <- './output/Alu_iden'
cores <- 6
setwd(working_dir)
dir.create(output_dir)
cl<-makeCluster(cores)


# blast_output_total
# blast_res<-read.table('./output/Alu_iden/blast_output.txt',sep='\t')
blast_res<-read.table('./ABS_event/Alu_iden/blast_output_total.txt',sep='\n')

filt_idx<-str_match(unlist(blast_res),'NA')
blast_res_filt<-data.frame(t(data.frame(strsplit(unlist(blast_res),' ',fixed = T))))
# blast_res_filt<-blast_res_filt[,c(1:3,13)]
colnames(blast_res_filt)<-c('circ','Up','Down','Bitscore')

blast_res_filt<-blast_res_filt[blast_res_filt$Bitscore!='NA',]
# extract loc of alu-like and loc
extract_loc<-function(x,id,flank){
  if (id == 'alu') {
    tmp<-strsplit(x,'::',fixed = T)[[1]][2]
    if (flank == 'up') {
      
      out<-strsplit(tmp,'-',fixed = T)[[1]][2]
    } else if (flank=='down') {
      start_tmp<-strsplit(tmp,':',fixed = T)[[1]][2]
      out<-strsplit(start_tmp,'-',fixed = T)[[1]][1]
    } else if (flank =='element') {
      # element
      out<-strsplit(x,'::',fixed = T)[[1]][1]
    }
    
  } else if (id == 'circ') {
    tmp<-strsplit(x,':',fixed = T)[[1]][2]
    if (flank == 'up') {
      out<-strsplit(tmp,'|',fixed = T)[[1]][1]
      
    } else {
      
      out<-strsplit(tmp,'|',fixed = T)[[1]][2]
    }
    
  }
  if (flank=='element') {
    return(out)
  } else {
    return(as.numeric(out))
    
  }
}
# 
blast_res_filt$upL<-unlist(lapply(blast_res_filt$circ,extract_loc,'circ','up'))-unlist(lapply(blast_res_filt$Up,extract_loc,'alu','up'))
blast_res_filt$downL<-unlist(lapply(blast_res_filt$Down,extract_loc,'alu','down'))-unlist(lapply(blast_res_filt$circ,extract_loc,'circ','down'))

blast_res_filt$sL<-apply(blast_res_filt[,c('upL','downL')],1,min);blast_res_filt$lL<-apply(blast_res_filt[,c('upL','downL')],1,max)
blast_res_filt$L<-rowSums(blast_res_filt[,c('upL','downL')])
blast_res_filt$CSI<-blast_res_filt$sL*(as.numeric(blast_res_filt$Bitscore))**2/(blast_res_filt$lL*(blast_res_filt$L)**2)
blast_res_filt<-blast_res_filt[blast_res_filt$CSI>0,]

# Distribution of Alu-like elements
blast_res_filt$pair_l<-unlist(lapply(blast_res_filt$Up,extract_loc,'alu','element'));blast_res_filt$pair_r<-unlist(lapply(blast_res_filt$Down,extract_loc,'alu','element'))
blast_res_filt$pair<-paste(blast_res_filt$pair_l,blast_res_filt$pair_r,sep='-')
blast_res_filt$count<-1
# write.table(blast_res_filt,'./ABS_event/Alu_iden/blast_res_filt.txt',quote = F,row.names = F,sep='\t')

blast_res_filt<-read.table('./output/Alu_iden/blast_res_filt.txt',header=T,sep='\t')
df<-blast_res_filt[,c('pair','count','CSI')]
df$label<-'Others';df$label[!is.na(str_match(df$pair,'Bov'))]<-'Alu_like';df$label[!is.na(str_match(df$pair,'BOV'))]<-'Alu_like'
df$label<-factor(df$label)
pairs<-unique(df$pair);df_show<-data.frame('type'=pairs)
for (i in 1:length(pairs)){
  df_show$freq[i]<-sum(df$count[df$pair==pairs[i]])/sum(df$count)

}
df_show$label<-paste0(df_show$type,'(',round(df_show$freq,3)*100,'%)')
df_show$label[df_show$freq<=0.01]<-paste0('Others','(',round(sum(df_show$freq[df_show$freq<=0.01]),3),'%)')

#####################################
## Fig 5
#####################################
library(ggsci)
## Fig 5A
ggplot(df_show,aes(x='',y=freq,fill=label))+geom_bar(stat = 'identity')+coord_polar(theta = 'y')+
  guides(fill=guide_legend(nrow = 4,byrow = T,title = 'Paired type'))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(face = 'italic',size = 12),
        legend.title = element_text(face = 'italic',size = 15),
        legend.position = 'bottom',
        # legend.
        axis.title = element_blank(),axis.ticks = element_blank())+
  
  scale_fill_manual(values=pal_igv()(23))
## Fig 5B
df$log2CSI<-log2(df$CSI)
ggplot(df)+geom_boxplot(aes(x=label,y=log2(CSI),fill=label))+theme_bw()+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))
## Fig 5C
# ggplot(df[df$label=='Alu_like',])+geom_boxplot(aes(x=pair,y=log2(CSI),fill=pair))+theme_bw()+
  # theme(axis.title.x = element_text(size=18),
        # axis.title.y = element_text(size=18),
        # axis.text.x = element_text(angle=45,hjust=1,size=18),
        # axis.text.y = element_text(size=15),
        # legend.title = element_text(size=15),
        # legend.text = element_text(size=12))
  
# hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(23)
# ## Fig 5C
# ggstatsplot::ggbetweenstats(
  # data = df,
  # x=label,y=log2CSI,
  # type='p',
  # package = 'ggsci',
  # palette = 'alternating_igv',
  # centrality.plotting=F,
  # messages = F
# )+guides(color=guide_legend(title = "Paired type"))  +
  # theme_bw()+
  # xlab('Alu-like pairs')+theme(axis.title.x = element_text(size=25),
                               # axis.title.y = element_text(size=25),
                               # axis.text.x = element_text(size=18),
                               # axis.text.y = element_text(size=18),
                               # legend.title = element_text(size=18),
                               # legend.text = element_text(size=15))

## Fig 5C
df$shape<-1
shapes<-seq(1,23)
names(shapes)<-unique(df$pair[df$label=='Alu_like'])
for (i in 1:length(shapes)){
  df$shape[df$label=='Alu_like'&df$pair==names(shapes[i])]<-shapes[i]
}


ggstatsplot::ggbetweenstats(
  data = df[df$label=='Alu_like',],
  x=pair,y=log2CSI,
  type='p',
  # point.args = list(shape=1),
  package = "ggsci", # package from which color palette is to be taken
  palette = "default_igv",
  centrality.plotting=F,
  messages = F
)+
  # scale_colour_manual(values = c(pal_npg('nrc')(10),pal_aaas('default')(10)))+
  theme_bw()+
  # scale_
  guides(color=guide_legend(title = "Paired type"))  +
  xlab('Alu-like pairs')+theme(axis.title.x = element_text(size=25),
                                axis.title.y = element_text(size=25),
                                axis.text.x = element_text(angle=60,hjust=1,size=15),
                                axis.text.y = element_text(size=18),
                                legend.title = element_text(size=18),
                                legend.text = element_text(size=15))
#####################################
## Fig 6
#####################################
# Significant circQTLs within Alu-like pairs 
total_snps<-read.table('./output/combined/combined.vcf',header = T)[,c('VariantID','Chrom','ARS_UCD1_2_Pos','ConsequenceType')]
intron_snps<-total_snps[total_snps$ConsequenceType=='intron_variant',]
intron_snps$ars_pos<-as.numeric(unlist(lapply(intron_snps$ARS_UCD1_2_Pos,function(x){strsplit(x,':',fixed = T)[[1]][2]})))
total_filt<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)[,c(1,3,4,7,10,11)]
intron_circqtls<-merge(total_filt,intron_snps,by='VariantID')

# Alu-like/Non
alu_info_final<-read.table('./output/Alu_iden/alu_info_final.txt',header=T,sep='\t')
alu_info_filt<-alu_info_final[setdiff(as.numeric(rownames(alu_info_final)),unique(filt_idx)),]

# stat number of each interaction part
stat_n_snps<-function(x,ars_pos,id,type){
  library(stringr)
  # x is a vector
  
  if (id == 'total') {
    x<-strsplit(x,',')[[1]]
    output<-unlist(lapply(x,function(x){chr<-strsplit(x,':',fixed=T)[[1]][1]
    tmp<-strsplit(x,':',fixed=T)[[1]][2];start<-strsplit(tmp,'-',fixed = T)[[1]][1];end<-strsplit(tmp,'-',fixed = T)[[1]][2]
    out<-sum(ars_pos$Chrom==chr&ars_pos$ars_pos>=start&ars_pos$ars_pos<=end,na.rm = T);return(out)}))
    
    return(sum(output))
  } else if (id == 'each'){
    chr<-strsplit(x,':',fixed=T)[[1]][1]
    tmp<-strsplit(x,':',fixed=T)[[1]][2];start<-strsplit(tmp,'-',fixed = T)[[1]][1];end<-strsplit(tmp,'-',fixed = T)[[1]][2]
    
    if (type == 'id') {
      output<-ars_pos$VariantID[ars_pos$Chrom==chr&ars_pos$ars_pos>=start&ars_pos$ars_pos<=end]
    } else {
      output<-sum(ars_pos$Chrom==chr&ars_pos$ars_pos>=start&ars_pos$ars_pos<=end,na.rm = T)
    }
  return(output)
  }
}
res_snps<-unlist(parLapply(cl,alu_info_filt$alu_loc_re,stat_n_snps,intron_snps,'total'))
alu_info_filt$n_alu_snps<-res_snps
alu_info_re<-data.frame('sine_name'=unlist(strsplit(alu_info_filt$sine_name,',',fixed = T)),'sine_loc'=unlist(strsplit(alu_info_filt$alu_loc_re,',',fixed = T)))
alu_info_re<-alu_info_re[!duplicated(alu_info_re$sine_loc),]
# stat number of snps within each alu-like element
alu_info_re$n_snps<-unlist(parLapply(cl,alu_info_re$sine_loc,stat_n_snps,intron_snps,'each','n'))
idx<-alu_info_re$n_snps>0;alu_info_re$snps<-'NA'
alu_info_re$snps[idx]<-parLapply(cl,alu_info_re$sine_loc[idx],stat_n_snps,intron_snps,'each','id')
alu_info_re$snps[idx]<-unlist(lapply(alu_info_re$snps[idx],function(x){paste(na.omit(unlist(x)),collapse = ',')}))
alu_info_re$snps<-unlist(alu_info_re$snps)

alu_info_re$type<-unlist(lapply(alu_info_re$sine_name,function(x){a1=!is.na(str_match(x,'Bov'))
   a2=!is.na(str_match(x,'BOV'))
   if (a1|a2){return('Alu_like')} else {return('Others')}}))
# write.table(alu_info_re,'./output/circRNA_mRNA/alu_info_re.txt',quote = F,row.names = F,sep='\t')
alu_info_re<-read.table('./Correlation in whole/alu_info_re.txt',header = T,sep='\t')
alu_info_re$QTL[idx]<-lapply(alu_info_re$snps[idx],function(x){query<-strsplit(x,',')[[1]]
  output<-unlist(lapply(query,function(x){
    if (x %in% intron_circqtls$VariantID) {
      return('circQTL')} else {return('Non-circQTL')}
    }))
  })
alu_info_re$QTL<-unlist(lapply(alu_info_re$QTL,function(x){paste(na.omit(unlist(x)),collapse = ',')}))
alu_like_circqtl<-unlist(strsplit(unlist(alu_info_re$snps[alu_info_re$type=='Alu_like'&alu_info_re$n_snps>0]),','))
other_circqtl<-unlist(strsplit(unlist(alu_info_re$snps[alu_info_re$QTL=='circQTL'&alu_info_re$n_snps>0]),','))
total_pop<-length(unique(unlist(strsplit(unlist(alu_info_re$snps[alu_info_re$n_snps>0]),','))))
alu_res<-dcast(alu_info_re,sine_name~QTL,value.var = 'n_snps',fun.aggregate = sum)
alu_res$circQTL<-alu_res$circQTL+alu_res$`circQTL,circQTL`
alu_res<-alu_res[alu_res$circQTL!=0|alu_res$`Non-circQTL`!=0,]
alu_res_filt<-melt(alu_res[,c(1,3,5)],id.vars = 'sine_name',variable.name = 'SINEs',value.name = 'Number')

## Fig 6B
# ggplot(alu_res_filt)+geom_bar(aes(x=sine_name,y=Number))
ggplot(alu_res_filt)+geom_bar(aes(x=reorder(sine_name, -Number),y=Number,fill=SINEs),stat = 'identity')+
  theme_bw()+xlab('SINEs')+guides(fill=guide_legend(title="SNPs"))+
  scale_fill_manual(values = c('#fc7171','cornflowerblue'))+
  theme(axis.title.x = element_text(size=18),
                   axis.title.y = element_text(size=18),
                   axis.text.x = element_text(angle=45,hjust=1,size=18),
                   axis.text.y = element_text(size=15),
                   legend.title = element_text(size=15),
                   legend.text = element_text(size=12))
## Fig 6A
# venn
venn.diagram(list('Alu-like element'=alu_like_circqtl,'circQTLs'=other_circqtl),
             filename = './ABS_event/Alu_iden/distribution of circQTLs within Alu-like elements.png',hyper.test = T,
             total.population=total_pop,sub.cex=2,cex=2,
             lwd = 4,cat.cex=2,cat.pos=c(-20,50),
             fill = c("cornflowerblue", "red"))
# p = 0.33 indicates no significant enrichment of circQTLs
# ggplot(alu_info_re[alu_info_re$n_snps>0,])+geom_boxplot(aes(x=type,y=n_snps,fill=type))
VennDiagram::draw.pairwise.venn(12,25,2,category = c("First", "Second"),
                                cat.pos = c(0, 180),
                                euler.d = TRUE,
                                sep.dist = 0.03,
                                rotation.degree = 45,fill=c('red','green'))

# SINEs in Alu-like
total_filt$ARS_pos<-unlist(lapply(total_filt$ARS_UCD1_2_Pos,function(x){strsplit(x,':',fixed = T)[[1]][2]}))

# ABS expr count matrix
circRNA_expr<-read.table('./output/ABS_event/ABS_events_readscount.txt',header=T,sep = '\t')

circRNA_expr$circ<-paste0(circRNA_expr$chrom,':',circRNA_expr$start,'|',circRNA_expr$end)
circRNA_expr<-circRNA_expr[,c('circ','ABS','ABS_class')]
# rownames(blast_res_filt)<-gsub(x=rownames(blast_res_filt),pattern = 'V',replacement = '')
circRNA_expr_csi<-merge(circRNA_expr,blast_res_filt[,c('circ','CSI')],by='circ')
# Differernce of CSI among different ABS events
circRNA_expr_csi_res<-melt(dcast(circRNA_expr_csi,ABS_class~ABS,value.var = 'CSI',fun.aggregate = sum),
                           id.vars = 'ABS_class',variable.name = 'ABS',value.name = 'CSI')
## Fig 6C
# set ABS colors
ABS=pal_nejm()(8)[1:4]
names(ABS)<-c("A5BS & A3BS","A5BS",'A3BS','not ABS')
# Cumulative freq of CSI of ABS events
circRNA_expr_csi_res$log2CSI<-log2(abs(circRNA_expr_csi_res$CSI))
types<-circRNA_expr_csi_res$ABS[!duplicated(circRNA_expr_csi_res$ABS)]
for (i in 1:length(types)) {
  # i=1
  tmp<-circRNA_expr_csi_res[circRNA_expr_csi_res$ABS==types[i],]
  tmp<-tmp[order(tmp$log2CSI),]
  idx<-rownames(tmp)
  circRNA_expr_csi_res[idx,'Frac']<-seq(0+1/length(tmp[,1]),1,1/length(tmp[,1]))
}
ggplot(circRNA_expr_csi_res[circRNA_expr_csi_res$ABS!='not ABS',])+geom_line(aes(x=log2CSI,y=Frac,color=ABS),size=1.5)+
  scale_color_manual(values = ABS[1:3])+
  theme_bw()+xlab('Sum of CSI for each ABS event')+ylab('Cumulative frequency of circRNAs')+ylim(0.6,1)+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

# Cumulative freq indicates the higher divergence of CSI for ABS event types in lower CSI region. and as a result of unstability of 
# A5BS & A3BS, it needs stronger CSI to form circRNA compared to other ABS types.


## Fig 6D
# ggplot(circRNA_expr_csi_res)+geom_violin(aes(x=ABS,y=log2(CSI),fill=ABS))+
#   geom_boxplot(aes(x=ABS,y=log2(CSI),fill=ABS))
# box and vilon plot
ggbetweenstats(
  data = circRNA_expr_csi_res[circRNA_expr_csi_res$CSI>0,],
  x=ABS,y=log2CSI,
  type='p',
  centrality.plotting=F,
  messages = F
)+theme_bw()+
  xlab('ABS events')+theme(axis.title.x = element_text(size=18),
                               axis.title.y = element_text(size=18),
                               axis.text.x = element_text(angle=45,hjust=1,size=18),
                               axis.text.y = element_text(size=15))
##
intron_len_dis<-read.table('./ABS_event/Alu_iden/intron_len_dis.txt',header=T,sep='\t')
colnames(intron_len_dis)[1]<-'circ'
tmp<-intron_len_dis[,c('circ','ABS')]
intron_len_dis$Complexity<-apply(tmp,1,function(x){max(intron_len_dis$Complexity[intron_len_dis$circ==x[1]&
                                                                                   intron_len_dis$ABS==x[2]])})
circRNA_expr_csi<-merge(circRNA_expr_csi,intron_len_dis[,c('circ','ABS','Complexity')],by=c('circ','ABS'))
circRNA_expr_csi<-circRNA_expr_csi[!duplicated(circRNA_expr_csi[,c('circ','ABS')]),]

## Fig 6D
# Difference of CSI in different Complexity of ABS events
circRNA_expr_csi$Complexity[circRNA_expr_csi$Complexity>=3]<-paste0(circRNA_expr_csi$ABS[circRNA_expr_csi$Complexity>=3],'(>=3)')
circRNA_expr_csi$Complexity[circRNA_expr_csi$Complexity==2]<-paste0(circRNA_expr_csi$ABS[circRNA_expr_csi$Complexity==2],'(=2)')
circRNA_expr_csi$Complexity[circRNA_expr_csi$Complexity==0]<-'control'
circRNA_expr_csi$Complexity<-factor(circRNA_expr_csi$Complexity,levels = c("A5BS(>=3)","A5BS(=2)", 'A3BS(>=3)',"A3BS(=2)","control"))
# ggplot(circRNA_expr_csi)+geom_boxplot(aes(x=Complexity,y=log2(CSI),fill=Complexity))
# set colors
ABS_complexity<-rep(ABS[2:3],each=2)
names(ABS_complexity)<-c("A5BS(>=3)","A5BS(=2)", 'A3BS(>=3)',"A3BS(=2)")
circRNA_expr_csi$log2CSI<-log2(circRNA_expr_csi$CSI)
ggbetweenstats(data = circRNA_expr_csi[circRNA_expr_csi$Complexity!='control',],
               x=Complexity,y=log2CSI,
               type='p',
               centrality.plotting=F,
               messages = F)+theme_bw()+
  scale_color_manual(values = ABS_complexity)+
  # guides(color=guide_legend())
  scale_fill_manual(values = ABS_complexity)+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(angle=30,hjust=1,size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

##Fig 15C
abs_class<-na.omit(unique(circRNA_expr_csi$ABS_class[circRNA_expr_csi$ABS!='A5BS & A3BS'&
                                                           is.na(str_match(circRNA_expr_csi$ABS_class,'NA'))]))
abs_isforms_expr<-data.frame('ABS_class'=rep(abs_class,each=2))
# col=colnames(circRNA_expr)[6:length(circRNA_expr)]
circRNA_expr_csi$chrom<-unlist(lapply(circRNA_expr_csi$circ,function(x){tmp<-strsplit(x,':')[[1]][1]
return(tmp)}))
circRNA_expr_csi$start<-as.numeric(unlist(lapply(circRNA_expr_csi$circ,function(x){tmp<-strsplit(x,':')[[1]][2]
  return(strsplit(tmp,'|',fixed = T)[[1]][1])})))
circRNA_expr_csi$end<-as.numeric(unlist(lapply(circRNA_expr_csi$circ,function(x){tmp<-strsplit(x,':')[[1]][2]
  return(strsplit(tmp,'|',fixed = T)[[1]][2])})))
##
circRNA_expr<-read.table('./output/ABS_event/ABS_events_readscount.txt',header=T,sep = '\t')
circRNA_expr$circ<-paste0(circRNA_expr$chrom,':',circRNA_expr$start,'|',circRNA_expr$end)

circRNA_expr<-subset(circRNA_expr,select = -c(ABS,ABS_class))
circRNA_expr_csi<-merge(circRNA_expr_csi,circRNA_expr,by=c('chrom','start','end'))
breeds_info<-read.table('./SRR_breeds.CSV',header = T,sep=',')
srr<-breeds_info$SRR_list
circRNA_expr_csi$mean_expr<-rowMeans(circRNA_expr_csi[,srr])
for (i in 1:length(abs_class)) {
  # i=2
  abs_isforms_expr[(2*i-1),'type']='predominant'
  
  idx<-!is.na(str_match(circRNA_expr_csi$ABS_class,abs_class[i]))
  tmp<-circRNA_expr_csi[idx,]
  abs_isforms_expr[(2*i-1),'log2CSI']=tmp[which(tmp$mean_expr==max(tmp$mean_expr)),'log2CSI']
  abs_isforms_expr[2*i,'type']='others'
  abs_isforms_expr[2*i,'log2CSI']=mean(tmp[which(tmp$mean_expr<max(tmp$mean_expr)),'log2CSI'])
  
}
## classification of circRNA isoforms by length
for (i in 1:length(abs_class)) {
  # i=2
  abs_isforms_expr[(2*i-1),'type']='shorter'
  
  idx<-!is.na(str_match(circRNA_expr_csi$ABS_class,abs_class[i]))
  tmp<-circRNA_expr_csi[idx,]
  abs_isforms_expr[(2*i-1),'log2CSI']=tmp[which(tmp$end-tmp$start==min(tmp$end-tmp$start)),'log2CSI']
  abs_isforms_expr[2*i,'type']='longer'
  abs_isforms_expr[2*i,'log2CSI']=mean(tmp[which(tmp$end-tmp$start>min(tmp$end-tmp$start)),'log2CSI'],na.rm = T)
  
}
# set color
abs_isforms_expr$type<-factor(abs_isforms_expr$type,levels=c('shorter','longer'))
isoform_color<-pal_jama()(10)[c(3,4)]
names(isoform_color)<-c('longer','shorter')
## Fig 6I
ggbetweenstats(data = abs_isforms_expr,
               x=type,y=log2CSI,
               type='p',
               centrality.plotting=T,
               messages = F)+theme_bw()+
  # scale_x_discrete(breaks = c('shorter','longer'),labels = c('shorter','longer'))+
  scale_color_manual(values = isoform_color)+
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size=15))

