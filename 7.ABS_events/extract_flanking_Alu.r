#!/bin/Rscript

library(ggplot2)
library(stringr)
library(parallel)
library(reshape2)
library(ggpubr)
library(ggstatsplot)
library(RColorBrewer)
library(VennDiagram)

############################################################################################################
working_dir <- "path_to_your_working_dir"
all_RNA_seq_file='./SRR_list.txt'
output_dir <- './output/Alu_iden'
cores <- 6
setwd(working_dir)
dir.create(output_dir)
cl<-makeCluster(cores)


# CSI calculate
alu_info<-read.table('./output/Alu_iden/alu_info_sine.txt',sep='\t',header = T)
alu_info$circ<-unlist(lapply(alu_info$name,function(x){strsplit(x,'_',fixed = T)[[1]][1]}))
filt_loc<-unlist(lapply(alu_info$name[alu_info$Alu_number=="0.0"],function(x){strsplit(x,'_',fixed = T)[[1]][1]}))
filt_idx<-unlist(lapply(filt_loc,function(x){which(alu_info$circ==x)}))
unique(filt_idx)
write.table(alu_info,'./output/Alu_iden/alu_info_final.txt',row.names=F,quote=F,sep='\t')

alu_info_final<-read.table('./output/Alu_iden/alu_info_final.txt',header=T,sep='\t')
alu_info_filt<-alu_info_final[setdiff(as.numeric(rownames(alu_info_final)),unique(filt_idx)),]
alu_info_filt<-alu_info_filt[alu_info_filt$chrom!='chrom',]
alu_info_filt$alu_loc<-paste(alu_info_filt$chrom,alu_info_filt$Alu_start,alu_info_filt$Alu_end,sep=':')

stats_sine<-function(x){
  chr<-strsplit(x,':',fixed = T)[[1]][1];chr<-str_replace(chr,'chr','')
  start_tmp<-strsplit(x,':',fixed = T)[[1]][2]; end_tmp<-strsplit(x,':',fixed = T)[[1]][3]
  if (start_tmp==''){
    return(NA)
  }else{
    tmp<-data.frame('start'=unlist(strsplit(start_tmp,', ')),'end'=unlist(strsplit(end_tmp,', ')))
    tmp$chr<-chr;output<-paste0(tmp$chr,':',tmp$start,'-',tmp$end)
    return(output)
  }
  
}

alu_info_filt$alu_loc_re<-lapply(alu_info_filt$alu_loc,stats_sine)
alu_info_filt$alu_loc_re<-unlist(lapply(alu_info_filt$alu_loc_re,function(x){str_c(x,collapse = ',')}))
alu_info_filt$Flanking<-unlist(lapply(alu_info_filt$name,function(x){strsplit(x,'_',fixed = T)[[1]][2]}))
df<-alu_info_filt[,c('circ','Flanking','alu_loc_re')]
a<-tidyr::pivot_wider(df,id_cols='circ',names_from='Flanking',values_from='alu_loc_re')
write.table(a,'./output/Alu_iden/alu_flank.txt',row.names=F,quote=F,sep='\t')
