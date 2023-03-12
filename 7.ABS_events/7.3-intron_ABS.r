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

## Basic statistic info of Alu-like elements within flanking introns of ABS circRNAs
# import related info
circRNA<-read.table('./output/ABS_event/ABS_events.txt',sep='\t',header = T)
total_complex_info<-circRNA[!is.na(circRNA$ABS_class),1:5]
# import intron_total & intron_info
intron_total<-read.table('./output/Alu_iden/introns.gtf',sep='\t')
colnames(intron_total)<-c('chrom','database','type','start','end','na','strand','b','feature')
intron<-read.table('./ABS_event/Alu_iden/introns_info.txt',sep='\t',header = T)
intron$name<-paste0(intron$chrom,':',intron$start,'|',intron$end)
intron$length<-intron$end_intron-intron$start_intron
abs_class<-unique(circRNA$ABS_class)[-2] # no NA(not ABS)
# DES of intron length and alu elements within up/downstream for ABS circRNAs
# intron length in control (total genome)
intron_total$length<-intron_total$end-intron_total$start
intron_total$geneid<-unlist(lapply(intron_total$feature,function(x){tmp<-strsplit(x,';',fixed = T)[[1]][1]
                                                            return(strsplit(tmp,' ',fixed = T)[[1]][2])}))
idx_up<-unique(unlist(parLapply(cl,intron_total$geneid,function(x,y){which(y==x)[1]},intron_total$geneid)))
idx_down<-unique(unlist(parLapply(cl,intron_total$geneid,function(x,y){which(y==x)[length(which(y==x))]},intron_total$geneid)))
intron_total$Flanking<-1
intron_total$Flanking[idx_up]<-'Upstream';intron_total$Flanking[idx_down-1]<-'Downstream'
intron_filt<-intron_total[intron_total$Flanking=='Upstream'|intron_total$Flanking=='Downstream',]
intron_len_control<-data.frame('name'=NA,'Flanking'=intron_filt$Flanking,'ABS'='control','length'=intron_filt$length,'Complexity'=0)
write.table(intron_total,'output/Alu_iden/intron_total.txt',quote = F,row.names = F,sep='\t')
############################################################################################################
#####################################
## Figure 4
#####################################
# intron length in ABS circRNAs
total_complex_data<-merge(intron,circRNA[circRNA$ABS!='not ABS',1:5],by.x=c('chrom','start','end'))
total_complex_info$count<-1
for (i in 1:length(abs_class)) {
  # i=1
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
  total_complex_data$A53BS[i]<-a53bs
}
# total_complex_data<-melt(total_complex_data,id.vars = 'ABS_class',variable.name = 'ABS_type',
#                          value.name = 'Complexity')
# head(total_complex_data)
intron_len_dis<-melt(total_complex_data[,c('name','Flanking','A5BS','A3BS','A53BS','length')],
                     id.vars = c('name','Flanking','length'),variable.name = 'ABS',value.name = 'Complexity')
write.table(intron_len_dis,'./output/ABS_event/Alu_iden/intron_len_dis.txt',quote = F,row.names = F,sep='\t')
# 
intron_len_dis<-read.table('./output/ABS_event/Alu_iden/intron_len_dis.txt',sep='\t')
head(intron_len_dis)
intron_len_combine<-rbind(intron_len_dis,intron_len_control)
# intron_len_combine<-intron_len_combine[intron_len_combine$Complexity!='not ABS',]
class<-unique(circRNA$ABS)[-2];tmp<-intron_len_combine[,c('name','ABS')]
intron_len_combine$Complexity<-apply(tmp,1,function(x){max(intron_len_combine$Complexity[intron_len_combine$name==x[1]&
                                                                                           intron_len_combine$ABS==x[2]])})
intron_len_combine<-intron_len_combine[!(intron_len_combine$Complexity==0&intron_len_combine$ABS!='control'),]
intron_len_combine$X<-'control';intron_len_combine$X[intron_len_combine$Complexity>=3]<-'≥3';
intron_len_combine$X[intron_len_combine$Complexity==2]<-'=2'
intron_len_combine$X[intron_len_combine$Complexity==1]<-'=1'
intron_len_combine$X<-paste0(intron_len_combine$ABS,' ',intron_len_combine$X)
intron_len_combine$X<-factor(intron_len_combine$X,levels = c('A5BS ≥3','A5BS =2', 'A3BS ≥3', 'A3BS =2','A53BS ≥3', 'A53BS =2','A53BS =1','control'))
intron_len_combine$Flanking<-factor(intron_len_combine$Flanking,levels = c('Upstream','Downstream'))
intron_len_combine$X[is.na(intron_len_combine$X)]<-'control'
# t test
flanking<-unique(intron_len_combine$Flanking);Com<-levels(intron_len_combine$X)
res<-list()
for (i in 1:length(flanking)) {
  # i=2
    tmp<-flanking[i]
    for (j in (1:(length(Com)-1))) {
      # j=2
      each<-list()
      for (h in (j+1):length(Com)) {
        # h=j+2
      sampleA<-Com[j];sampleB<-Com[h];
      contrast<-paste0(sampleA,' vs ',sampleB)
      x<-intron_len_combine$length[intron_len_combine$Flanking==tmp&intron_len_combine$X==sampleA]
      y<-intron_len_combine$length[intron_len_combine$Flanking==tmp&intron_len_combine$X==sampleB]
      x_noout <- x[-which(x%in%boxplot.stats(x)$out)];y_noout<-y[-which(y%in%boxplot.stats(y)$out)]
      
      if (is.na(x[1])){
        next
      }
      result<-t.test(x_noout,y_noout,alternative = 'two.sided')
      # result<-t.test(x,y,alternative = 'two.sided')
      if (h==(j+1)) {
        each[[1]]<-result$p.value
        # name[[i]]<-tmp
      }
      each[[1]][(h-j)]<-result$p.value
      names(each[[1]])[(h-j)]<-contrast
      
      }
      if (j ==1) {
        res[[j]]<-each
      }
      res[[j]]<-each
    }
    
}
names(res)<-flanking
write.table(intron_len_combine,'./figures_data/intron_len_combine.txt',quote = F,row.names = F)
# plot the boxplot of intron length within ABS circRNAs
intron_len_combine<-read.table('./figures_data/intron_len_combine.txt')
## Figure 4D
ggplot(intron_len_combine)+geom_boxplot(aes(x=X,y=length,fill=Flanking),outlier.alpha = 0)+ylim(0,6e4)+
  theme_bw()+xlab('Complexity')+ylab('intron length (bp)')+
  theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(angle=45,hjust=1,size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15))+
  annotate('text',x = 1.3,y=4.2e4,label=paste0('p =',format(res[[1]][[1]][1],digits = 3,scientific = T)),size=5)+
  annotate('segment', x=0.8, xend=1.8, y=4.0e4, yend=4.0e4,size=1,
           arrow=arrow(angle = 0,length = unit(0,'cm')))


# alu elements
# import control alu elements
alu_total<-read.table('./output/Alu_iden/bt9_SINE.gtf',sep='\t')
colnames(alu_total)<-c('chrom_alu','start_alu','end_alu','name','strand_alu','chrom','database','type','start','end','na','strand','b','feature')

# import alu info within flanking introns of ABS circRNA
alu_info<-read.table('./output/Alu_iden/alu_info.bed',sep='\t',header = T)
alu_info<-alu_info[,-1]

# statistic the number of different SINEs
sine_types<-unique(alu_total$name)
sine_type_stat<-alu_info
# analyze the number of all SINEs for each one
# stats_sine<-function(chr,start,end){
  
stats_sine<-function(X,alu_total,id){
  library(stringr)
  # start and end is vectors, start: c(start1,start2,..), end: just like start.
  # output is vector.
  chr<-X[1];start=X[2];end=X[3]
  if (start=='') {
    return(NA)
  }
  start<-as.numeric(strsplit(start,', ',fixed = T)[[1]]);end<-as.numeric(strsplit(end,', ',fixed = T)[[1]])
  chr<-stringr::str_sub(chr,4,str_count(chr))
 
  for (i in 1:length(start)) {
    if (id == 0) {
      res<-alu_total$name[alu_total$chrom_alu==chr&alu_total$start_alu==start[i]&alu_total$end_alu==end[i]][1]
    } else {
      res<-alu_total$strand[alu_total$chrom_alu==chr&alu_total$start_alu==start[i]&alu_total$end_alu==end[i]][1]
    }

    output<-append(output,res)
  }
  return(unlist(output))
}
chr<-alu_info$chrom[1];start<-alu_info$Alu_start[1];end<-alu_info$Alu_end[1]
# stats_sine(chr,start,1)
library(compiler)
com.stats_sine<-cmpfun(stats_sine)

alu_info$sine_name<-parApply(cl,alu_info[,c('chrom','Alu_start','Alu_end')],FUN = stats_sine,1,alu_total,0)

alu_info$sine_name<-unlist(lapply(alu_info$sine_name,function(x){paste(x,collapse = ',')}))

# add strand info of SINEs 
alu_info$sine_strand<parApply(cl,alu_info[,c('chrom','Alu_start','Alu_end')],FUN = stats_sine,1,alu_total,1)
write.table(alu_info,'./output/Alu_iden/alu_info_sine.txt',sep='\t',header = T,quote=F,row.names=F)
write.table(alu_total,'./output/Alu_iden/alu_total.txt',sep='\t',header = T,quote=F,row.names=F)



