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

#####################################
## Fig 4
#####################################
alu_info_final<-read.table('./output/Alu_iden/alu_info_final.txt',header=T,sep='\t')

# stat number of different sines
stats_n<-function(x){
  # x is a list
  x<-strsplit(x,',',fixed=T)[[1]]
  output<-list()
  # which(x==sine_types)
  for (i in 1:length(sine_types)){
    res<-which(unlist(x)==sine_types[i])
    if (length(res)==0){
      output<-append(output,0)
    }else{
      output<-append(output,length(unlist(res)))
    }
  }
  return(output)
}
# stats_n(alu_info_final$sine_name[1])
b<-lapply(alu_info_final$sine_name,stats_n)
c<-matrix(unlist(b),ncol=length(sine_types),byrow = T)
colnames(c)<-sine_types
alu_info_final<-cbind(alu_info,c)

## Fig 4A-B
# Frequency of different SINEs in up/downstream
df<-alu_info_final[,c('name','Alu_number')]
df<-cbind(df,c)
df$name<-unlist(lapply(df$name,function(x){strsplit(x,'_',fixed = T)[[1]][2]}))
df<-melt(df,id.vars = c('name','Alu_number'),variable.name = 'SINE_types',value.name = 'Number')
df<-df[df$Number!=0,]
df_show<-data.frame()
for (j in 1:2) {
  flank_type<-unique(df$name)[j]
for (i in 1:length(sine_types)) {
  df_show[((j-1)*length(sine_types)+i),'flanking']<-flank_type;df_show[((j-1)*length(sine_types)+i),'SINE_types']<-sine_types[i]
  df_show[((j-1)*length(sine_types)+i),'Alu_number']<-sum(as.numeric(df[df$name==flank_type&df$SINE_types==sine_types[i],'Alu_number']))/
    sum(as.numeric(df[df$name==flank_type,'Alu_number']))
  df_show[((j-1)*length(sine_types)+i),'Number']<-sum(df[df$name==flank_type&df$SINE_types==sine_types[i],'Number'])/
    sum(as.numeric(df[df$name==flank_type,'Number']))
}
}
library(ggsci)
df_show$SINE<-df_show$SINE_types
df_show$SINE[df_show$Number<0.01]<-'Others'
df_show$SINE<-paste0(df_show$SINE,'(',round(df_show$Number,3)*100,'%)')
df_show_filt<-df_show[df_show$Number>=0.01,]
df_show$SINE[df_show$Number<0.01]<-'Others'
df_show$SINE[df_show$SINE=='Others']<-paste0('Others','(',round(sum(df_show$Number[df_show$Number<0.01]),3)*100,'%)')

# ggpie(df_show[df_show$flanking=='Upstream',],"Number",
      # label='SINE',lab.font=c(0, "bold", "red"),
      # fill = "SINE")
ggplot(df_show[df_show$flanking=='Upstream',],aes(x='',y=Number,fill=SINE))+geom_bar(stat = 'identity')+coord_polar(theta = 'y')+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(face = 'italic',size = 12),
        legend.title = element_text(face = 'italic',size = 15),
        legend.position = 'bottom',
        axis.title = element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(12)))

ggplot(df_show[df_show$flanking=='Downstream',],aes(x='',y=Number,fill=SINE))+geom_bar(stat = 'identity')+coord_polar(theta = 'y')+
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(face = 'italic',size = 12),
        legend.title = element_text(face = 'italic',size = 15),
        legend.position = 'bottom',
        axis.title = element_blank(),axis.ticks = element_blank())+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(12)))
# dev.off()


############################################################################################################
# deal with the control
intron_total<-read.table('output/Alu_iden/intron_total.txt',header=T,sep='\t')
alu_total<-read.table('output/Alu_iden/alu_total.txt',header=T,sep='\t')

alu_total$length<-alu_total$end-alu_total$start
introns_flanking<-intron_total[intron_total$Flanking!=1,c('chrom','start','end','Flanking')]
introns_flanking$name<-paste0(introns_flanking$chrom,':',introns_flanking$start,'-',introns_flanking$end)
alu_total_filt<-merge(introns_flanking,alu_total,by=c('chrom','start','end'))

# calculate the alu number of introns
introns_class<-unique(alu_total_filt$name.x)
intron_alu_data<-data.frame('intron_class'=introns_class)
alu_total_filt$count<-1
alu_like<-c('Bov-tA1','Bov-tA2','Bov-tA3','BOV-A2')
analyze<-function(){
  for (i in 1:length(introns_class)) {
    print(i)
    intron_info<-introns_class[i]
    up<-list();down<-list()
    for (j in 1:length(alu_like)) {
      a5bs<-sum(alu_total_filt$count[alu_total_filt$name.x==intron_info&
                                       alu_total_filt$Flanking=='Upstream'&alu_total_filt$name.y==alu_like[j]])
      a3bs<-sum(alu_total_filt$count[alu_total_filt$name.x==intron_info&
                                       alu_total_filt$Flanking=='Downstream'&alu_total_filt$name.y==alu_like[j]])
      up<-append(up,a5bs);down<-append(down,a3bs)
    }
    
    intron_alu_data$Upstream[i]<-list(unlist(up));intron_alu_data$Downstream[i]<-list(unlist(down))
  }
}
library(compiler)
com.analyze<-cmpfun(analyze)
com.analyze()

# left<-matrix(unlist(intron_alu_data$Upstream),ncol=length(alu_like),byrow = T);alu_len_control_left<-cbind(intron_alu_data[,1],left)
right<-matrix(unlist(intron_alu_data$Downstream),ncol=length(alu_like),byrow = T)
left<-matrix(unlist(intron_alu_data$Upstream),ncol=length(alu_like),byrow = T)
alu_len_control_left<-as.data.frame(cbind(intron_alu_data[,1],left,right))
colnames(alu_len_control_left)<-c('intron_class',alu_like,paste0(alu_like,'.Down'))
# add total number info of tRNA-Core-RTE (Bov-tA family)
alu_len_control_left$`Bov-tA family`<-rowSums(apply(alu_len_control_left[,2:4],2,as.numeric));alu_len_control_left$`Bov-tA family.Down`<-rowSums(apply(alu_len_control_left[,6:9],2,as.numeric))
alu_len_control_left$`Total`<-rowSums(apply(alu_len_control_left[,2:5],2,as.numeric));alu_len_control_left$`Total.Down`<-rowSums(apply(alu_len_control_left[,6:8],2,as.numeric))

alu_len_control<-melt(alu_len_control_left,id.vars = 'intron_class',variable.name = 'Alu-like',value.name = 'Alu_number')
alu_len_control$Flanking<-'Upstream';
a=str_match(string=alu_len_control$`Alu-like`,pattern = '.Down')
alu_len_control$Flanking[which(a[,1]=='.Down')]='Downstream'
alu_len_control$`Alu-like`[which(a[,1]=='.Down')]<-str_replace(alu_len_control$`Alu-like`[which(a[,1]=='.Down')],".Down",'')

# alu_len_control<-melt(intron_alu_data,id.vars = 'intron_class',variable.name = 'Flanking',value.name = 'Alu_number')
alu_len_control$ABS<-'control';alu_len_control$Complexity=0;colnames(alu_len_control)[1]<-'name'

# alu_info
alu_len_dis<-alu_info_final[alu_info_final$chrom!='chrom',c('name',alu_like)]
alu_len_dis$Flanking<-unlist(lapply(alu_len_dis$name,function(x){strsplit(x,'_',fixed = T)[[1]][2]}))
alu_len_dis$name<-unlist(lapply(alu_len_dis$name,function(x){strsplit(x,'_',fixed = T)[[1]][1]}))
alu_len_dis$`Bov-tA family`<-rowSums(alu_len_dis[,2:4]);alu_len_dis$`Total`<-rowSums(alu_len_dis[,2:5])

alu_len_dis_abs=merge(alu_len_dis,intron_len_dis,by.x=c('name','Flanking'))
alu_len_dis_abs<-alu_len_dis_abs[,-9]
alu_len_dis_abs<-melt(alu_len_dis_abs,id.vars = c('name','Flanking','ABS','Complexity'),variable.name = 'Alu-like',value.name = 'Alu_number')
alu_len_dis_abs<-alu_len_dis_abs[,c(1,5,6,2,3,4)]
# output data
alu_len_combine<-rbind(alu_len_dis_abs,alu_len_control)
alu_len_combine<-alu_len_combine[alu_len_combine$Complexity!='not ABS',]
class<-unique(circRNA$ABS)[-2];tmp<-alu_len_combine[,c('name','ABS')]
alu_len_combine$Complexity<-apply(tmp,1,function(x){max(alu_len_combine$Complexity[alu_len_combine$name==x[1]&
                                                                                     alu_len_combine$ABS==x[2]])})
alu_len_combine<-alu_len_combine[!(alu_len_combine$Complexity==0&alu_len_combine$ABS!='control'),]
alu_len_combine$X<-'control';alu_len_combine$X[alu_len_combine$Complexity>=3]<-'≥3';
alu_len_combine$X[alu_len_combine$Complexity==2]<-'=2'
alu_len_combine$X[alu_len_combine$Complexity==1]<-'=1'
alu_len_combine$X<-paste0(alu_len_combine$ABS,' ',alu_len_combine$X)
alu_len_combine$X<-factor(alu_len_combine$X,levels = c('A5BS ≥3','A5BS =2', 'A3BS ≥3', 'A3BS =2','A53BS ≥3', 'A53BS =2','A53BS =1','control'))
alu_len_combine$Flanking<-factor(alu_len_combine$Flanking,levels = c('Upstream','Downstream'))
alu_len_combine$X[is.na(alu_len_combine$X)]<-'control'
alu_len_combine$Alu_number<-as.numeric(alu_len_combine$Alu_number)
# write.table(alu_len_combine,'./output/Alu_iden/alu_len_combine.txt',quote = F,row.names = F,sep='\t')
# alu_len_combine<-alu_len_combine[alu_len_combine$Alu_number!=0,]
# t test
alu_len_combine<-read.table('./output/Alu_iden/alu_len_combine.txt',header = T,sep = '\t')
# alu_len_combine<-alu_len_combine[alu_len_combine$Alu_number!=0,]
colnames(alu_len_combine)[2]<-'Alu-like'
alu_len_combine$`Alu-like`<-factor(alu_len_combine$`Alu-like`,levels = unique(alu_len_combine$`Alu-like`))

flanking<-unique(alu_len_combine$Flanking);Com<-unique(alu_len_combine$X)
res<-list()
k=2
for (i in 1:length(flanking)) {
  # i=1
  tmp<-flanking[i]
  for (j in (1:(length(Com)-1))) {
    # j=2
    each<-list()
    for (h in (j+1):length(Com)) {
      # h=j+2
      sampleA<-Com[j];sampleB<-Com[h];
      contrast<-paste0(sampleA,' vs ',sampleB)
      x<-alu_len_combine$Alu_number[alu_len_combine$Flanking==tmp&alu_len_combine$X==sampleA&alu_len_combine$`Alu-like`==alu_like[k]]
      y<-alu_len_combine$Alu_number[alu_len_combine$Flanking==tmp&alu_len_combine$X==sampleB&alu_len_combine$`Alu-like`==alu_like[k]]
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
res
names(res)<-flanking
# plot the boxplot of Alu_number within ABS circRNAs
## Figure 4C
ggstatsplot::ggbetweenstats(
  data = alu_len_combine,
  x=X,y=Alu_number,grouping.var = Flanking,
  messages = FALSE
)

############################################################################################################
################################
## Figure S14
################################
ggplot(alu_len_combine)+geom_boxplot(aes(x=X,y=Alu_number,fill=Flanking),outlier.alpha = 0)+ylim(0,10)+facet_wrap(~`Alu-like`)+
  theme_bw()+xlab('Complexity')+ylab('Alu-like number')+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=18),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        strip.text = element_text(size = 12))

