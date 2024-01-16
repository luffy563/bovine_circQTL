#!/bin/Rscript

library(MatrixEQTL)
library(CMplot)
library(ggsci)
library(ggplot2)
library(parallel)

############################################################################################################
working_dir <- "path_to_your_working_dir"
# working_dir <- "D:\\work\\master\\coderepoistory\\bovine_circQTL"

base.dir='./output/circQTLs'
setwd(working_dir)

## Create the parllel handle
cl<-makeCluster(3)

base.dir='./output/circQTL'
##Set the analysis model
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 
# set the location of SNP 
SNP_file_name = paste(base.dir, "/snp.txt", sep="");
snps_location_file_name = paste(base.dir, "/snp_pos.txt", sep="");

# set the location of expression data for circRNA
expression_file_name = paste(base.dir, "/circRNA.txt", sep="");
gene_location_file_name = paste(base.dir, "/anno_info.txt", sep="");

# set the location of output file
output_file_name_cis = './cis-circqt.txt'
output_file_name_tra = './tra-circqtl.txt'

# set the covariate
covariates_file_name = paste(base.dir, "/cov.txt", sep="");
output_file_name = 'circqtl_pre'

# pvOutputThreshold = 1e-2
pvOutputThreshold_tra = 1e-2
pvOutputThreshold_cis = 0.05
# set the error covariate
errorCovariance = numeric()

# Distance for local gene-SNP pairs
cisDist = 1e6;

# create SlicedData for storing matrix data
snps = SlicedData$new();
# set the delimiter to tab
snps$fileDelimiter = "\t";      # the TAB character
# set the missing data to NA
snps$fileOmitCharacters = "NA"; # denote missing values;

snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
# cvrt=read.csv('./output/matrixeqtl/cov.csv')
# cvrt = as.matrix(cvrt)
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;  # one column of row labels
cvrt$LoadFile(covariates_file_name)

## Run the analysis
snpspos = na.omit(read.table(snps_location_file_name, header = T, stringsAsFactors = FALSE))
genepos = na.omit(read.table(gene_location_file_name, header = T, stringsAsFactors = FALSE))

me = Matrix_eQTL_main(
   snps = snps, 
   gene = gene, 
   cvrt = cvrt,
   output_file_name = output_file_name_tra,
   pvOutputThreshold = pvOutputThreshold_tra,
   useModel = useModel, 
   errorCovariance = errorCovariance, 
   verbose = TRUE, 
   output_file_name.cis = output_file_name_cis,
   pvOutputThreshold.cis = pvOutputThreshold_cis,
   snpspos = snpspos, 
   genepos = genepos,
   cisDist = cisDist,
   pvalue.hist = "qqplot",
   min.pv.by.genesnp = FALSE,
   noFDRsaveMemory = FALSE);

## Results
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values
# plot(me)

############################################################################################################
bovine=me$trans$eqtls
bovine_filt=bovine[,c(1,4)]
colnames(bovine_filt)=c('SNP','transCircRNA')
colnames(snpspos)=c('SNP','Chromosome','Position')
bovine_filt = merge(bovine_filt,snpspos,by='SNP')
# c_maxes = bovine_filt.groupby(c('SNP','Chromosome','Position')).transCircRNA.transform(min)
# df = df.loc[df.C == c_maxes]
bovine_filt_1=bovine_filt[!duplicated(bovine_filt[,1]),]
# bovine_filt_1=bovine_filt
bovine_filt_1 = bovine_filt_1[,c(1,3,4,2)]
highlight = bovine_filt_1[,1][(bovine_filt_1[,4]<=1e-6)&(bovine_filt_1[,4]>=1e-10)]
highlight_1 = bovine_filt_1[,1][(bovine_filt_1[,4]<=1e-23)&(bovine_filt_1[,4]>=1e-24)]
highlight_2 = bovine_filt_1[,1][(bovine_filt_1[,4]<=1e-307)]

data=read.table('./output/combined/combined_data.vcf',sep=' ',header = T)
snp_info=data[,1:15]
significant = bovine_filt_1[,c(1,4)]
colnames(significant)=c('VariantID','pvalue')
significant = merge(significant,snp_info,by='VariantID')
significant=significant[order(significant$pvalue),]
# write.table(significant,'./significant_circQTL_log_quant_hfilt.txt',row.names = F,quote = F)

type <- significant$ConsequenceType[!duplicated(significant$ConsequenceType)]
nums = list()
for (i in 1:length(type)){
   nums[i]=length(significant$ConsequenceType[significant$ConsequenceType==type[i]])
}
df <- data.frame(type = type, nums = as.integer(nums))
label_value <- paste('(', round(df$nums/sum(df$nums) * 100, 1), '%)', sep = '')
# label_value
label <- paste(df$type, label_value, sep = '')

#####################################
## Figure 1A
#####################################
CMplot(bovine_filt_1,ylim=c(2,10),threshold = 1e-6,highlight=highlight,highlight.text=highlight,file="tiff", dpi=300)
CMplot(bovine_filt_1,ylim=list(c(23,24)),highlight=highlight_1,highlight.text=highlight_1,file="tiff", dpi=300,chr.border=F,plot.type = 'm')
CMplot(bovine_filt_1,ylim=list(c(307,308)),highlight=highlight_2,highlight.text=highlight_2,file="tiff", dpi=300,chr.border=F,plot.type = 'm',ylab='')

CMplot(bovine_filt_1, col=c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
       bin.size=1e6, pch=19, band=1, cir.band=0.5, H=1.5,
       ylim=c(0,10), cex.axis=1, plot.type="b", multracks=FALSE, cex=c(0.5,1,1),
       r=0.3, ylab=expression(-log[10](italic(p))), outward=FALSE, threshold = 1e-4, threshold.lty=2, threshold.lwd=2,threshold.col='black',
       amplify= FALSE, chr.labels=NULL,highlight=highlight,highlight.text=highlight,
       signal.cex = 1.5, signal.pch = 19, signal.col="red", signal.line=1,
       cir.chr=TRUE, cir.chr.h=1.5, chr.den.col=c("darkgreen", "yellow", "red")
       , cir.legend=TRUE, cir.legend.cex=0.6, cir.legend.col="black",
       LOG10=TRUE, box=FALSE, conf.int.col="grey", file.output=TRUE,
       file="tiff", dpi=300, memo="")

#####################################
## Figure 1B
#####################################
## Percentage of different functional regions in trans-circQTLs
data<-df
data<-data[-c(12:15),]
data$Percentage<-(round(data$nums/sum(data$nums),4))*100
colnames(data)<-c('Group','nums','Percentage')
data$Group <- factor(data$Group, levels=data$Group)
mylabel<-paste(data[,3],"%")
mylabel<-rev(mylabel)
percent<-rev(data$Percentage)
p<-ggplot(data,aes(x="",y=Percentage,fill=Group)) +
   geom_bar(stat = "identity",color="white")+
   coord_polar(theta = "y") +
   theme(axis.text.x = element_blank(),
         axis.ticks = element_blank(),
         panel.grid = element_blank())
write.csv(data,'./output/circQTLs/Percentage of different functional regions in trans-circQTLs.csv')
# # pie plot
# library(ggrepel)
# colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(5)))
# p <- ggplot(data = df, mapping = aes(x = 'Content', y = nums, fill = type)) + geom_bar(stat = 'identity', position = 'stack',width = 1)
# p + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme_bw()+
   # theme(panel.grid = element_blank(),
         # rect = element_blank(),
         # axis.ticks = element_blank(),
         # axis.text = element_blank())+
   # scale_fill_manual(values = colpalettes)+
   # theme(legend.title = element_text(size=15),
         # legend.text = element_text(size=12))+
   # geom_label_repel(data = df, aes(x = 'Content', y = nums, label = type), 
                    # label.padding = unit(0.1, "lines"), 
                    # size = 8, 
                    # max.overlaps=21,
                    # fill = alpha(c("white"),0.5),
                    # show.legend = FALSE, 
                    # inherit.aes = FALSE)
#####################################
## Figure 1D
#####################################
## Beta and p value showed by voclano plot pvalue<1e-4
transcircQTL<-me$trans$eqtls
data<-transcircQTL
res<-transcircQTL[transcircQTL$pvalue<1e-4,]
# write.table(res,'hfilt_circQTL.txt',quote = F,row.names = T)
# Beta <- transcircQTL$beta
# FDR <- res$padj
res$group='not significant'
res$group[which((res$beta>=1) & (res$pvalue<1e-4))]='Positive correlation'
res$group[which((res$beta<=-1) & (res$pvalue<1e-4))]='Negative correlation'
# png(paste(output, '/',"Volcano plot.png", sep=""),width=1920,height=1920,res=300)
ggplot(res,aes(x=beta,y=-1*log10(pvalue),color=group))+geom_point(size=2)+
   geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8)+
   geom_hline(yintercept = -log10(1e-4),lty=4,col="black",lwd=0.8)+
   scale_color_manual(values = c("lightblue","grey","orangered"))+
   xlab('Beta')+ylab('-log10(pvalue)')+theme_bw()+xlim(-50,500)+ylim(0,100)+
   theme(axis.text.x = element_text(color="black", size=15,hjust = 1),
         axis.text.y = element_text(color="black", size=15),
         axis.title.x = element_text(size = 25, colour = 'black'),
         axis.title.y = element_text(size = 25, colour = 'black'),
         legend.title=element_text(size=15),
         legend.text = element_text(size=15))

# Filt the negative-positive circQTLs
high_vaild_res<-res[res$pvalue>=1e-20,]
write.table(high_vaild_res,'./output/circQTLs/hvaild_circQTLs.txt',quote = F,row.names = T)

# Filt the high vaild circQTL info
high_vaild_res<-read.table('./output/circQTLs/hvaild_circQTLs.txt',header = T)
significant<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
idx<-lapply(high_vaild_res$snps,function(x){which(significant$VariantID==x)})
high_vaild_sign<-significant[unlist(idx),]
write.table(high_vaild_sign,'./output/circQTLs/high_vaild_significant_log_quant.txt',quote = F,row.names = F)

#####################################
## Figure 1C
#####################################
## Constructing the circos data file for the interaction relationship between circQTLs and circRNAs
res<-read.table('./output/circQTLs/hvaild_circQTLs.txt',header = TRUE)
res_hvfilt<-res[,'snps'][!duplicated(res$snps)]
res_hvfilt<-data.frame(''=res_hvfilt)
prep_data <- res[,1:2]
head(snpspos)
colnames(prep_data)<-c('SNP','circRNA')
prep_data <- merge(prep_data,snpspos,by='SNP')
prep_data<-prep_data[,c(1,3,4,2)]
start=list();end=list();circRNA_chr=list()
for (i in 1:length(unlist(prep_data['circRNA']))) {
   circRNA_chr[i]<-strsplit(prep_data['circRNA'][[1]][i],':',fixed = T)[[1]][1]
   temp<-strsplit(unlist(prep_data['circRNA'][[1]][i]),':',fixed = T)[[1]][2]
   start[i]<-strsplit(temp,'-',fixed = T)[[1]][1]
   end[i]<-strsplit(temp,'-',fixed = T)[[1]][2]
   
}
circos_data<-prep_data
circos_data['circRNA_chr']<-unlist(circRNA_chr)
circos_data['start']<-unlist(start)
circos_data['end']<-unlist(end)
circos_data[,1]<-circos_data$Chromosome
circos_data$circRNA<-circos_data$Position
circos_data<-circos_data[,c(1,3,4,2,5,6,7)]
circos_data=data.frame(circos_data)
circos_data$SNP<-paste0('chr',circos_data$SNP)
circos_data$snps<-paste0(circos_data$SNP,':',circos_data$Position)
circos_data$circs<-paste0(circos_data$circRNA_chr,':',circos_data$start,'|',circos_data$end)
## make colinear map by using ggplot
library(ggforce)
ggplot() +
   geom_diagonal_wide(data=circos_data,
                      aes(snps, circs, group = Chromosome,fill=Chromosome),
                      alpha=0.5,color="black")+
   geom_diagonal_wide(data=circos_data,
                      aes(circs, snps, group = circRNA_chr),
                      alpha=0.2,color="black",
                      fill="blue")+
   theme_minimal()+
   theme(panel.grid = element_blank(),
         axis.title = element_blank(),
         axis.text = element_blank())

write.table(circos_data,'./output/circQTLs/hfilt_circQTL-1.data',quote = F,row.names = F)

## Example of the effect of rs134248710 on chr12:19479415-19482423
# significant<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
breed_info = read.csv('./config/SRR_breeds.csv')
circRNA_express<-read.csv('./output/circQTLs/log_anno_array_hfilt_3soft_anno.csv')

#####################################
## Figure S3C&F
#####################################
snp_id<-'rs109001072'
snp_id<-'rs136591890'
circ_id<-'chr25:3818938-3819240'
circ_id<-'chr9:69799938-69815551'

circRNA_express$chr_pos<-paste(circRNA_express$chrom,':',circRNA_express$start,'-',circRNA_express$end,sep='')
specific_circRNA<-circRNA_express[circRNA_express$chr_pos==circ_id,-c((length(circRNA_express)-3):length(circRNA_express))]
effect_data<-data.frame(t(specific_circRNA))
colnames(effect_data)<-'circRNA'
effect_data$SRR_list<-rownames(effect_data)
effect_data<-merge(effect_data,breed_info,by='SRR_list')
# read the freq of mutation allies in snps
specific_snps<-data.frame(t(snps$FindRow(snp_id)$row))
specific_snps$Breeds<-rownames(specific_snps)
# combine the allies and circRNA_data
total_data<-merge(effect_data,specific_snps,by='Breeds')
# Calibration of expression of circRNAs by covariate
age<-data.frame(cvrt$getSlice(1))
colnames(age)<-cvrt$columnNames
age<-data.frame('Age'=t(age))
age$Breeds<-rownames(age)
total_data<-merge(total_data,age,by='Breeds')
# Calibration
lm_res<-lm(total_data$circRNA~total_data$Age+total_data$rs109001072)
total_data$circRNA_calib<-total_data$circRNA-lm_res$coefficients[1]+total_data$Age*lm_res$coefficients[2]
total_data$Genotype<-'NA'
allies<-significant$Alleles[significant$VariantID==snp_id]
wt<-substr(allies,1,1);mu<-substr(allies,3,3)
total_data$Genotype[total_data[,snp_id]<1/3]<-paste(wt,wt,sep='')
total_data$Genotype[(total_data[,snp_id]>=1/3)&(total_data[,snp_id]<2/3)]<-paste(wt,mu,sep='')
total_data$Genotype[(total_data[,snp_id]>=2/3)]<-paste(mu,mu,sep='')
# reorder the order of genotypes
WT<-paste(wt,wt,sep='');Heter<-paste(wt,mu,sep='');MU<-paste(mu,mu,sep='')
total_data$Genotype<-factor(total_data$Genotype,levels = c(WT,Heter,MU))
# plot the boxplot of effect on circRNAs for rs134248710
ggplot(total_data)+geom_boxplot(aes(x=Genotype,y=circRNA_calib,fill=Genotype))+
   labs(title=paste(snp_id,': ',wt,'>',mu,'  ',circ_id,sep=''))+ylim(0,max(total_data$circRNA_calib)+2)+
   theme_bw()+theme(axis.text.x = element_text(color="black", size=18),
                    axis.text.y = element_text(color="black", size=18),
                    axis.title.x = element_text(size = 25, colour = 'black'),
                    axis.title.y = element_text(size = 25, colour = 'black'),
                    legend.title=element_text(size=18),
                    legend.text = element_text(size=15),
                    plot.title = element_text(hjust = 0.5))+
   annotate('text',x = 2,y=max(total_data$circRNA_calib)+1,parse=TRUE,
            label=paste("~ italic(p) ==",round(significant$pvalue[significant$VariantID==snp_id],3),
                        sep = ''),size=5)

############################################################################################################
#####################################
## Figure S5
#####################################
## Get the distribution of snps in the flanking seq of splice sites
library(ChIPseeker)
library(ggplot2)
library(zoo)
library(BayesPeak)
library(stringr)
library(GenomicRanges)

flanking = 1e3
window_size = 200
range = c(-flanking,flanking)
ss<-read.table('./reference/genome/btau9.ss',col.names = c('Chrom','start','end','strand'))
# circqtl<-read.table('significant_circQTL_log_quant.txt',header = T)
circqtl<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
# filt the splice site
circqtl_filt<-circqtl[circqtl$ConsequenceType=='intron_variant',] #'splice_region_variant'
circqtl_filt<-circqtl[circqtl$ConsequenceType=='splice_region_variant',]
# Total circQTLs or Filt the splice sites 
circqtl_filt<-subset(circqtl_filt,select = c(1,11))

# circqtl_filt$Chrom<-apply(X=data.frame(circqtl_filt$UMD3_1_1_Pos),1,FUN=function(x){strsplit(x,':',fixed=T)[[1]][1]})
# circqtl_filt$Pos<-apply(X=data.frame(circqtl_filt$UMD3_1_1_Pos),1,FUN=function(x){strsplit(x,':',fixed=T)[[1]][2]})
circqtl_filt$Chrom<-apply(X=data.frame(circqtl_filt$ARS_UCD1_2_Pos),1,FUN=function(x){strsplit(x,':',fixed=T)[[1]][1]})
circqtl_filt$Pos<-apply(X=data.frame(circqtl_filt$ARS_UCD1_2_Pos),1,FUN=function(x){strsplit(x,':',fixed=T)[[1]][2]})
circqtl_filt<-circqtl_filt[,-2]
peak<-circqtl_filt[,-1]
#5_splice_site
ss_5p<-ss[,-3]
index<-str_detect(ss$Chrom,'NKLS')
ss_5p<-ss_5p[!index,]
ss_5p$start<-ss_5p$start-flanking
ss_5p$end<-ss_5p$start+flanking
ss_5p<-GRanges(Rle(ss_5p[,1]),IRanges(start=ss_5p$start,end=ss_5p$end),Rle(strand(ss_5p$strand)))
peak<-GRanges(Rle(peak[,1]),IRanges(peak[,2]))
tagmatrix<-getTagMatrix(peak, weightCol = NULL, windows=ss_5p, flip_minor_strand = TRUE)
distri_ss5p<-colSums(tagmatrix)
#3_splice_site
ss_3p<-ss[,-3]
# index<-str_detect(ss$Chrom,'NKLS')
ss_3p<-ss_3p[!index,]
ss_3p$start<-ss_3p$start-flanking
ss_3p$end<-ss_3p$start+flanking
ss_3p<-GRanges(Rle(ss_3p[,1]),IRanges(start=ss_3p$start,end=ss_3p$end),Rle(strand(ss_3p$strand)))
tagmatrix_3p<-getTagMatrix(peak, weightCol = NULL, windows=ss_3p, flip_minor_strand = TRUE)
distri_ss3p<-colSums(tagmatrix_3p)

# Distribution of snps in bins
bins=50
xaxis<-seq(-1000,1000,2*bins)
g<-rep(1:((length(distri_ss3p)-1)/bins),each=bins)
distri_ss5p_sum<-tapply(distri_ss5p[-1],g,sum)
distri_ss3p_sum<-tapply(distri_ss3p[-1],g,sum)
distri_data_bins<-data.frame(t(rbind(c(distri_ss5p_sum[(length(distri_ss5p_sum)/2+1):length(distri_ss5p_sum)],
                                   distri_ss3p_sum[1:(length(distri_ss5p_sum)/2)]),xaxis)))
colnames(distri_data_bins)<-c('Number','Distance')
distri_data_bins$Splice_type<-'Splice_site_5p'
distri_data_bins$Splice_type[distri_data_bins$Distance<0]<-'Splice_site_3p'

ggplot(distri_data_bins)+geom_histogram(aes(x=Distance,y=Number,fill=Splice_type),stat = 'identity')+
   theme_bw()+xlab('Distance to splice sites (bp)')+ylab('Number of circQTLs')+
   theme(axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20),
         axis.text = element_text(size=15),
         legend.title = element_text(size=15),
         legend.text = element_text(size=12))+
   # annotate("text", x = -600, y = -0.5, label = "Exons",size=5)+
   # annotate("text", x = 0, y = -0.5, label = "Introns",size=5)+
   # annotate("text", x = 600, y = -0.5, label = "Exons",size=5)+
   geom_vline(xintercept=c(-250,250), linetype="dotted")
# Distribution of snps in sliding windowsize
distri_data<-data.frame(t(rbind(c(rollsum(distri_ss5p,window_size/2)[(length(rollsum(distri_ss5p,window_size/2))/2+1):length(rollsum(distri_ss5p,window_size/2))],
                                rollsum(distri_ss3p,window_size/2)[1:(length(rollsum(distri_ss5p,window_size/2))/2)]),
                                seq(-1000+window_size/2-2,1000-window_size/2,2))))
colnames(distri_data)<-c('value','Distance')
distri_data$Splice_type<-'Splice_site_5p'
distri_data$Splice_type[distri_data$Distance<0]<-'Splice_site_3p'
# distri_data<-melt(distri_data,id.vars = 'Distance',variable_name = 'Splice_type',value.name='Mean per windowsize')
distri_data$value<-scale(distri_data$value,center = T,scale = F)


ggplot(distri_data,aes(x=Distance))+geom_histogram(aes(y=value,fill=Splice_type),stat="identity")+
   theme_bw()+xlab('Distance to splice sites (bp)')+ylab('Mean per window size (200bp)')+
   theme(axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20),
         axis.text = element_text(size=15),
         legend.title = element_text(size=15),
         legend.text = element_text(size=12))+
   # annotate("text", x = -600, y = -0.5, label = "Exons",size=5)+
   # annotate("text", x = 0, y = -0.5, label = "Introns",size=5)+
   # annotate("text", x = 600, y = -0.5, label = "Exons",size=5)+
   geom_vline(xintercept=c(-250,250), linetype="dotted")
   # annotate("text", x = 350, y = -1.5, label = "Exons")+
   # annotate("text", x = 820, y = -1.5, label = "Introns")+
   # annotate("rect", xmin = -850, xmax = -250, ymin = -1.1, ymax = -1.2, alpha = .99, colour = "black")+
   # annotate("rect", xmin = -250, xmax = 0, ymin = -1.0, ymax = -1.3, alpha = .2, colour = "black")+
   # annotate("rect", xmin = 0, xmax = 700, ymin = -1.1, ymax = -1.2, alpha = .99, colour = "black")+
   # annotate("rect", xmin = 700, xmax = 950, ymin = -1.0, ymax = -1.3, alpha = .2, colour = "black")

#####################################
## Figure S9 
#####################################
##Functional enrichment analysis
BiocManager::install('org.Bt.eg.db')
BiocManager::install('createKEGGdb')
library(org.Bt.eg.db)
library(clusterProfiler)
library(createKEGGdb)
library(RColorBrewer)
total_gene<-significant$Gene
filt_gene<-significant$Gene[significant$pvalue<=1e-4]
intron_gene<-significant$Gene[significant$ConsequenceType=='intron_variant'&(significant$pvalue<=1e-4)]
ss_gene<-significant$Gene[significant$ConsequenceType=='splice_region_variant']

# total_gene<-significant$[significant$ConsequenceType=='']
eg_total = bitr(total_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Bt.eg.db")
eg_total = bitr(intron_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Bt.eg.db")
eg_total = bitr(ss_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Bt.eg.db")

results=eg_total$ENTREZID
id=na.omit(results)  #??ȡ????NA??ENTREZID

#GO enrichment analysis
ego_BP <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis

BP <- subset(ego_BP@result)
BP <- BP[1:10,]

ego_MF <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "MF",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis

MF <- subset(ego_MF@result)
MF <- MF[1:10,]
ego_CC <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "CC",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis

CC <- subset(ego_CC@result)
CC <- CC[1:10,]

# GO enrichment result of Top10
GO_result<-rbind(BP,MF,CC)
GO_result$GO_term<-rep(c('BP','MF','CC'),rep(10,3))
# 
write.table(GO_result,'./output/circQTLs/Go enrichment analysis for host gene with intron-circQTLs.txt',quote = F,row.names = F,sep = "\t")
GO_result$'-log10(qvalue)'= (-1)*log10(GO_result$qvalue)
# pdf(paste(output, '/',name," Go enrichment analysis.pdf", sep=""),width=12,height=12)
ggplot(GO_result,mapping = aes(x=reorder(Description,-log10(qvalue)),y=-log10(qvalue),fill=GO_term))+
   geom_bar(stat = "identity",position = position_dodge(width = 0.1),width = 0.5)+coord_flip()+scale_x_discrete(limits=GO_result$Description)+theme_bw()+xlab('')+
   theme(axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
         axis.text.x = element_text(size=10,angle = 80,hjust=1),axis.text.y = element_text(size=15))
tiff('C:/Users/luffy/Desktop/figures-2021.6.29/Go analysis of introns mRNA within circQTLs.tiff',width = 1500,height = 1000,res = 150)
ggplot(GO_result)+
   geom_point(aes(y=GO_term,x=Description,size=-log10(qvalue),color=GO_term))+coord_flip()+
   scale_x_discrete(limits=GO_result$Description)+theme_bw()+xlab('')+
   theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
         axis.text.x = element_text(size=12,angle = 80,hjust=1),axis.text.y = element_text(size=12))

dev.off()

kk <- enrichKEGG(gene = id,
                 organism     = 'bta',
                 pAdjustMethod = "BH",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05)
KEGG_result<-subset(kk@result)[1:10,]
write.table(KEGG_result,'./output/circQTLs/KEGG enrichment analysis for host gene with intron-circQTLs.txt',quote = F,row.names = F,sep = "\t")

KEGG_result$'-log10(qvalue)'= (-1)*log10(KEGG_result$qvalue)

# pdf(paste(output, '/',name," KEGG enrichment analysis.pdf", sep=""),width=12,height=12)
tiff('C:/Users/luffy/Desktop/figures-2021.6.29/KEGG analysis of ss mRNA within circQTLs.tiff',width = 1500,height = 1000,res = 150)
ggplot(KEGG_result)+
   geom_bar(stat = "identity",aes(x=reorder(Description,-log10(qvalue)),y=-log10(qvalue),fill=-log10(qvalue)))+coord_flip()+theme_bw()+xlab('')+
   scale_fill_distiller(palette = 'Blues',trans = 'reverse')+
   theme(axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
         axis.text.x = element_text(size=18,angle = 80,hjust=1),axis.text.y = element_text(size=15),
         legend.title = element_text(size = 15),legend.text = element_text(size=12))
dev.off()
