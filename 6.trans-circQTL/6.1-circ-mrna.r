#!/bin/Rscript

library(ggplot2)
# library(ballgown)
library(rtracklayer)
library(parallel)
library(clusterProfiler)
library(org.Bt.eg.db)
library(reshape2)

############################################################################################################
working_dir <- "path_to_your_working_dir"
all_RNA_seq_file='./SRR_list_CIRI.txt'
output_dir <- paste(working_dir,'./output/circRNA_mRNA',sep='/')
setwd(circRNAprofiler_dir)
dir.create(output_dir)
# Detect cores and initialize that in local computer
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

# get the all gtf files of mRNA expression level from CIRI2
all_gtf<-unlist(read.table(all_RNA_seq_file))
for (i in 1:length(all_gtf)) {
  gtf_each<-all_gtf[i]
  gtf_each_path<-paste0(working_dir,'/output/circRNA_identification/CIRI2/',gtf_each,'/gene/',gtf_each,'_quant_genes.list')
  gene_expr_each<-read.table(gtf_each_path,sep='\t',header = T)
  gene_expr_info<-gene_expr_each[,1:6];
  
  if (i == 1){
    gene_id<-gene_expr_info$Gene.ID
    total_gene_expr<-cbind(gene_expr_info,gene_expr_each$TPM)
    colnames(total_gene_expr)[7]<-gtf_each
    # write.table(gene_expr_info,'./combined_mrna.txt',sep='\t',quote = F,row.names = F)
  } else {
    idx<-parLapply(cl,gene_id,function(x,gene_expr_each){which(gene_expr_each$Gene.ID==x)},gene_expr_each)
    total_gene_expr<-cbind(total_gene_expr,gene_expr_each$TPM[unlist(idx)])

    sample<-7+i-1
    colnames(total_gene_expr)[sample]<-gtf_each
  }
}
# import other single-end mRNAseq TPM matrix
gene_tpm<-read.table('./output/circRNA_identification/CIRI2/gene_tpm_matrix.csv',sep=',',header = T)
head(gene_tpm)
gene_tpm$Gene.ID<-unlist(lapply(gene_tpm$gene_id,function(x){strsplit(x,'|',fixed = T)[[1]][1]}))
gene_tpm<-gene_tpm[,-1]
total_gene_expr<-merge(total_gene_expr,gene_tpm,by='Gene.ID')
total_gene_expr<-total_gene_expr[,-c(2:5)]
write.table(total_gene_expr,'./combined_mrna.txt',sep='\t',quote = F,row.names = F)
stopCluster(cl)

# Filter the result of rowsum>=1
total_gene_expr<-read.table('./output/circRNA_identification/combined_mrna.txt',sep='\t',header = T)
total_gene_expr_filt<-total_gene_expr[rowSums(total_gene_expr[,7:length(total_gene_expr)])>=1,]

total_gene_expr_mean<-cbind(total_gene_expr_filt[,1:6],
                            rowMeans(total_gene_expr_filt[,7:length(total_gene_expr_filt)]))
breed_info = read.csv('./SRR_breeds.csv')
breeds = breed_info$Breeds[!duplicated(breed_info$Breeds)]
total_gene_convert=data.frame(total_gene_expr_filt[,1:6])
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  total_gene_expr_breed = subset(total_gene_expr_filt,select = srr)
  total_gene_convert[,as.character(breeds[i])]=apply(total_gene_expr_breed,1,mean)
}
# read the circQTLs of intersection or total significant
total<-read.table('./output/circQTLs/high_vaild_significant_log_quant.txt',header = T)
intersect<-read.table('./output/coloc/intersect.txt',header=T,sep='\t')
idx<-lapply(intersect$SYMBOL, 
            function(x,gene_expr_each){which(gene_expr_each$Gene.Name==x)},total_gene_convert)
filt_gene_convert<-total_gene_convert[unique(unlist(idx)),]
# total
idx_total<-lapply(total$Gene, 
            function(x,gene_expr_each){which(gene_expr_each$Gene.Name==x)},total_gene_convert)
filt_gene_convert<-total_gene_convert[unique(unlist(idx_total)),]

#####################################
## Figure S6
#####################################
## Figure S6A Relationship of mRNAs among different breeds
## relationship of these intersection genes (gwas-circqtls) among different breeds (NO adjusting)
# relationship of mean among breeds
library(pheatmap)
library(RColorBrewer)
data<-filt_gene_convert[,7:length(filt_gene_convert)]
pcr_matrix_mean<-cor(data,data,method = 'pearson')
cols<-colnames(pcr_matrix_mean)
cols<-data.frame(cols)
colnames(cols)<-'Breeds'
rownames(cols)<-cols$Breeds
pheatmap(pcr_matrix_mean,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         cluster_rows = T,cluster_cols = T,show_colnames = F,legend=T,annotation_col = cols,
         annotation_row = cols,fontsize = 15,
         show_rownames=F)

## Figure S6B
# relationship of measure among samples
data1<-total_gene_expr_filt[unlist(idx),7:length(total_gene_expr_filt)]
data1<-total_gene_expr_filt[unlist(idx_total),7:length(total_gene_expr_filt)]
pcr_matrix<-cor(data1,data1,method = 'pearson')
cols<-colnames(pcr_matrix)
cols<-data.frame(cols)
colnames(cols)<-'SRR_list'
annotation_col<-merge(cols,breed_info,by.y ='SRR_list')
rownames(annotation_col)<-annotation_col$SRR_list
annotation_col<-subset(annotation_col,select = 'Breeds')
pheatmap(pcr_matrix,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         cluster_rows = T,cluster_cols = T,show_colnames = F,legend=T,
         annotation_col=annotation_col,annotation_legend = T,
         annotation_names_col = T,annotation_row=annotation_col,
         show_rownames=F,fontsize = 15)

#####################################
## Figure S7
#####################################
# relationship of mRNA and circRNA among breeds for each pair
# read the pairs
hvaild_circqtl<-read.table('./output/circQTLs/hvaild_circQTLs.txt',header = T)
# read the expression of circRNAs
circ_expr<-read.table('./output/circQTLs/circRNA.txt',header = T,sep='\t')
colnames(intersect)[4]<-'snps'
colnames(intersect)[1]<-'Gene'
colnames(total)[1]<-'snps'
intersect_pairs<-merge(intersect,hvaild_circqtl,by='snps')
intersect_pairs<-merge(total,hvaild_circqtl,by='snps')

snps<-unique(intersect_pairs$snps)

pcr_data<-intersect_pairs[,c('snps','gene')]
# effect of circQTL on circRNAs and mRNA
library(ComplexHeatmap)
library(circlize)
# create circQTLoncomap and mRNA/circRNA heatmap
circ_mrna_qtldata<-intersect_pairs[,c('snps','gene','pvalue.x','ConsequenceType','Chrom','ARS_UCD1_2_Pos')]
circ_mrna_qtldata$mRNA<-intersect_pairs$Gene

circ_mrna_qtldata$ARS_UCD1_2_Pos<-unlist(lapply(circ_mrna_qtldata$ARS_UCD1_2_Pos,function(x){
  strsplit(x,split = ':',fixed = TRUE)[[1]][2]}))
cytoband = read.table("./reference/cytoBed/cytoBandIdeo.txt", sep = "\t")
cytoband = cytoband[grepl('chr[0-9,X,M](\\d)?',cytoband$V1),]
cytoband$V5[seq(1,length(cytoband$V1))] = 'gpos25'
cytoband$V5[seq(2,length(cytoband$V1),3)] = 'gpos75'
cytoband$V5[seq(3,length(cytoband$V1),3)] = 'gpos100'
# circos.initializeWithIdeogram(cytoband,chromosome.index = paste0("chr", c(1:29, "X")),plotType = c('axis','labels',"cytobands"))
circos.genomicInitialize(cytoband)
# add chromIdenogram track
circos.genomicTrackPlotRegion(cytoband, ylim = c(0, 1), bg.border = NA, track.height = 0.02,
                              panel.fun = function(region, value, ...) {
                                col = cytoband.col(value[[2]])
                                circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = col, border = NA, ...)
                                xlim = get.cell.meta.data("xlim")
                                circos.rect(xlim[1], 0, xlim[2], 1, border = "black")
                              }, cell.padding = c(0, 0, 0, 0)
)

# add circQTL track
circqtl_bed<-circ_mrna_qtldata[,c('Chrom','ARS_UCD1_2_Pos','ARS_UCD1_2_Pos','pvalue.x')]
circqtl_bed$Chrom<-paste0('chr',circqtl_bed$Chrom)
circqtl_bed$pvalue.x=-log10(circqtl_bed$pvalue.x)
colnames(circqtl_bed)<-c('chr','start','end','value1')
circqtl_bed$start<-as.numeric(circqtl_bed$start)
circqtl_bed$end<-as.numeric(circqtl_bed$end)
# circqtl_bed$end<-as.numeric(circqtl_bed$start)+100
# bed = generateRandomBed(nr = 200)
circos.genomicTrack(circqtl_bed, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16,cex = 0.5,col = 1, ...)}
                    )
##
circqtl_bed<-circqtl_bed[!duplicated(circqtl_bed),]
# circqtl_bed$value2<-1
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
circos.initializeWithIdeogram(cytoband[cytoband$V1!="chrM",])
circos.genomicHeatmap(circqtl_bed, side = "inside",col = col_fun)
circos.genomicIdeogram()
##

# relationship between circRNA and mRNAs
for (i in 1:length(snps)) {
  # i=
  circ_name<-intersect_pairs[intersect_pairs$snps==snps[i],'gene']
  # a<-log2(filt_gene_convert[filt_gene_convert$Gene.Name==intersect_pairs$SYMBOL[i],7:length(filt_gene_convert)]+1)
  mrna_name<-intersect_pairs$Gene[intersect_pairs$snps==snps[i]]
  a<-log2(filt_gene_convert[filt_gene_convert$Gene.Name==intersect_pairs$Gene[intersect_pairs$snps==snps[i]],7:length(filt_gene_convert)]+1)
  
  if (nrow(a)==0){
    print('opps')
    print(a)
    next
  } else if (nrow(a)>=1) {
    a<-a[1,]
  }
  for (j in 1:length(circ_name)) {
    # pcr<-list()
    b<-circ_expr[circ_expr$circ_name==circ_name[j],2:length(circ_expr)]
    res<-cor.test(unlist(a),unlist(b),method = 'pearson')
    pcr_data$PCR[pcr_data$snps==snps[i]&pcr_data$gene==circ_name[j]]<-res$estimate
    pcr_data$pvalue[pcr_data$snps==snps[i]&pcr_data$gene==circ_name[j]]<-res$p.value
    # pcr_data$PCR[pcr_data$snps==snps[i]&pcr_data$gene==circ_name[j]]<-c
    pcr_data$mRNA[pcr_data$snps==snps[i]&pcr_data$gene==circ_name[j]]<-mrna_name
    # pcr<-append(pcr,c)
    
  }
  # pcr<-unlist(pcr)
  # d<-pcr[which(abs(pcr)==max(abs(pcr)))]
  
  
}
# pcr_data<-pcr_data[!is.na(pcr_data$PCR)&abs(pcr_data$PCR)>=0.6,]
pcr_data<-pcr_data[!is.na(pcr_data$PCR),]
pcr_data$type<-'Negative';pcr_data$type[pcr_data$PCR>0]<-'Positive'
# Annotate mRNA genes by KEGG pathway
up_genes<-pcr_data$mRNA[pcr_data$type=='Positive'];down_genes<-pcr_data$mRNA[pcr_data$type=='Negative']

gene.df <- bitr(pcr_data$mRNA[abs(pcr_data$PCR>0.5)], fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Bt.eg.db)
go_res<-enrichGO(gene.df$ENTREZID, 'org.Bt.eg.db', ont="BP", pvalueCutoff=0.01)
go_res_filt<-go_res@result[go_res@result$pvalue<=0.05,]
anno_gene<-function(x){
  # x is a char
  library(stringr)
  idx<-!is.na(str_match(go_res_filt$geneID,x))
  tmp<-go_res_filt[idx,]
  output<-tmp$Description[tmp$pvalue==min(tmp$pvalue)][1]
  return(output)
}
colnames(gene.df)[1]<-'mRNA'
pcr_data$Go_term<-unlist(lapply(pcr_data$ENTREZID, anno_gene))
write.table(pcr_data,'./output/circRNA_mRNA/mrna_circ.txt',quote = F,row.names = F,sep='\t')

pcr_data<-read.table('./output/circRNA_mRNA/mrna_circ.txt',header = T,sep='\t')
## Figure S7B-C
ggplot(pcr_data)+geom_violin(aes(x=type,y=PCR,fill=type),alpha=0.7)+theme_bw()+
  # geom_jitter(aes(x=type,y=PCR,color=type),width = 0.2,alpha=0.1)+
  # scale_color_manual(values = c("darkblue","darkred"))+
  scale_fill_manual(values = c("darkblue","darkred"))+
  xlab('Beta type')+ylab('Pearson correlation')+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))

## Calibration of circRNA expression
# genotype
library(GenVisR)
snps_geno<-read.table('./output/circQTLs/snp.txt',header = T)

# cvrt
age<-read.table('./output/matrixeqtl/cov.txt',header = T)
age<-data.frame(t(age)[-1,])
colnames(age)<-'Age';age$Breeds<-rownames(age);age$Age<-as.numeric(age$Age)
# circRNA expression matrix
# significant<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
breed_info = read.csv('./SRR_breeds.CSV')
circRNA_express<-read.csv('./output/circQTLs/log_anno_array_hfilt_3soft_anno.csv')
# total
freq_grad<-names(sort(colMeans(snps_geno[,2:length(snps_geno)])))
idx<-unlist(lapply(freq_grad, function(x){which(breed_info$Breeds==x)}))
srr<-breed_info[idx,]
circRNA_expr_filt<-data.frame()
for (i in 1:length(snps)){
  # i=1
  snp_id<-snps[i]
  # read the freq of mutation allies in snps
  specific_snps<-data.frame(snp_id=t(snps_geno[snps_geno$SNP_ID==snp_id,2:length(snps_geno)]))
  colnames(specific_snps)<-"snp_id";specific_snps$Breeds<-rownames(specific_snps);
  circ_name<-intersect_pairs[intersect_pairs$snps==snps[i],'gene']
  # mrna_name<-intersect_pairs$Gene[intersect_pairs$snps==snps[i]]
  
  for (j in 1:length(circ_name)){
    # j=1
    circ_id<-circ_name[j]
    # get the related circRNA expression
    circRNA_express$chr_pos<-paste(circRNA_express$chrom,':',circRNA_express$start,'-',circRNA_express$end,sep='')
    specific_circRNA<-circRNA_express[circRNA_express$chr_pos==circ_id,-c((length(circRNA_express)-3):length(circRNA_express))]
    effect_data<-data.frame(t(specific_circRNA))
    colnames(effect_data)<-'circRNA'
    effect_data$SRR_list<-rownames(effect_data)
    effect_data<-merge(effect_data,breed_info,by='SRR_list')
    # combine the allies and circRNA_data
    total_data<-merge(effect_data,specific_snps,by='Breeds')
    total_data<-merge(total_data,age,by='Breeds')
    # LM model
    lm_res<-lm(total_data$circRNA~total_data$Age+total_data$snp_id)
    total_data$circRNA_calib<-total_data$circRNA-lm_res$coefficients[1]+total_data$Age*lm_res$coefficients[2]
    total_data$Genotype<-'NA'
 
    total_data$Genotype[total_data[,"snp_id"]<1/3]<-'WT'
    total_data$Genotype[(total_data[,"snp_id"]>=1/3)&
                          (total_data[,"snp_id"]<2/3)]<-'HT'
    total_data$Genotype[(total_data[,"snp_id"]>=2/3)]<-'MU'
    # reorder the order of genotypes
    # WT<-paste(wt,wt,sep='');Heter<-paste(wt,mu,sep='');MU<-paste(mu,mu,sep='')
    # total_data$Genotype<-factor(total_data$Genotype,levels = c(WT,Heter,MU))
    # total_data$Genotype<-factor(total_data$Genotype,levels = c("WT","HT","MU"))
    circRNA_expr_filt[circ_name,total_data$SRR_list]<-total_data$circRNA
  }
    
}
write.table(circRNA_expr_filt,'./circRNA_mRNA/circRNA_calibrated_expr.txt',quote = F,row.names = F,sep='\t')

## Figure S7A
## circosplot of circRNA and mRNA
combine<-function(a,b,col1,col2){
  colnames(b)[col2]<-colnames(a)[col1]
  return(merge(a,b,by=colnames(a)[col1]))
}
circRNA_expr_heatmap<-circRNA_expr_filt
circRNA_expr_heatmap$circ<-rownames(circRNA_expr_heatmap)
circRNA_expr_heatmap$chrom<-unlist(lapply(circRNA_expr_heatmap$circ,function(x){strsplit(x,':',fixed = T)[[1]][1]}))
circRNA_expr_heatmap<-combine(circRNA_expr_heatmap,circ_mrna_qtldata[,c('gene','mRNA')],93,1)
circRNA_expr_heatmap<-circRNA_expr_heatmap[!duplicated(circRNA_expr_heatmap$circ),]

#mRNA
mRNA_expr_heatmap<-total_gene_expr_filt[,7:length(total_gene_expr_filt)]
mRNA_expr_heatmap$mRNA<-total_gene_expr_filt$Gene.Name
mRNA_expr_heatmap$circ_mRNA<-''
idx<-unlist(lapply(circRNA_expr_heatmap$mRNA, function(x){
  y<-strsplit(x,',',fixed = T)[[1]]
  for (i in 1:length(y)){
    idx<-c(which(mRNA_expr_heatmap$mRNA==y[i]))
    if (length(idx)!=0){names(idx)<-x;break}
  }
  
  return(idx)}))
mRNA_expr_heatmap$circ_mRNA[idx]<-names(idx)
library(ComplexHeatmap)
## circos plot in samples level
mRNA_expr_heatmap<-combine(mRNA_expr_heatmap,circRNA_expr_heatmap[,c('chrom','mRNA','circ')],97,2)
circRNA_expr_heatmap_filt<-combine(circRNA_expr_heatmap,mRNA_expr_heatmap[,c('chrom','circ_mRNA')],95,2)
circRNA_expr_heatmap_filt<-circRNA_expr_heatmap_filt[!duplicated(circRNA_expr_heatmap_filt$circ),]
# tapply(circRNA_expr_heatmap_filt[,srr$SRR_list],INDEX=circRNA_expr_heatmap_filt$mRNA,FUN=mean)
# rownames(circRNA_expr_heatmap_filt)<- paste0("R", 1:length(circRNA_expr_heatmap_filt$mRNA))
# colnames(circRNA_expr_heatmap_filt)<- paste0("C", 1:length(circRNA_expr_heatmap_filt))
circos.clear()
tiff(filename = './figures/circosplot.tiff',width = 1800,height = 1600,res = 300)
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
library(gridBase)
circos.initializeWithIdeogram(cytoband[cytoband$V1!="chrM",])
# circos.genomicHeatmap(circqtl_bed, side = "inside",col = col_fun)

# par(omi = gridOMI(), new = TRUE)
circRNA_expr_heatmap_filt[,srr$SRR_list]<-scale(log2(circRNA_expr_heatmap_filt[,srr$SRR_list]+1),center = T,scale = T)
col_fun1 = colorRamp2(c(min(circRNA_expr_heatmap_filt[,srr$SRR_list]), median(apply(circRNA_expr_heatmap_filt[,srr$SRR_list],FUN=median,MARGIN = 1)),
                        max(circRNA_expr_heatmap_filt[,srr$SRR_list])),
                      c("#337abf", "#fafcfe", "#932523"))
# circos.heatmap(mat1, split = split, col = col_fun1)
split =factor(circRNA_expr_heatmap_filt$chrom.x,levels=unique(circRNA_expr_heatmap_filt$chrom.x))
circRNA_expr_heatmap_filt$chr<-unlist(lapply(circRNA_expr_heatmap_filt$circ,function(x){strsplit(x,':',fixed = T)[[1]][1]}))
circRNA_expr_heatmap_filt$start<-unlist(lapply(circRNA_expr_heatmap_filt$circ,function(x){temp<-strsplit(x,':',fixed = T)[[1]][2];
  return(as.numeric(strsplit(temp,'-')[[1]][1]))}))
circRNA_expr_heatmap_filt$end<-unlist(lapply(circRNA_expr_heatmap_filt$circ,function(x){temp<-strsplit(x,':',fixed = T)[[1]][2];
return(as.numeric(strsplit(temp,'-')[[1]][2]))}))
circRNA_expr_heatmap_filt<-circRNA_expr_heatmap_filt[circRNA_expr_heatmap_filt$chr%in%cytoband$V1,]

circos.genomicHeatmap(circRNA_expr_heatmap_filt[,c("chr",'start',"end",srr$SRR_list)],col=col_fun1,side = 'inside',
                      heatmap_height=0.2,
                      connection_height = mm_h(2))

# format mRNA expression matrix
mRNA_expr_heatmap$chr<-unlist(lapply(mRNA_expr_heatmap$circ,function(x){strsplit(x,':',fixed = T)[[1]][1]}))
mRNA_expr_heatmap$start<-unlist(lapply(mRNA_expr_heatmap$circ,function(x){temp<-strsplit(x,':',fixed = T)[[1]][2];
return(as.numeric(strsplit(temp,'-')[[1]][2]))}))
mRNA_expr_heatmap$end<-unlist(lapply(mRNA_expr_heatmap$circ,function(x){temp<-strsplit(x,':',fixed = T)[[1]][2];
return(as.numeric(strsplit(temp,'-')[[1]][2]))}))

mRNA_expr_heatmap[,srr$SRR_list]<-scale(log2(mRNA_expr_heatmap[,srr$SRR_list]+1),center = T,scale = T)
dim(mRNA_expr_heatmap)

col_fun = colorRamp2(c(min((mRNA_expr_heatmap[,srr$SRR_list])), median(apply(mRNA_expr_heatmap[,srr$SRR_list],MARGIN = 1,FUN = median)),
                       max((mRNA_expr_heatmap[,srr$SRR_list]))), c("#91bfdb", "#ffffbf","#fc8d59"))

mRNA_expr_heatmap<-mRNA_expr_heatmap[mRNA_expr_heatmap$chr%in%cytoband$V1,]
circos.genomicHeatmap(mRNA_expr_heatmap[,c("chr",'start',"end",srr$SRR_list)],heatmap_height=0.2,
                      col=col_fun,side = 'inside',connection_height = NULL)

# calculate the PCR
pcr_data<-circRNA_expr_heatmap_filt[,c('mRNA','chrom.x')]
for (i in 1:length(circRNA_expr_heatmap_filt$mRNA)){
  res<-cor.test(unlist(circRNA_expr_heatmap_filt[i,srr$SRR_list]),unlist(mRNA_expr_heatmap[i,srr$SRR_list]))
  
  pcr_data$chr[i]<-circRNA_expr_heatmap_filt$chr[i]
  pcr_data$start[i]<-circRNA_expr_heatmap_filt$start[i]
  pcr_data$end[i]<-circRNA_expr_heatmap_filt$end[i]
  pcr_data$PCR[i]<-res$estimate
  pcr_data$p[i]<-res$p.value
  
}
# pcr_data$res<-pcr_data$PCR;pcr_data$res[pcr_data$p>0.05]<-0.4
# add PCR track
circos.initializeWithIdeogram(cytoband[cytoband$V1!="chrM",])
pcr_data<-pcr_data[,-c(1:2)]
# region = circRNA_expr_heatmap_filt[,c("chr",'start',"end")]
circos.genomicTrackPlotRegion(pcr_data,  ylim = c(-1, 1),numeric.column = c(4,5),
                              panel.fun = function(region, value,...) {

  y = value[[1]]
  circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")

  circos.genomicPoints(region, y, pch=16,cex=0.5,col = ifelse(abs(y) > 0.5, ifelse(y>0.5,"red","blue"), "grey"))
}, cell.padding = c(0.02, 0, 0.02, 0))

# add the legend of colorbars
lgd_mRNAexpr = Legend(title = "mRNA", col_fun = col_fun)
lgd_expr = Legend(title = "circRNA", col_fun = col_fun1)
circos.clear()
# lgd_PCR = Legend(title = "PCR", col_fun = col_fun)
h=dev.size()[2]
lgd_list = packLegend(lgd_expr, lgd_mRNAexpr,max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = unit(1.1, "snpc"), just = "right")
dev.off()

############################
## Figure S8
############################
## circos plot in breeds level
# # circos.heatmap(mat1, split = split, col = col_fun1)
circRNA_expr_heatmap_breeds<-circRNA_expr_heatmap_filt[,c("chr",'start',"end")]
mRNA_expr_heatmap_breeds<-mRNA_expr_heatmap[,c("chr",'start',"end")]
for (i in 1:length(breeds)){
  idx = srr$SRR_list[srr$Breeds==breeds[i]]
  print(idx)
  circRNA_expr_heatmap_breeds[,breeds[i]]<-rowMeans(circRNA_expr_heatmap_filt[,idx])
  mRNA_expr_heatmap_breeds[,breeds[i]]<-rowMeans(mRNA_expr_heatmap[,idx])
}
# library(gridBase)
circos.clear()
circos.initializeWithIdeogram(cytoband[cytoband$V1!="chrM",])
# scale
# circRNA_expr_heatmap_breeds[,breeds]<-scale(circRNA_expr_heatmap_breeds[,breeds],center = T,scale = T)
col_fun1 = colorRamp2(c(min(circRNA_expr_heatmap_breeds[,breeds]), median(apply(circRNA_expr_heatmap_breeds[,breeds], 1, FUN = median)),
                        max(circRNA_expr_heatmap_breeds[,breeds])),
                      c("#337abf", "#fafcfe", "#932523"))
# split =factor(circRNA_expr_heatmap_breeds$chrom.x,levels=unique(circRNA_expr_heatmap_breeds$chrom.x))

#
tiff(filename = 'circosplot-2.tiff',width = 1800,height = 1600,res = 300)
circle_size = unit(1, "snpc") # snpc unit gives you a square region
# pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
#                       just = c("left", "center")))
# par(omi = gridOMI(), new = TRUE)
# circRNA_expr_heatmap_breeds<-data.frame(circRNA_expr_heatmap_breeds)
circos.genomicHeatmap(circRNA_expr_heatmap_breeds[,c("chr",'start',"end",breeds)],col=col_fun1,side = 'inside',
                      heatmap_height=0.2,
                      connection_height = mm_h(2))
# circos.heatmap(as.matrix(circRNA_expr_heatmap_breeds[,breeds]),show.sector.labels = F,
               # split = split,col=col_fun1)
# max(scale(mRNA_expr_heatmap[,srr$SRR_list]))), c("green", "white", "red"))
# mRNA_expr_heatmap_breeds<-scale(log2(mRNA_expr_heatmap_breeds[,breeds]+1),center = T,scale = T)
col_fun = colorRamp2(c(min((mRNA_expr_heatmap_breeds[,breeds])), median(apply(mRNA_expr_heatmap_breeds[,breeds], 1, FUN = median)),
                       max((mRNA_expr_heatmap_breeds[,breeds]))), c("#91bfdb", "#ffffbf","#fc8d59"))
# split=factor(mRNA_expr_heatmap_breeds$chrom,levels=unique(mRNA_expr_heatmap_breeds$chrom))
circos.genomicHeatmap(mRNA_expr_heatmap_breeds[,c("chr",'start',"end",breeds)],col=col_fun,side = 'inside',
                      heatmap_height=0.2,
                      connection_height = NULL)
# circos.heatmap(mRNA_expr_heatmap_breeds[,breeds],
#                col=col_fun)
# pcr_data<-cor(mRNA_expr_heatmap_breeds[,breeds],circRNA_expr_heatmap_breeds[,breeds],)
pcr_data<-circRNA_expr_heatmap_breeds[,c("chr",'start',"end")]
for (i in 1:length(circRNA_expr_heatmap_breeds$chr)){
  # pcr_data$PCR[i]<-cor(unlist(circRNA_expr_heatmap_breeds[i,breeds]),mRNA_expr_heatmap_breeds[i,breeds])
  res<-cor.test(unlist(circRNA_expr_heatmap_breeds[i,breeds]),unlist(mRNA_expr_heatmap_breeds[i,breeds]))
  
  pcr_data$chr[i]<-circRNA_expr_heatmap_breeds$chr[i]
  pcr_data$start[i]<-circRNA_expr_heatmap_breeds$start[i]
  pcr_data$end[i]<-circRNA_expr_heatmap_breeds$end[i]
  pcr_data$PCR[i]<-res$estimate
  pcr_data$p[i]<-res$p.value
}

##
circos.genomicTrackPlotRegion(pcr_data,  ylim = c(-1, 1),numeric.column = c(4,5),
  panel.fun = function(region, value,...) {
    
    y = value[[1]]
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
    
    circos.genomicPoints(region, y, pch=16,cex=0.5,col = ifelse(abs(y) > 0.5, ifelse(y>0.5,"red","blue"), "grey"))
  
  }, cell.padding = c(0.02, 0, 0.02, 0))

lgd_expr = Legend(title = "circRNA", col_fun = col_fun1)
lgd_mRNAexpr = Legend(title = "mRNA", col_fun = col_fun)
# lgd_expr = Legend(title = "circRNA", col_fun = col_fun1)

# lgd_PCR = Legend(title = "PCR", col_fun = col_fun)
h=dev.size()[2]
lgd_list = packLegend(lgd_expr, lgd_mRNAexpr,max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = unit(1.1, "snpc"), just = "left")
dev.off()

# # boxplot+ violin plot
# circRNA_expr_filt_1<-data.frame(row.names = rownames(circRNA_expr_filt))
# for (i in 1:length(freq_grad)){
  # # i=1
  # circRNA_expr_filt_1[,freq_grad[i]]<-rowMeans(circRNA_expr_filt[,srr$SRR_list[srr$Breeds==freq_grad[i]]])
# }

# circRNA_expr_filt_1<-melt(circRNA_expr_filt_1,variable.name = 'Breeds',value.name = 'circRNA')
# circRNA_expr_filt_line<-data.frame('Breeds'=unique(circRNA_expr_filt_1$Breeds))
# circRNA_expr_filt_line$Median<-unlist(lapply(unique(circRNA_expr_filt_1$Breeds),
                                    # function(x){median(circRNA_expr_filt_1[circRNA_expr_filt_1$Breeds==x,2])}));
# circRNA_expr_filt_line$sd<-unlist(lapply(unique(circRNA_expr_filt_1$Breeds),
                                  # function(x){sd(circRNA_expr_filt_1[circRNA_expr_filt_1$Breeds==x,2])}));
# # ANOVA test
# res<-aov(circRNA~Breeds,circRNA_expr_filt_1)
# pvalue<-format(summary(res)[[1]]['Breeds','Pr(>F)'],digits = 3)
# ggplot()+geom_violin(data=circRNA_expr_filt_1,aes(x=Breeds,y=circRNA,fill=Breeds))+
  # geom_boxplot(data=circRNA_expr_filt_1,aes(x=Breeds,y=circRNA,fill=Breeds),width=0.2,position=position_dodge(0.9))+
  # geom_line(data=circRNA_expr_filt_line,aes(x=Breeds,y=Median,group=1),linetype='dashed',size=1)+
  # geom_point(data=circRNA_expr_filt_line,
    # aes(x=Breeds,y=Median))+
  # ylim(0,17)+ylab('Clabriation experission \nof circRNAs')+
  # theme_bw()+theme(axis.title.x = element_text(size=18),
                   # axis.title.y = element_text(size=18),
                   # axis.text.x = element_text(angle=45,hjust=1,size=18),
                   # axis.text.y = element_text(size=15),
                   # legend.title = element_text(size=15),
                   # legend.text = element_text(size=12))+
  # annotate('text',x=3.5,y=16.5,parse=TRUE,size=7,
           # label=paste("~italic(P)==",pvalue,sep=''))
# ggplot(circRNA_expr_filt_line,aes(x=Breeds,y=Mean))+geom_line(aes(group = 1))+geom_point()+
  # geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(0.05))+
  # theme_bw()

################################################################################################
## Fig 15A-B
## Distribution of ABS isoforms among breeds
total_filt_abs_reads<-read.table('./ABS_event/ABS_events_readscount.txt',sep='\t',header = T)
# hfilt
total_filt_abs_reads<-read.table('./ABS_event/ABS_events_readscount_hfilt.txt',sep='\t',header = T)
rownames(total_filt_abs_reads)<-paste(total_filt_abs_reads$chrom,':',total_filt_abs_reads$start,'-',total_filt_abs_reads$end,sep='')

# introns
intron_snps<-read.table('./output/circQTLs/high_vaild_significant_log_quant.txt',header = T)
intron_snps<-intron_snps[intron_snps$ConsequenceType=='intron_variant',]
snps<-intersect(intron_snps$VariantID,hvaild_circqtl$snps)
idx<-unlist(lapply(snps, function(x){which(hvaild_circqtl$snps==x)}))
circ_intron<-hvaild_circqtl$gene[idx]
idx_intron<-unlist(lapply(circ_intron, function(x){which(rownames(total_filt_abs_reads)==x)}))
# total_filt_info<-total_filt_info[unique(idx_intron),]
total_filt_abs_reads<-total_filt_abs_reads[unique(idx_intron),]
library(stringr)

# Calculate the expression of two isoforms (longer and shorter) for each circRNA in breeds
abs_class<-na.omit(unique(total_filt_abs_reads$ABS_class[total_filt_abs_reads$ABS!='A5BS & A3BS'&
                                                           is.na(str_match(total_filt_abs_reads$ABS_class,'NA'))]))

# abs_class<-na.omit(unique(total_filt_abs_reads$ABS_class))
abs_isforms_expr<-data.frame('ABS_class'=rep(abs_class,each=2))
col=colnames(total_filt_abs_reads)[6:length(total_filt_abs_reads)]
breeds_info<-read.table('./output/matrixeqtl/SRR_breeds.CSV',header = T,sep=',')
srr<-breeds_info$SRR_list
total_filt_abs_reads$mean_expr<-rowMeans(total_filt_abs_reads[,srr])
##
for (i in 1:length(abs_class)) {
  # i=1
  abs_isforms_expr[(2*i-1),'type']='predominant'
  
  idx<-!is.na(str_match(total_filt_abs_reads$ABS_class,abs_class[i]))
  tmp<-total_filt_abs_reads[idx,]
  abs_isforms_expr[(2*i-1),col]=tmp[which(tmp$mean_expr==max(tmp$mean_expr)),col]
  abs_isforms_expr[2*i,'type']='others'
  abs_isforms_expr[2*i,col]=colMeans(tmp[which(tmp$mean_expr<max(tmp$mean_expr)),col])
  
}
##
for (i in 1:length(abs_class)) {
  # i=2
  abs_isforms_expr[(2*i-1),'type']='shorter'
  idx<-!is.na(str_match(total_filt_abs_reads$ABS_class,abs_class[i]))
  tmp<-total_filt_abs_reads[idx,]
  abs_isforms_expr[(2*i-1),col]=tmp[which(tmp$end-tmp$start==min(tmp$end-tmp$start)),col]
  abs_isforms_expr[2*i,'type']='longer'
  abs_isforms_expr[2*i,col]=colMeans(tmp[which(tmp$end-tmp$start>min(tmp$end-tmp$start)),col])
  
}

# import breed info and combine the expr info according to breeds
breed_info<-read.table('./output/matrixeqtl/SRR_breeds.CSV',header = T,sep=',')
breeds<-unique(breed_info$Breeds)
abs_isforms_bred<-abs_isforms_expr[,1:2]
for (i in 1:length(breeds)) {
  # i=1
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  abs_isforms_bred[,breeds[i]]<-rowMeans(abs_isforms_expr[,srr],na.rm = T)
}
abs_isforms_bred<-melt(abs_isforms_bred,id.vars = c('ABS_class','type'),
                       variable.name = 'Breeds',value.name = 'MeanExpr')
# abs_isforms_bred<-abs_isforms_bred[!is.na(abs_isforms_bred$type)]
abs_isforms_bred$type<-factor(abs_isforms_bred$type,levels=c('shorter','longer'))
# set isoform color
isoform_color<-pal_jama()(10)[c(3,4)]
names(isoform_color)<-c('longer','shorter')
# plot the boxplot
p<-ggplot(abs_isforms_bred)+geom_boxplot(aes(x=Breeds,y=MeanExpr,fill=type))+
  ylab('Mean expression of circRNAs')+
  scale_fill_manual(values = isoform_color)+
  theme_bw()+theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text.x = element_text(angle=45,hjust=1,size=18),
                   axis.text.y = element_text(size=18),
                   legend.title = element_text(size=18),
                   legend.text = element_text(size=15))

# Student's t test
p_values<-list()
for (i in 1:length(breeds)) {
  # i=1
  isoform_each<-abs_isforms_bred[abs_isforms_bred$Breeds==breeds[i],]
  p_values[[i]]<-format(t.test(MeanExpr~type,isoform_each,alternative = 'greater')$p.value,digits = 3)
  p<-p+annotate('text',x=i,y=16.5,parse=TRUE,size=5,
               label=paste("~italic(P)==",p_values[[i]],sep=''))+
    annotate('segment', x=i-0.2, xend=i+0.2, y=16.0, yend=16.0,color='black',
             arrow=arrow(angle = 0,length = unit(0,'cm')))+
    annotate('segment', x=i-0.2, xend=i-0.2, y=15.5, yend=16.0,color='black',
             arrow=arrow(angle = 0,length = unit(0,'cm')))+
    annotate('segment', x=i+0.2, xend=i+0.2, y=15.5, yend=16.0,color='black',
             arrow=arrow(angle = 0,length = unit(0,'cm')))
}
print(p)

## Correlation of PCU (Predominant circRNA unit) among breeds
abs_isform_cor<-data.frame('ABS_class'=abs_class);col<-colnames(abs_isforms_expr)[3:length(abs_isforms_expr)]
for (i in 1:length(abs_class)) {
  # i=1
  abs_isform_cor[i,col]<-abs_isforms_expr[abs_isforms_expr$ABS_class==abs_class[i]&abs_isforms_expr$type=='shorter',3:length(abs_isforms_expr)]/
    (abs_isforms_expr[abs_isforms_expr$ABS_class==abs_class[i]&abs_isforms_expr$type=='shorter',3:length(abs_isforms_expr)]+
       abs_isforms_expr[abs_isforms_expr$ABS_class==abs_class[i]&abs_isforms_expr$type=='longer',3:length(abs_isforms_expr)])
  
}

# Pearson correlation calcualte
abs_isform_cor[is.na(abs_isform_cor)]<-0
cor_mat<-cor(abs_isform_cor[,-1],abs_isform_cor[,-1],method = 'pearson')
breed_info = readxl::read_xlsx('./output/matrixeqtl/SRR_breeds.xlsx')
cols<-colnames(cor_mat)
cols<-data.frame(cols)
colnames(cols)<-'SRR_list'
annotation_col<-merge(cols,breed_info,by.y ='SRR_list');
rownames(annotation_col)<-annotation_col$SRR_list
annotation_col<-subset(annotation_col,select = 'Breeds');annotation_col['SRR87031987',]<-'Hereford'
pheatmap(cor_mat,annotation_row = annotation_col,annotation_col = annotation_col,fontsize=18,
         show_rownames = F,show_colnames = F)

