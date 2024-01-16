#!/bin/Rscript
# BiocManager::install("snpStats")
# BiocManager::install("coloc")
library(coloc)
devtools::install_github("boxiangliu/locuscomparer")
library(locuscomparer)
library(RColorBrewer)

############################################################################################################
working_dir <- "path_to_your_working_dir"

setwd(working_dir)
setwd('G:/LHF/bovine circQTL/SNP data')

# eqtl <- read.table(file="./coloc/eqtl.txt", header=T, as.is=T)
# head(eqtl)
# gwas <- read.table(file="./coloc/gwas.tsv", header=T, as.is=T)
# head(gwas)

# # import gwas\exper qtl datasets
# gwas_raw<-read.table('./output/coloc/QTL_ARS_UCD1.gff.txt',comment.char = '#',sep='\t',fill = T)
# gwas_filt<-gwas_raw
# colnames(gwas_filt)<-c('chrom','Source','Association','start','end','un1','un2','un3','Details')
# # from chr.X to X
# # gwas_filt$chrom<-t(a)[,2]
# gwas_filt$chrom<-t(data.frame('chrom'=strsplit(gwas_filt$chrom,'.',fixed = T)))[,2]
# write.table(gwas_filt,'./output/temp.txt',sep='\t',quote = F,row.names = F)
# # Split details
# gwas_filt<-read.table('./output/temp.txt',sep='\t',header = T)
# gwas_filt<-gwas_filt[!duplicated(gwas_filt),]
# 
# temp<-strsplit(gwas_filt$Details[200],';',fixed = T)[[1]][1:14]
# tmp<-t(data.frame(strsplit(temp,'=',fixed = T)))
# filt_colnames<-tmp[,1][-9]
# gwas_filt<-gwas_filt[-48628,]
# for (i in 1:length(gwas_filt[,1])) {
#   ## only extract shared 14 cols
#   # [1] "QTL_ID=75091"                   "Name=Scrotal circumference"     "Abbrev=SCRCIR"
#   # [4] "PUBMED_ID=23785023"             "trait_ID=1270"                  "trait=Scrotal circumference"
#   # [7] "FlankMarker=rs134177902"        "breed=Tropical Composite"       "VTO_name=scrotum circumference"
#   # [10] "Map_Type=Genome"                "Model=Mendelian"                "Test_Base=Experiment-wise"
#   # [13] "Significance=Significant"       "P-value=<0.05"
#   line=gwas_filt$Details[i]
#   line=gsub('\xa0','',line)
#   detail<-strsplit(line,';',fixed = T)[[1]]
#   temp<-t(data.frame(strsplit(detail,'=',fixed = T)))
#   detail_colnames<-temp[,1]
#   detail<-data.frame(t(data.frame(temp[,2])))
#   colnames(detail)<-detail_colnames
#   check_res<-filt_colnames %in% colnames(detail)
#   index<-which(check_res==FALSE)
#   if (length(index)>0) {
#     index<-sort(index)
#     for (j in 1:length(index)) {
#       if (filt_colnames[index[j]]=='P-value'){
#         tmp<-data.frame(0.05)
#       } else {
#         tmp<-data.frame(NA) }
#       colnames(tmp)<-filt_colnames[index[j]]
#       detail<-cbind(detail,tmp) }
#     if (detail[,'P-value']=='<0.05') {
#       detail[,'P-value']=0.05 }
#     else if (detail[,'P-value']=='<0.01') {
#       detail[,'P-value']=0.01 }
#   }
#   detail_filt<-data.frame(detail[,filt_colnames])
#   if (i==1){
#     write.table(detail_filt,'./output/coloc/details.txt',quote = F,row.names = F,col.names = T,sep='\t')
#   } else {
#     write.table(detail_filt,'./output/coloc/details.txt',quote = F,row.names = F,col.names = F,append = T,sep='\t')
#   }
# } 

# details<-read.table('./output/coloc/details.txt',header = T,sep='\t')
# gwas_final<-cbind(gwas_filt[,1:5],details)
# write.table(gwas_final,'./coloc/gwas_final.txt',quote = F,row.names = F,sep='\t')


## import SNP datasets
snpspos<-read.table('./output/matrixeqtl/snp_pos.txt',header = F,sep = '\t')
colnames(snpspos)<-c('snps','chrom','loc')
head(snpspos)
gwas_final<-read.table('./output/coloc/gwas_final.txt',header = T,sep='\t')
gwas_final_filt <- gwas_final[(!is.na(gwas_final$start) & !is.na(gwas_final$end) 
                               & (gwas_final$Association=='Production_Association')),]
#19094 lead QTLs

# ARS_UCD1.2-BTAU5.0 location
total_filt<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)[,c(1,3,4,7,10,11)]
# total_sing_chipid2snpid
total_filt$chrom<-data.frame(t(data.frame(strsplit(total_filt$UMD3_1_1_Pos,':',fixed = T))))[,1]
total_filt$loc<-data.frame(t(data.frame(strsplit(total_filt$UMD3_1_1_Pos,':',fixed = T))))[,2]
GSE95358_map<-read.table('./output/remap/GSE95358.map',sep='\t',header = F)
colnames(GSE95358_map)<-c('chrom','Chip_ID','Group','loc')
chip_circqtl_sign=merge(GSE95358_map,total_filt,by=c('chrom','loc'))
# write.table(chip_circqtl_sign,'./coloc/chip_circqtl_sign.txt',sep='\t',quote = F,row.names = F)
#
colnames(total_filt)[1]<-'snps'
eqtl <- read.table(file="./output/circQTLs/hvaild_circQTLs.txt", header=T, as.is=T)
head(eqtl)
eqtl_filt<-merge(eqtl,total_filt,by='snps')
#UMD3_1_1_Pos
eqtl_filt$chrom<-data.frame(t(data.frame(strsplit(eqtl_filt$UMD3_1_1_Pos,':',fixed = T))))[,1]
eqtl_filt$loc<-data.frame(t(data.frame(strsplit(eqtl_filt$UMD3_1_1_Pos,':',fixed = T))))[,2]
#ars_ucd1.2
eqtl_filt$chrom<-data.frame(t(data.frame(strsplit(eqtl_filt$ARS_UCD1_2_Pos,':',fixed = T))))[,1]
eqtl_filt$loc<-data.frame(t(data.frame(strsplit(eqtl_filt$ARS_UCD1_2_Pos,':',fixed = T))))[,2]
GSE95358_map<-read.table('./output/remap/GSE95358.map',sep='\t',header = F)
colnames(GSE95358_map)<-c('chrom','Chip_ID','Group','loc')
chip_circqtl=merge(GSE95358_map,eqtl_filt,by=c('chrom','loc'))
# write.table(chip_circqtl,'./coloc/chip_circqtl.txt',sep='\t',quote = F,row.names = F)
# write.table(chip_circqtl$Chip_ID,'./coloc/chipid.txt',sep='\t',quote = F,row.names = F)
write.table(eqtl_filt,'./output/coloc/eqtl_filt.txt',quote = F,row.names = F,sep='\t')

## Figure 4A
##Venn plot
library(VennDiagram)
library(org.Bt.eg.db)
library(clusterProfiler)
library(ggplot2)
eqtl_filt<-read.table('./output/coloc/eqtl_filt.txt',header = T)
gwas_final<-read.table('./output/coloc/gwas_final.txt',header=T,sep='\t')
gwas_final_filt <- gwas_final[(!is.na(gwas_final$start) & !is.na(gwas_final$end) 
                               & (gwas_final$Association=='Production_Association')),]
circqtl_list<-eqtl_filt$snps;gwas_list<-gwas_final_filt$FlankMarker;
venn.diagram(x=list('circQTL'=circqtl_list,"GWAS-loci"=gwas_list),resolution =300, 
             filename = 'Venn plot of circQTLs and GWAS-loci.png',fill=c("cornflowerblue","green"),
             alpha = 0.50, cex=4, cat.cex=4,cat.pos=c(-20,20),cat.dist=0,main.cex = 4,
             ext.pos=0,ext.dist=-0.1)


total_filt<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
inter_set<-data.frame('VariantID'=intersect(circqtl_list,gwas_list))
## Enrichment analysis of gwas-circQTL compared to eQTL in different production traits
idx<-unlist(lapply(inter_set$VariantID, function(x){grep(x,gwas_final_filt$FlankMarker)}))
intersect_gwas_filt<-gwas_final_filt[idx,]
head(intersect_gwas_filt)
colnames(intersect_gwas_filt)
write.table(intersect_gwas_filt,'./output/coloc/intersect_gwas_filt.txt',quote = F,
            row.names = F,sep = "\t")
intersect_gwas_circqtl<-intersect_gwas_filt[,c("FlankMarker","trait")]
intersect_gwas_circqtl$isCircQTL<-T
intersect_gwas_circqtl<-intersect_gwas_circqtl[!duplicated(intersect_gwas_circqtl),]
# freqency of different traits for circQTLs
ggplot(intersect_gwas_circqtl)+geom_bar(aes(x='',fill=trait),stat = 'count',position='stack')+coord_polar(theta = 'y')+
  guides(fill=guide_legend(title = "Production traits",nrow = 8))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  
  theme(legend.text = element_text(face = 'italic',size = 12),
        legend.title = element_text(face = 'italic',size = 15),
        legend.position = 'bottom',
        axis.title = element_blank(),axis.ticks = element_blank(),
        
        axis.text.y = element_text(size=15))+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(10)))

# read eQTLs
eQTLs<-read.table('./output/matrixeqtl/eqtl_filt_res.txt',header = T,sep='\t')
intersect_eQTLs<-intersect(eQTLs$snps,gwas_final_filt$FlankMarker)
eQTL_idx<-unlist(lapply(intersect_eQTLs, function(x){grep(x,gwas_final_filt$FlankMarker)}))
intersect_gwas_eqtl<-gwas_final_filt[eQTL_idx,]
intersect_gwas_eqtl<-intersect_gwas_eqtl[,c("FlankMarker","trait")]
intersect_gwas_eqtl$isCircQTL<-F
# combine
intersect_gwas_qtls<-rbind(intersect_gwas_circqtl,intersect_gwas_eqtl)
# get all production traits
type<-unique(intersect_gwas_qtls$trait)
pvalue<-list();OR<-list()
for (i in 1:length(type)) {
  # i=1
  # i=1
  #beta>0 & isin the type
  a<-nrow(intersect_gwas_qtls[intersect_gwas_qtls$isCircQTL==T&intersect_gwas_qtls$trait==type[i],])
  b<-nrow(intersect_gwas_qtls[intersect_gwas_qtls$isCircQTL==T&intersect_gwas_qtls$trait!=type[i],])
  c<-nrow(intersect_gwas_qtls[intersect_gwas_qtls$isCircQTL==F&intersect_gwas_qtls$trait==type[i],])
  d<-nrow(intersect_gwas_qtls[intersect_gwas_qtls$isCircQTL==F&intersect_gwas_qtls$trait!=type[i],])
  table<-data.frame('Is'=c(a,c),'Not'=c(b,d))
  print(table)
  res<-fisher.test(table,alternative = 'greater')
  pvalue[i]<-res$p.value
  OR[i]<-res$estimate
  
}
## plot
enrich_res<-data.frame('Type'=type,'pvalue'=unlist(pvalue),'OR'=unlist(OR))
enrich_res<-enrich_res[enrich_res$OR!=Inf,]
highlight<-enrich_res[(enrich_res$pvalue<=10**(-0.25))&(enrich_res$OR>=1),]
rand<-runif(length(highlight[,1]),0,0.2)
rand1<-runif(length(highlight[,1]),-0.1,0.1)
highlight$ORend<-highlight$OR-rand;highlight$pvalueend<-highlight$pvalue+rand1
n=length(enrich_res[,1])-length(highlight[,1])
blank<-data.frame('Type'=rep('',n),'pvalue'=rep(0,n),'OR'=rep(0,n),'pvalueend'=0,'ORend'=0)
highlight<-rbind(highlight,blank)
highlight$ORend[highlight$Type=='non_coding_transcript_exon_variant']=1.28
# bubble plot
ggplot(enrich_res)+geom_point(aes(x=OR,y=-log10(pvalue)),color='lightblue',size=10,alpha=0.8)+
  xlab('Odds ratio #(circQTL)/#(eQTL)')+
  ylab(expression(paste(-log[10],'(',italic("P")," value)",sep='')))+
  geom_vline(xintercept = 1,color='red',linetype='dashed')+
  geom_hline(yintercept = -log10(0.05),color='red',linetype='dashed')+
  # geom_text(aes(x=highlight$ORend,y=-log10(highlight$pvalueend),label=highlight$Type),size=5)+
  # annotate('segment', x=highlight$OR, xend=highlight$OR+(highlight$ORend-highlight$OR)*0.8, 
  # y=-log10(highlight$pvalue), 
  # yend=-log10(highlight$pvalue)+(-log10(highlight$pvalueend)+log10(highlight$pvalue))*0.8,
  # arrow=arrow(angle = 30,length = unit(0,'cm')))+
  ggrepel::geom_text_repel(data = enrich_res[-log10(enrich_res$pvalue)>0.25&enrich_res$OR>1,],
                           aes(x=OR,y=-log10(pvalue),label=Type),size=8)+
  theme_bw()+theme(axis.title.x = element_text(size=25),
                   axis.title.y = element_text(size=25),
                   axis.text = element_text(size=18),
                   legend.title = element_text(size=15),
                   legend.text = element_text(size=12))

## Figure 4A-go enrichment
## Intersect go analysis
inter_temp<-mergeCustom(total_filt,inter_set,xcol='VariantID',ycol="VariantID")
# read breed freq matrix
breed_freq<-read.table('./output/matrixeqtl/snp.txt',header = T)
colnames(breed_freq)[1]<-"VariantID"
inter_temp_frep<-mergeCustom(inter_temp,breed_freq,xcol="VariantID",ycol='VariantID')
write.table(inter_temp_frep,'./output/coloc/inter_gwas_loci_circQTL.txt',sep='\t',row.names = F,quote = F)
## import correlation of qtl with circRNA
hfilt_circQTL<-read.table('./output/circQTLs/hfilt_circQTL.txt',sep='',header = T)
temp<-data.frame("snps"=inter_temp[,"VariantID"])
hfilt_circQTL_gwas<-merge(hfilt_circQTL,temp,by="snps")
write.table(hfilt_circQTL_gwas,'./output/coloc/hfilt_circQTL_gwas.csv',sep=',',quote = F,row.names = F)

inter_symbol<-unlist(lapply(inter_temp$Gene,function(x){strsplit(x,',',fixed = T)[[1]][1]}))
inter_snp<-data.frame('VariantID'=inter_temp$VariantID,'SYMBOL'=inter_symbol)

eg_total = bitr(inter_symbol, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Bt.eg.db")
eg_total = eg_total[!duplicated(eg_total$SYMBOL),]
eg_total = merge(eg_total,inter_snp,by='SYMBOL')
write.table(eg_total,'./output/coloc/intersect.txt',quote = F,row.names = F,sep='\t')

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
# GO_result$GO_term_1<-rep(c('BP','MF','CC'),rep(10,3))
write.table(GO_result,'./output/coloc/Go enrichment analysis for GWAS-related circQTLs.txt',
            quote = F,row.names = F,sep = "\t")
GO_result$'-log10(qvalue)'= (-1)*log10(GO_result$qvalue)
# pdf(paste(output, '/',name," Go enrichment analysis.pdf", sep=""),width=12,height=12)
# p<-ggplot(GO_result,mapping = aes(x=reorder(Description,log10(qvalue)),y=-log10(qvalue),fill=GO_term))+
#   geom_bar(stat = "identity")+scale_x_discrete(limits=GO_result$Description)+theme_bw()+xlab('')+
#   theme(axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
#         axis.text.x = element_text(size=18,angle = 80,hjust=1),axis.text.y = element_text(size=15))
# print(p)
# Point plot
ggplot(GO_result,mapping = aes(y=reorder(Description,log10(qvalue)),x=GO_term,size=-log10(qvalue),color=GO_term))+
  geom_point()+theme_bw()+xlab('')+
  theme(axis.title.x = element_text(size=25),axis.title.y = element_blank(),
        axis.text.x = element_text(size=18,angle = 80,hjust=1),axis.text.y = element_text(size=15))

kk <- enrichKEGG(gene = id,
                 organism     = 'bta',
                 pAdjustMethod = "BH",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05)
KEGG_result<-subset(kk@result,select=c('Description','qvalue'))[1:10,]
KEGG_result$'-log10(qvalue)'= (-1)*log10(KEGG_result$qvalue)
# pdf(paste(output, '/',name," KEGG enrichment analysis.pdf", sep=""),width=12,height=12)
p<-ggplot(KEGG_result,mapping = aes(x=reorder(Description,log10(qvalue)),y=-log10(qvalue)))+
  geom_bar(stat = "identity",fill='blue')+scale_x_discrete(limits=KEGG_result$Description)+theme_bw()+xlab('')+
  theme(axis.title.x = element_text(size=18),axis.title.y = element_blank(),
        axis.text.x = element_text(size=18,angle = 80,hjust=1),axis.text.y = element_text(size=15))
print(p)

# diff_set<-setdiff(union(circqtl_list,gwas_list),inter_set)

## Figure 4B
## Plot LDblockheatmap
library(gpart)
library(genetics)
# read the sign_circqtls
# total_filt<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
## circQTLs with duplicated because of different circQTL-circRNA pairs
total_filt<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt_duplicated.txt',header = T)
# 
sign_circqtl<-total_filt[,c(1,2,11)]
sign_circqtl$chrom<-as.numeric(data.frame(t(data.frame(strsplit(sign_circqtl$ARS_UCD1_2_Pos,':',fixed = T))))[,1])
sign_circqtl$loc<-as.numeric(data.frame(t(data.frame(strsplit(sign_circqtl$ARS_UCD1_2_Pos,':',fixed = T))))[,2])
unique_circqtls<-unique(sign_circqtl$VariantID)
sign_circqtl_nodup<-data.frame('VariantID'=unique_circqtls)
for (i in 20108:length(unique_circqtls)){
  print(i)
  sign_circqtl_nodup[i,'pvalue']=min(sign_circqtl$pvalue[sign_circqtl$VariantID==unique_circqtls[i]])
  sign_circqtl_nodup[i,'ARS_UCD1_2_Pos']=sign_circqtl$ARS_UCD1_2_Pos[sign_circqtl$VariantID==unique_circqtls[i]][1]
}
# read the sign_circqtls
total_filt_nodup<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
total_filt_nodup<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)

# sign_circqtl_nodup<-total_filt_nodup[,c(1,2,11)]
sign_circqtl_nodup$chrom<-as.numeric(data.frame(t(data.frame(strsplit(sign_circqtl_nodup$ARS_UCD1_2_Pos,':',fixed = T))))[,1])
sign_circqtl_nodup$loc<-as.numeric(data.frame(t(data.frame(strsplit(sign_circqtl_nodup$ARS_UCD1_2_Pos,':',fixed = T))))[,2])

# read the GSE100038 genotype datasets
chip_total<-read.table('./output/imputation/GSE100038/GSE100038_imp.vcf', sep='')
colnames(chip_total)[3]<-'Chip_ID'
chip_total<-chip_total[,-c(1:2,4:9)]
# read the chipid-snpid info
chip2snp<-read.table('./coloc/chip_circqtl_sign.txt',header = T)[,c(3,5)]
chip_total<-merge(chip_total,chip2snp,by='Chip_ID')
chip_total<-subset(chip_total,select = -Chip_ID)
# combine the circqtls and gwas-loci
inputs_gwas<-data.frame('rsid'=gwas_final_filt$FlankMarker,'pval'=gwas_final_filt$P.value)
# inputs_eqtl<-data.frame('rsid'=eqtl_filt$snps,'pval'=eqtl_filt$pval)
inputs_eqtl<-data.frame('rsid'=eqtl_filt$snps,'pval'=eqtl_filt$pvalue,'MAF'=eqtl_filt$MAF)
inputs<-merge(inputs_eqtl, inputs_gwas, by="rsid", suffixes=c("_eqtl","_gwas"))
inputs<-inputs[!duplicated(inputs$rsid),]
colnames(inputs)[1]<-'snps'
# add the chr and pos to inputs
inputs_prep<-merge(inputs,eqtl_filt[,c(1,12,13)],by='snps')
# inputs_prep<-merge(inputs,eqtl_filt,by='snps')

inputs_prep<-inputs_prep[!duplicated(inputs_prep$snps),]
#chr levels
chr_levels<-sort(inputs_prep$chrom[!duplicated(inputs_prep$chrom)]) # character

# gene_info
bos_ref <- read.table('./bos_ref.txt',header = T,sep='')
geneinfo <- data.frame(bos_ref$name2, bos_ref$chrom,bos_ref$txStart,bos_ref$txEnd)
colnames(geneinfo)<-c('genename','chrN','start.bp','end.bp')
chrom<-lapply(list(geneinfo$chrN), function(x){ res <- substr(x,4,length(x))})
geneinfo$chrN <- unlist(chrom)
geneinfo<-geneinfo[!duplicated(geneinfo),]
# unique_genename<-unique(geneinfo$genename)
# geneinfo_nodup<-data.frame()
# for (i in 1:length(unique_genename)) {
#   geneinfo_nodup[i,'genename']<-unique_genename[i]
#   geneinfo_nodup[i,'chrN']<-unique(geneinfo$chrN[geneinfo$genename==unique_genename[i]])
#   geneinfo_nodup[i,'start.bp']<-min(geneinfo$start.bp[geneinfo$genename==unique_genename[i]])
#   geneinfo_nodup[i,'end.bp']<-max(geneinfo$end.bp[geneinfo$genename==unique_genename[i]])
# }

# cM/Mb=1.25
# 5% crossover=5cM=4Mb
flanking = 2e6

# convert A/A to 0,1,2 
conv_AB2number <- function(x){
  if ((is.na(x)) | (x == './.')) {
    res <- NA } else if (x=='A/A') {
      res <- 0 } else if (x=='A/B') {
        res<-1 } else {
          res<-2
        }
  return (res)
}
# convert 0|0 to A/A
conv_AB <- function(x){
  if ((is.na(x)) | (x == '.|.')) {
    res <- NA } else if (x=='0|0') {
      res <- 'A/A' } else if (x=='0|1') {
        res<-'A/B' } else {
          res<-'B/B'
        }
  return (res)
}

#variant_ID(chip_ID) list type
chip_ID_total<-list()
for (i in 1:length(chr_levels)){
  # i=3
  # get the length of each chrom
  chr_length<-max(geneinfo$end.bp[geneinfo$chrN==chr_levels[i]])
  
  inputs_each<-inputs_prep[inputs_prep$chrom==chr_levels[i],]
  
  loc<-sort(as.numeric(inputs_each$loc))
  for (j in 1:length(loc)) {
    # j=1
    # if ((i ==1 | 2) & (j == 2 | 3)) {
    #   next
    # }
    # print(j)
    # j=1
    start<-max(loc[j]-flanking,1);end<-min(loc[j]+flanking,chr_length)
    circqtls<-sign_circqtl_nodup[sign_circqtl_nodup$chrom==chr_levels[i]&sign_circqtl_nodup$loc>=start&sign_circqtl_nodup$loc<=end,]
    circqtls<-na.omit(circqtls)
    circqtls_each_loc<-circqtls$VariantID
    variant_each_loc<-data.frame('VariantID'=circqtls_each_loc)
    variant<-merge(chip_total,variant_each_loc,by='VariantID')
    id<-data.frame('rsID'=variant$VariantID)
    chip_each<-subset(variant, select = -c(VariantID))
    # print(VariantID)
    print('ok')
    # SNPinfo <- data.frame(circqtls$chrom,variant$VariantID,circqtls$loc)
    # colnames(SNPinfo)<-c('chrN','rsID','bp')
    # chip_each<-data.frame(t(chip_each))
    # colnames(chip_each)<-SNPinfo$rsID
    # rownames(chip_each)<-1:nrow(chip_each)
    # chip_each
    # calculate LD
    for(h in 1:ncol(chip_each)){
      chip_each[,h]<-unlist(lapply(chip_each[,h],conv_AB))
      chip_each[,h]<-unlist(lapply(chip_each[,h],conv_AB2number))
      # chip_each[,h]<-as.genotype(chip_each[,h])
    }
    # 
    # rownames(chip_each)<-unlist(c(1:nrow(chip_each)))
    # chip_each<-gsub(chip_each,pattern = c("A:A/","B:G/"))
    # 
    # 
    # ld_each<-LDheatmap(chip_each,LDmeasure="r")$LDmatrix
    
    # export the variant of each loc in eahc chrom
    # chip_ID_total[[i]][[j]]<-circqtl
    
    # LDblock
    SNPinfo <- data.frame(circqtls$chrom,circqtls$VariantID,circqtls$loc)
    colnames(SNPinfo)<-c('chrN','rsID','bp')
    SNPinfo<-merge(SNPinfo,id,by='rsID')[,c(2,1,3)]
    # SNPinfo<-SNPinfo[!duplicated(SNPinfo$rsID),]
    chip_each<-data.frame(t(chip_each))
    colnames(chip_each)<-SNPinfo$rsID
    LDblockHeatmap(geno = chip_each,SNPinfo = SNPinfo,geneinfo = geneinfo,
                   LD = "r2", filename = paste("./coloc/LDblock-3/heatmap_",'chr',chr_levels[i],"_",loc[j],sep=''),
                   blocktype = c("bigld"),CLQshow=T,ensbversion=101,
                   geneshow=F, MAFcut = 0.05,res = 300)
  }
  # # coloc
  # result_each <- coloc.abf(
  #   dataset1=list(pvalues=as.numeric(inputs_each$pval_gwas), type="cc", s=0.33, N=nrow(inputs_each),
  #                 snp=inputs_each$snps,MAF=inputs_each$MAF),
  #   dataset2=list(pvalues=inputs_each$pval_eqtl, type="quant", N=nrow(inputs_each),
  #                 snp=inputs_each$snps,MAF=inputs_each$MAF))
}


### locuszoom
## circQTL pvalue with single pvalue or mulit-pvalue
head(eqtl_filt)
circqtl_locuszoom<-sign_circqtl_nodup[,c("VariantID","pvalue")]

colnames(circqtl_locuszoom)<-c("MarkerName","P-value")
write.table(circqtl_locuszoom,'./output/coloc/circqtl_locuszoom.txt',sep='\t',row.names = F,quote = F)

snp_pos_db_file<-sign_circqtl_nodup[,c(1,4,5)]
colnames(snp_pos_db_file)<-c("snp",'chr','pos')
write.table(snp_pos_db_file,'./output/coloc/snp_pos_db_file.txt',sep='\t',row.names = F,quote = F)
## refFlat file
ref_flat_db_file<-bos_ref[,c("name2","name","chrom",        "strand",       "txStart",      "txEnd",
                             "cdsStart",     "cdsEnd","exonCount",    "exonStarts",   "exonEnds" )]
colnames(ref_flat_db_file)[1]<-c("geneName")
write.table(ref_flat_db_file,'./output/coloc/ref_flat_db_file.txt',sep='\t',row.names = F,quote = F)
## Figure 4B-bottom
library(LDheatmap)
coloc_sites<-data.frame('coloc_sites'=inter_set,'Marked'='*')
# heatmap_chr23-2001057 heatmap_chr7-48680543 heatmap_chr7-50794543
# chr<-29;loc<-46456814;chr_length<-max(geneinfo$end.bp[geneinfo$chrN==chr])
# chr<-23;loc<-2001057;chr_length<-max(geneinfo$end.bp[geneinfo$chrN==chr])
# chr<-7;loc<-48680543;chr_length<-max(geneinfo$end.bp[geneinfo$chrN==chr])
chr<-7;loc<-50794543;chr_length<-max(geneinfo$end.bp[geneinfo$chrN==chr])

start<-max(loc-flanking,1);end<-min(loc+flanking,chr_length)
circqtls<-sign_circqtl_nodup[sign_circqtl_nodup$chrom==chr&
                               sign_circqtl_nodup$loc>=start&sign_circqtl_nodup$loc<=end,]
circqtls<-na.omit(circqtls)
circqtls_each_loc<-circqtls$VariantID
variant_each_loc<-data.frame('VariantID'=circqtls_each_loc)
variant<-merge(chip_total,variant_each_loc,by='VariantID')
id<-data.frame('rsID'=variant$VariantID)
chip_each<-subset(variant, select = -c(VariantID))

SNPinfo <- data.frame(circqtls$chrom,variant$VariantID,circqtls$loc)
colnames(SNPinfo)<-c('chrN','rsID','bp')
chip_each<-data.frame(t(chip_each))
colnames(chip_each)<-SNPinfo$rsID
rownames(chip_each)<-1:nrow(chip_each)

# calculate LD
for(h in 1:ncol(chip_each)){
  chip_each[,h]<-unlist(lapply(chip_each[,h],conv_AB))
  chip_each[,h]<-as.genotype(chip_each[,h])
}
ld_each<-LDheatmap(chip_each,LDmeasure="r")$LDmatrix
diag(ld_each)<-1

key_loci<-intersect(colnames(ld_each),inter_set$VariantID)
# key_loci="rs29015854"
loci<-key_loci[[1]];index<-which(colnames(ld_each)==loci)
print(loci)
ld_r2=data.frame('R2'=c(ld_each[1:index,index],ld_each[index,(index+1):nrow(ld_each)]))
ld_r2$VariantID = rownames(ld_r2)
locus_data<-subset(sign_circqtl_nodup,select=c(VariantID,pvalue,chrom,loc))
# locus_data<-locus_data[!duplicated(locus_data$VariantID),]
locus_data<-merge(locus_data,ld_r2,by='VariantID',all=FALSE)
locus_data$Type<-'circQTLs'
locus_data$Type[locus_data$VariantID==loci]='Gwas/Experiment based-loci'
locus_data_filtout<-locus_data
ld_data_within_locus<-locus_data_filtout[,c("VariantID","R2")]
ld_data_within_locus$snp2<-loci;ld_data_within_locus$dprime<-NA
ld_data_within_locus<-ld_data_within_locus[,c("VariantID","snp2","dprime","R2")]
colnames(ld_data_within_locus)<-c("snp1","snp2","dprime","rsquare")
write.table(ld_data_within_locus,'./coloc/ld_data_within_locus.txt',quote = F,row.names = F,sep=' ')
# locus_data_filtout$VariantID<-runif()
palette<-brewer.pal(7, "Spectral")
ggplot(locus_data)+geom_point(aes(x=loc/1e6,y=-log10(pvalue),color=R2,shape=Type),size=5,alpha=0.8)+
  scale_color_gradientn(colours = rev(palette))+
  scale_x_continuous(breaks=seq(round(min(locus_data$loc)/1e6,3),round(max(locus_data$loc)/1e6,3),1e6/1e6),
                     limits=c(min(locus_data$loc)/1e6,max(locus_data$loc)/1e6))+theme_bw()+
  xlab(paste('Genomic location (Mb) in ','chr',chr,sep=''))+
  guides(shape=guide_legend(order = 1,title="Type"),color = guide_colourbar(order = 2,title=expression("R"^{2})))+
  theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),
        axis.text.x = element_text(size=18,angle = 80,hjust=1),axis.text.y = element_text(size=18),
        legend.text = element_text(size=15),legend.title = element_text('R^2',size=18))+
  annotate('text',x=locus_data$loc[locus_data$VariantID==loci][1]/1e6,
           y=-log10(min(locus_data$pvalue[locus_data$VariantID==loci]))+0.2,label=loci,size=7)
