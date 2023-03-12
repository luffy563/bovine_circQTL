#!/bin/Rscript

library(circRNAprofiler)
library(sankeyNetwork)
library(tidyverse)
library(tidyr)
library(tibble)
library(viridis)
library(patchwork)
library(networkD3)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggstatsplot)
library(sankeyD3)
library(magrittr)
library(pheatmap)
library(ggpubr)
library(VennDiagram)
library(gridExtra)
library(biomaRt)
library(stringr)
library(parallel)

############################################################################################################
working_dir <- "path_to_your_working_dir"
circRNAprofiler_dir<-paste(working_dir,'output/RBPpred/predCirc',sep='/')
all_RNA_seq_file='./SRR_list_CIRI.txt'
output_dir <- paste(working_dir,'./output/RBPpred',sep='/')
setwd(circRNAprofiler_dir)
dir.create(output_dir)
dir.create(circRNAprofiler_dir)
cores<-3
cl<-makeCluster(cores)

############################################################################################################
## Prediction of RBP sites and miRNA sponge sites

mart<-useMart("ensembl")
mart <- useMart("ensembl","btaurus_gene_ensembl")

## import circRNA coordinate
circRNA<-read.table('./output/anno_info_geneBankId.csv',sep=',',header = T)
head(circRNA)
colnames(circRNA)[1]<-'circRNA_ID'

temp<-circRNA$start[circRNA$strand=='-']
circRNA$start[circRNA$strand=='-']<-circRNA$end[circRNA$strand=='-']
circRNA$end[circRNA$strand=='-']<-temp
circRNA<-circRNA[-6]
circRNA$genenBank_ID<-unlist(lapply(circRNA$genenBank_ID,function(x){strsplit(x,'.',fixed = T)[[1]][1]}))
id_list<-c('refseq_mrna','refseq_mrna_predicted','refseq_ncrna')
names(id_list)<-c('NM','XM','XR')

for (i in 1:3){
	idx<-!is.na(str_match(circRNA$genenBank_ID,names(id_list[i])))
	geneid<-circRNA$genenBank_ID[idx]

	hg_symbols<- getBM(attributes=c(id_list[i],'ensembl_gene_id','external_gene_name'), 
						   filters= id_list[i], values = geneid, mart = mart)
	colnames(hg_symbols)[1]<-'genenBank_ID'
	hg_symbols<-hg_symbols[!duplicated(hg_symbols$genenBank_ID),]
	a<-lapply(circRNA$genenBank_ID[idx],function(x){which(hg_symbols$genenBank_ID==x)})
	index<-as.numeric(unlist(lapply(a,function(x){if(length(x)==0){return('NA')}else{return(x)}})))
	circRNA$ensembl_geneid[idx]<-hg_symbols$ensembl_gene_id[index]
	# merge(circRNA,hg_symbols,by='genenBank_ID',all.x = T)
	circRNA$gene[idx]<-hg_symbols$external_gene_name[index]
}
circRNA<-subset(circRNA,select = -c(genenBank_ID,ensembl_geneid))
circRNA$coverage<-1
circRNA<-circRNA[-1]
colnames(circRNA)<-c('chrom','startUpBSE', 'endDownBSE','strand', 'gene','coverage')
write.table(circRNA,'./output/RBPpred/predCirc/other/circRNA_001.txt',sep='\t',quote = F,row.names = F)

############################################################################################################
### Initialization
setwd('./output/RBPpred/predCirc')
check <- checkProjectFolder()
check
# formats the given annotation file 
gtf <- formatGTF("./refrence/genome/Bos_taurus.ARS-UCD1.2.101.gtf")
# reads the circRNAs_X.txt 
backSplicedJunctions<-getBackSplicedJunctions()

# mergedBSJunctions <- mergeBSJunctions(backSplicedJunctions, gtf)
mergedBSJunctions<-backSplicedJunctions
# filteredCirc <-
#   filterCirc(mergedBSJunctions, allSamples = FALSE, min = 5)

filteredCirc<-mergedBSJunctions

## Annotate circRNAs internal structure and flanking introns
# annotatedBSJs<-read.table('./output/RBPpred/predCirc/annotatedBSJs.txt',header = T,sep='\t')
annotatedBSJs <- annotateBSJs(filteredCirc, gtf, isRandom = FALSE)
head(annotatedBSJs)
write.table(annotatedBSJs,'./output/RBPpred/predCirc/annotatedBSJs.txt',quote = F,row.names = F,sep='\t')

f <- sum((annotatedBSJs$exNumUpBSE == 1 |
         annotatedBSJs$exNumDownBSE == 1) ,
      na.rm = TRUE) / (nrow(annotatedBSJs) * 2)

# Retrieve random back-spliced junctions
randomBSJunctions <-
  getRandomBSJunctions(gtf, n = nrow(annotatedBSJs), f = f*100, setSeed = 123)
head(randomBSJunctions)

# annotatedRBSJs<-read.table('./output/RBPpred/predCirc/annotatedRBSJs.txt',header = T,sep='\t')
annotatedRBSJs <- annotateBSJs(randomBSJunctions, gtf, isRandom = TRUE)
write.table(annotatedRBSJs,'./output/RBPpred/predCirc/annotatedRBSJs.txt',quote = F,row.names = F,sep='\t')

# annotatedRBSJs<-read.table('./output/RBPpred/predCirc/annotatedRBSJs.txt',header = T,sep='\t')

# Length of flanking introns
p1 <- plotLenIntrons(
  annotatedBSJs,
  annotatedRBSJs,
  title = "Length flanking introns",
  df1Name = "predicted",
  df2Name = "random",
  setyLim = TRUE, 
  ylim = c(0,7)
)

# Length of back-splided exons
p2 <- plotLenBSEs(
  annotatedBSJs,
  annotatedRBSJs,
  title = "Length back-splided exons",
  df1Name = "predicted",
  df2Name = "random",
  setyLim = TRUE, 
  ylim = c(0,7)
)

# No. of circRNAs produced from the host genes
p3 <-
  plotHostGenes(annotatedBSJs, title = "# CircRNAs produced from host genes")

# No. of exons in between the back-spliced junctions
p4 <-
  plotExBetweenBSEs(annotatedBSJs, title = "# Exons between back-spliced junctions")

# Position of back-spliced exons within the host transcripts
p5 <-
  plotExPosition(annotatedBSJs,
                 n = 1,
                 title = "Position back-spliced exons in the transcripts")

# Total no. of exons within the host transcripts
p6 <-
  plotTotExons(annotatedBSJs, title = " Total number of exons in the host transcripts")

#####################################
## Figure
#####################################
# Combine plots
ggarrange(p1,
          p2,
          p3,
          p4,
          p5,
          p6,
          ncol = 2,
          nrow = 3)

# Get genome
genome <- BSgenome::getBSgenome("BSgenome.Btaurus.UCSC.bosTau9")

# Foreground target sequences
# load('./output/RBPpred/predCirc/targetsFTS_circ.Rdata')
targetsFTS_circ <-
  getCircSeqs(annotatedBSJs, gtf, genome)
meanlength<-mean(na.omit(annotatedBSJs$lenCircRNA))

idx<-is.na(str_match(annotatedRBSJs$chrom,'chrNKLS'))
annotatedRBSJs<-annotatedRBSJs[idx,]

# load('./output/RBPpred/predCirc/targetsRFTS_circ.Rdata')
targetsRFTS_circ <-
  getCircSeqs(annotatedRBSJs, gtf, genome)
meanRlength<-mean(na.omit(annotatedRBSJs$lenCircRNA))
## random seq
library(Biostrings)
library(BSgenome)
# random_seq<-read.table('./reference/genome/random_seq.fasta')
bta_random_file <- system.file("extdata", "./reference/genome/random_seq.fasta", package="BSgenome")
bta_random_seq<-readDNAStringSet("./reference/genome/random_seq.fasta", "fasta")
bta_random_seq<-data.frame(bta_random_seq)
targetsRFTS_circ$circ$seq<-bta_random_seq$bta_random_seq[1:length(targetsRFTS_circ$circ$seq)]
targetsRFTS_circ$circ$seq<-str_replace_all(targetsRFTS_circ$circ$seq,'T','U')
save(targetsFTS_circ,file='./targetsFTS_circ.Rdata')
save(targetsRFTS_circ,file='./targetsRFTS_circ-1.Rdata')
save(targetsRFTS_circ,file='./targetsRFTS_circ-random_genome.Rdata')
load('./targetsRFTS_circ-random_genome.Rdata')
write.table(targetsFTS_circ$circ,'./output/RBPpred/predCirc/targetsFTS_circ_seq.txt',quote = F,row.names = F,sep='\t')
write.table(targetsRFTS_circ$circ,'./output/RBPpred/predCirc/targetsRFTS_circ_seq-1.txt',quote = F,row.names = F,sep='\t')
write.table(targetsRFTS_circ$circ,'./output/RBPpred/predCirc/targetsRFTS_circ_seq-random_genome.Rdata.txt',quote = F,row.names = F,sep='\t')


# Atract the BSJ seq of Foreground target sequences
targetsFTS_bsj <-
  getSeqsAcrossBSJs(annotatedBSJs, gtf, genome)

## extract the BSJ region seqs
# Foreground target sequences
targetsFTS_gr <-
  getSeqsFromGRs(
    annotatedBSJs,
    genome,
    lIntron = 500,
    lExon = 9,
    type = "ie"
  )
# Background target sequences.
targetsBTS_gr <-
  getSeqsFromGRs(
    annotatedRBSJs,
    genome,
    lIntron = 500,
    lExon = 9,
    type = "ie")

## save the BSJ region seqs of 500bp intron and 9bp exon
save(targetsFTS_gr,file='./output/RBPpred/predCirc/targetsFTS_gr_500bp.Rdata')
save(targetsBTS_gr,file='./output/RBPpred/predCirc/targetsBTS_gr_500bp.Rdata')

load('./targetsFTS_gr.Rdata')
load('./targetsBTS_gr.Rdata')
save(targetsFTS_gr,file='./output/RBPpred/predCirc/targetsFTS_gr.Rdata')
save(targetsBTS_gr,file='./output/RBPpred/predCirc/targetsBTS_gr.Rdata')

# Find motifs in the predicted target sequences
motifsFTS_gr <-
  getMotifs(targetsFTS_gr,
            width = 6,
            database = 'ATtRACT',
            species = "Btaurus",
            rbp = TRUE,
            reverse = FALSE)
# Find motifs in the background target sequences
motifsBTS_gr <-
  getMotifs(targetsBTS_gr,
            width = 6,
            database = 'ATtRACT',
            species = "Btarurs",
            rbp = TRUE,
            reverse = FALSE)

load('./output/RBPpred/predCirc/motifsBTS_gr.Rdata')
save(motifsFTS_gr,file='./output/RBPpred/predCirc/motifsFTS_gr.Rdata')
save(motifsBTS_gr,file='./output/RBPpred/predCirc/motifsBTS_gr.Rdata')
## 
BiocManager::install('ggseqlogo')
BiocManager::install('seqLogo')
library(ggseqlogo)
library(seqLogo)
# Find motifs in the predicted internal sequences
motifsFTS_circ <-
  getMotifs(targetsFTS_circ,
            width = 6,
            database = 'ATtRACT',
            species = "Btarurs",
            rbp = TRUE,
            reverse = FALSE) 
motifsRFTS_circ <-
  getMotifs(targetsRFTS_circ,
            width = 6,
            database = 'ATtRACT',
            species = "Btarurs",
            rbp = TRUE,
            reverse = FALSE)

save(motifsFTS_circ,file='./output/RBPpred/predCirc/motifsFTS_circ.Rdata')
save(motifsRFTS_circ,file='./output/RBPpred/predCirc/motifsRFTS_circ_random.Rdata')
save(motifsRFTS_circ,file='./output/RBPpred/predCirc/motifsRFTS_circ.Rdata')

#
load('./output/RBPpred/predCirc/motifsFTS_circ.Rdata')
# Merge all the motifs by the same RBP 
mergedMotifsFTS_circ <- mergeMotifs(motifsFTS_circ)
mergedMotifsRFTS_circ <- mergeMotifs(motifsRFTS_circ)
write.table(mergedMotifsFTS_circ,'./output/RBPpred/predCirc/mergedMotifsFTS_circ.txt',quote = F,row.names = F,sep='\t')
write.table(mergedMotifsRFTS_circ,'./output/RBPpred/predCirc/mergedMotifsRFTS_circ.txt',quote = F,row.names = F,sep='\t')
write.table(mergedMotifsRFTS_circ,'./output/RBPpred/predCirc/mergedMotifsRFTS_circ_random.txt',quote = F,row.names = F,sep='\t')


mergedMotifsFTS_gr <- mergeMotifs(motifsFTS_gr)
mergedMotifsBTS_gr <- mergeMotifs(motifsBTS_gr)
write.table(mergedMotifsFTS_gr,'./output/RBPpred/predCirc/mergedMotifsFTS_gr_500bp.txt',quote = F,row.names = F,sep='\t')
write.table(mergedMotifsFTS_gr,'./output/RBPpred/predCirc/mergedMotifsBTS_gr_500bp.txt',quote = F,row.names = F,sep='\t')

mergedMotifsFTS_circ<-read.table('./output/RBPpred/predCirc/mergedMotifsFTS_circ.txt',header = T,sep='\t')
mergedMotifsRFTS_circ <- read.table('./output/RBPpred/predCirc/mergedMotifsRFTS_circ.txt',header = T,sep='\t')

## barplot of BSJ region motif
mergedMotifsFTS_gr<-mergedMotifsFTS_gr[order(mergedMotifsFTS_gr$count,decreasing = T),]
# ggplot(mergedMotifsFTS_gr[1:10,])+geom_bar(aes(x=count,y=reorder(id,count)),stat = 'identity')+
#   xlab('Count')+ylab('RBP')+
#   theme_classic()+theme(axis.title.x = element_text(size=18),
#                         axis.title.y = element_text(size=18),
#                         axis.text.x = element_text(angle=45,hjust=1,size=18),
#                         axis.text.y = element_text(size=15),
#                         legend.title = element_text(size=15),
#                         legend.text = element_text(size=12))

p <-
  plotMotifs(
    mergedMotifsFTS_gr,
    mergedMotifsBTS_gr,
    nf1 = (210*2)* nrow(annotatedBSJs) , 
    nf2 = (210*2)* nrow(annotatedRBSJs), 
    log2FC = 1,
    removeNegLog2FC = TRUE,
    df1Name = "all cricRNAs",
    df2Name = "random circRNAs",
    angle = 45
  )
ggarrange(p[[1]],
          p[[2]],
          labels = c("", ""),
          ncol = 2,
          nrow = 1)

## RBP of internal seq within CircRNA
# read related seqlogos
tf_names <- c('CTCF', 'SP1')
rbp_name<-rev(p[[1]]$data$id)[3]
rbp_db<-read.table('reference/ATtRACT/ATtRACT_db.txt',sep='\t',header = T)
head(rbp_db)
rbp_db_filt<-rbp_db[rbp_db$Organism=="Homo_sapiens",c('Gene_name',"Motif","Len")]
rbp_db_filt<-rbp_db_filt[!duplicated(rbp_db_filt),]
motifs<-read.table('./output/RBPpred/predCirc/motifs.txt',header = T,sep='\t')
colnames(rbp_db_filt)<-colnames(motifs)
write.table(rbp_db_filt,'./output/RBPpred/predCirc/motifs.txt',quote = F,row.names = F,sep='\t')

ppms<-read.table('reference/ATtRACT/pwm.txt',header=F, row.names=NULL,sep=' ')
query_ppm<-function(rbp_name,organism){
  organism='Homo_sapiens'
  
  rbp_id<-rbp_db[rbp_db$Gene_name==rbp_name&rbp_db$Organism==organism,]
  rbp_id<-rbp_id$Matrix_id[which(rbp_id$Len==max(rbp_id$Len))][1]
  print(rbp_id)
  
  idx<-which(ppms$V1==ppms$V1[!is.na(str_match(ppms$V1,paste0(">",rbp_id)))])
  for (i in (idx+1):length(ppms$V1)){
    # i=1
    if (substr(ppms$V1[i],1,1) == ">") {
      
      break
    }
  }
  idx_end<-i-1
  rbp_seq<-ppms$V1[(idx+1):idx_end]
  
  pwm<-data.frame(strsplit(rbp_seq,'\t',fixed=T))
  pwm<-apply(pwm,2,as.numeric)
  rownames(pwm) <- c("A", "G", "C", "T")
  
  # pwm<-makePWM(pwm)
  # seqLogo(pwm)
  return(ggseqlogo(pwm))
}
plot_selogos<-function(rbp_names){
  P<-list()
  organism='Homo_sapiens'
  for (i in 1:length(rbp_names)){
    P[[i]]<-query_ppm(rbp_names[i],organism)+theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
  }
  P[[i]]<-P[[i]]+theme(
    axis.title.x = element_text(size = 18)
  )
  
  return(P)
}

#####################################
## Figure S15A
#####################################
# barplot ranked by count
mergedMotifsFTS_circ<-mergedMotifsFTS_circ[order(mergedMotifsFTS_circ$count,decreasing = T),]
p<-ggplot(mergedMotifsFTS_circ[1:10,])+geom_bar(aes(x=count,y=reorder(id,count)),stat = 'identity')+
  xlab('Count')+ylab('')+
  theme_classic()+theme(axis.title.x = element_text(size=18),
                   axis.title.y = element_text(size=18),
                   axis.text.x = element_text(angle=45,hjust=1,size=18),
                   axis.text.y = element_text(size=15),
                   legend.title = element_text(size=15),
                   legend.text = element_text(size=12))
p1<-plot_selogos(p$data$id)
p2<-ggarrange(plotlist=p1,ncol = 1,nrow = length(p1),align = 'v',heights=c(rep(10,length(p1))))
ggarrange(p2,p,align = 'h',widths=c(20,80),nrow = 1,ncol=2)


p <-
  plotMotifs(
    mergedMotifsFTS_circ,
    mergedMotifsRFTS_circ,
    nf1 = (meanlength)*length(na.omit(annotatedBSJs$lenCircRNA)), 
    nf2 = (meanlength)*length(na.omit(annotatedBSJs$lenCircRNA)),
    log2FC = 1.4,
    # n=10,
    removeNegLog2FC = T,
    df1Name = "allCircRNA",
    df2Name = "random circRNAs",
    angle = 45
  )
p
p1<-plot_selogos(rev(p[[1]]$data$id))
p2<-ggarrange(plotlist=p1,ncol = 1,nrow = length(p1),align = 'v',heights=c(rep(10,length(p1))))
p2
p3<-ggarrange(p[[1]]+xlab('')+
                theme(axis.title.x = element_text(size=18),
                       axis.title.y = element_text(size=18),
                       axis.text.x = element_text(angle=0,hjust=1,size=18),
                       axis.text.y = element_text(size=15),
                       legend.title = element_text(size=15),
                       legend.text = element_text(size=12)),
          p[[2]]+xlab('RBP')+
            theme(axis.title.x = element_text(size=18),
                       axis.title.y = element_text(size=18),
                       axis.text.x = element_text(angle=30,hjust=1,size=18),
                       axis.text.y = element_text(size=15),
                       legend.title = element_text(size=15),
                       legend.text = element_text(size=12)),
          
          labels = c("", ""),
          align = 'h',
          ncol = 2,
          nrow = 1)
ggarrange(p2,p3,align = 'h',widths=c(20,80),nrow = 1,ncol=2)

## Find miRNA binding sites within internal circRNA seq (all circRNAs)
setwd('./output/RBPpred/predCirc')
miRsites_circ <- read.table('./output/RBPpred/predCirc/all_miRNA_binding_sites_internal_by_targetscan.txt',header=T,sep='\t')
head(miRsites_circ)
site_type<-unique(miRsites_circ$Site_type)
# pie chart of different sites_type
library(ggsci)
library(ggrepel)
df_show<-miRsites_circ
df_show$count<-1
#####################################
## Figure S16A
#####################################
ggplot()+geom_bar(df_show,aes(x='',fill=Site_type),stat = 'count',position='stack')+coord_polar(theta = 'y')+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(12)))+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.text = element_text(face = 'italic',size = 12),
        legend.title = element_text(face = 'italic',size = 15),
        legend.position = 'bottom',
        axis.title = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size=15))
# distribution of miRNA binding sites in all circRNAs
miRNA_count<-df_show

df_show$a_Gene_ID<-unlist(lapply(df_show$a_Gene_ID,function(x){strsplit(x,':',fixed = T)[[1]][1]}))

df_show<-dcast(df_show, a_Gene_ID~Site_type, fun.aggregate = sum,	
               subset = NULL, drop = TRUE, value.var = 'count')	
df_show$total<-rowSums(df_show[,2:5])
df_show<-df_show[order(df_show$total,decreasing = T),]
df_show<-subset(df_show,select = -c(total))

df_show_filt<-melt(df_show[1:10,],id.vars = 'a_Gene_ID',variable.name = 'sites_type',value.name = 'number')
df_show_filt$a_Gene_ID<-factor(df_show_filt$a_Gene_ID,levels=df_show[1:10,]$a_Gene_ID)
#####################################
## Figure S16B
#####################################
ggplot(df_show_filt)+geom_col(aes(x=a_Gene_ID,y=number,fill=sites_type),position = 'stack')+
  scale_fill_manual(values = c(pal_npg("nrc")(10),pal_aaas("default")(12)))+
  xlab("Host genes")+ylab("Number")+
  theme_bw()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=15),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))

############################################################################################################
#### Top 10 GO enrichment analysis for circRNA host genes with abundant miRNA binding sites
top10_host_genes<-c(unique(df_show_filt$a_Gene_ID))
# convert gene name to entrez gene id
# top10_host_gene_id_df <- bitr(top10_host_genes, fromType = "SYMBOL",
#                 toType = c("ENTREZID"),
#                 OrgDb = org.Bt.eg.db)
top10_host_gene_id_df<-getBM(attributes = c('external_gene_name','entrezgene_id'),filters = c('external_gene_name'),
      values = list(external_gene_name=top10_host_genes),mart = ensembl)
top10_host_geneid<-as.character(top10_host_gene_id_df$entrezgene_id)
gene <-  mapIds(org.Bt.eg.db, top10_host_gene_id_df$external_gene_name, 'ENTREZID', 'SYMBOL')    


#####################################
## Figure S24
#####################################
#GO enrichment analysis
ego_BP <- enrichGO(
  gene  = top10_host_geneid,
  keyType = "ENTREZID",
  OrgDb   = org.Bt.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  # qvalueCutoff  = 0.05,
  readable      = TRUE) #GO enrichment analysis

BP <- subset(ego_BP@result,select=c('Description','qvalue'))
BP <- BP[1:10,]

ego_MF <- enrichGO(
  gene  = top10_host_geneid,
  keyType = "ENTREZID",
  OrgDb   = org.Bt.eg.db,
  ont     = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE) #GO enrichment analysis

MF <- subset(ego_MF@result,select=c('Description','qvalue'))
MF <- MF[1:10,]
ego_CC <- enrichGO(
  gene  = top10_host_geneid,
  keyType = "ENTREZID",
  OrgDb   = org.Bt.eg.db,
  ont     = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE) #GO enrichment analysis

CC <- subset(ego_CC@result,select=c('Description','qvalue'))
CC <- CC[1:10,]

# GO enrichment result of Top10
GO_result<-rbind(BP,MF,CC)
GO_result$GO_term<-rep(c('BP','MF','CC'),rep(10,3))
GO_result$GO_term_1<-rep(c('BP','MF','CC'),rep(10,3))
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

## miRNA binding predict
miRbind_data<-read.table('./output/RBPpred/predCirc/miRNA_bySW.txt',sep='\t',header = T)
miRsites_circ$miRNA_name<-paste0('bta-',miRsites_circ$miRNA_family_ID)
colnames(miRbind_data)<-c('miRNA_name','end2','subseqDP','a_Gene_ID','mfe','alignedScore')
miRbind_data_filt<-merge(miRbind_data,miRsites_circ[,c('a_Gene_ID','miRNA_name','Site_type','UTR_start','UTR_end')],
                         by=c('miRNA_name','a_Gene_ID'))
df_show<-melt(miRbind_data_filt,id.vars = c(colnames(miRbind_data_filt)[1:4],'Site_type','UTR_start','UTR_end'),
              variable.name = 'Index',value.name = 'Score')
# df_show_filt<-dcast(df_show, miRNA_name~+Index, mean,value.var = 'Score')
df_show$Score<-as.numeric(df_show$Score)
# ggstatsplot::ggbetweenstats(df_show,aes(x=Site_type,y=Score))+facet_wrap(~Index)
#####################################
## Figure S16D
#####################################
ggplot(df_show)+geom_boxplot(aes(x=Site_type,y=Score,fill=Site_type))+facet_wrap(~Index,scales = "free_y")+
  theme_bw()+
  stat_compare_means(aes(x=Site_type,y=Score),bracket.size='20',paired=T,method='anova',size = 5)+
  scale_fill_manual(values = c(pal_npg("nrc")(10),pal_aaas("default")(12)))+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=15),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        
        legend.text = element_text(size=12),
        strip.text = element_text(size=18)
        )

head(miRbind_data_filt)
dim(miRbind_data_filt)
miRbind_data_filt<-read.table('./output/RBPpred/predCirc/miRbind_data_filt.txt',header = T,sep = '\t')
## Distribution of different breeds of RBP and miRNA binding sites
# RBP motifs
circRNA_expr<-read.table('./output/ABS_event/ABS_events_readscount.txt',header=T,sep = '\t')
circRNA_expr$circ<-paste0(circRNA_expr$chrom,':',circRNA_expr$start,'-',circRNA_expr$end)
# import all snps data and circQTL, GWAS-loci related circQTLs
# all snps
all_snps<-read.table('./output/combined/combined.vcf',header = T)
head(all_snps)
# import the correct ref and alt Alleles 
all_snps_ref_alt<-read.table('./output/combined/combined_data_ref_alt.vcf',sep="\t",header=F,comment.char = "#")
colnames(all_snps_ref_alt)[1:5]<-c("Chrom","ARS_UCD1_2_Pos","VariantID","Ref","Alt")
all_snps_ref_alt<-all_snps_ref_alt[,1:5]
all_snps_ref_alt$ARS_UCD1_2_Pos<-paste0(all_snps_ref_alt$Chrom,":",all_snps_ref_alt$ARS_UCD1_2_Pos)
# head(all_snps_ref_alt)
# circQTLs
circQTLs<-read.table('./output/circQTLs/high_vaild_significant_log_quant.txt',header = T)
# head(circQTLs)

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

############################################################################################################
## Distribution of all SNP and circQTL, GWAS-related circQTLs within RBP or miRNA binding sites
#RBP binding sites
library(zoo)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
setwd('./output/RBPpred/predCirc')
txdb<-makeTxDbFromGFF('./reference/genome/Bos_taurus.ARS-UCD1.2.101.gtf',format = 'gtf')

# get the location of RBP binding sites
# allRBPlist-->motifsFTS_circ$circ$locations
allexons<-exons(txdb)
rbp_binding_loc<-motifsFTS_circ$circ$locations
getMotifLocation<-function(id,pos,annotatedBSJs,allexons,type){
  
  #x:list of start locations for one Motif in specific circRNA 
  #id:"CREBRF:+:chr20:4799730:4804703"
  #pos:'1013, 1178, 1430, 1431, 1432, 1433' or NA
  # get the startUpBSE and endDownBSE according to the id
  
  # id="USP16:-:chr1:7149994:7145494";pos="48"
  # pos<-as.numeric(strsplit(rbp_binding_loc$AAAAAA[rbp_binding_loc$AAAAAA!="NA"][955],',')[[1]])
  annotatedBSJs<-annotatedBSJs;allexons<-allexons
  if (is.na(pos)){
    return(NA)
  }
  if (type=='RBP'){
    pos<-as.numeric(strsplit(pos,',')[[1]])
    # length<-6
  } else {
    temp<-strsplit(id,"|",fixed = T)[[1]]
    id<-temp[1]
    pos<-as.numeric(temp[2])
    # length<-pos[2]-pos[1]
  }
  
  chr<-str_sub(annotatedBSJs$chrom[annotatedBSJs$id==id],4,str_count(annotatedBSJs$chrom[annotatedBSJs$id==id]))
  startUpBSE<-annotatedBSJs$startUpBSE[annotatedBSJs$id==id]
  endDownBSE<-annotatedBSJs$endDownBSE[annotatedBSJs$id==id]
  strand<-annotatedBSJs$strand[annotatedBSJs$id==id]
  
  if (strand=="-"){
    id_exon_loc<-allexons[seqnames(allexons)==chr&start(allexons)<=startUpBSE&end(allexons)>=endDownBSE]
    id_exon_loc$width_sum<-rev(cumsum(rev(width(id_exon_loc))))
    motif_loc<-unlist(lapply(pos,function(x){poss_idx<-which(id_exon_loc$width_sum>=x)
    if (length(poss_idx)==0){
      return(NA)
    } else {
      idx<-max(poss_idx)
      return(end(id_exon_loc[idx])-x+1)}}))
  }else{
    id_exon_loc<-allexons[seqnames(allexons)==chr&start(allexons)>=startUpBSE&end(allexons)<=endDownBSE]
    id_exon_loc$width_sum<-cumsum(width(id_exon_loc))
    motif_loc<-unlist(lapply(pos,function(x){poss_idx<-which(id_exon_loc$width_sum>=x)
    if (length(poss_idx)==0){
      return(NA)
    } else {
      idx<-min(poss_idx)
      return(start(id_exon_loc[idx])+x-1)}}))
  }
  
  
}
# motif_loc_data<-data.frame()
motif_loc_data<-read.table('./output/RBPpred/predCirc/motif_loc_data.txt',header=T,sep='\t')
loc_data<-list()

for (i in 1:length(motifsFTS_circ$circ$locations$id)){
  # i=1
  id<-motifsFTS_circ$circ$locations$id[i]
  motif_row<-unlist(parLapply(cl,motifsFTS_circ$circ$locations[i,2:dim(motifsFTS_circ$circ$locations)[2]], 
                    function(x,y,z,d,e){library(stringr);getMotifLocation<-z;annotatedBSJs<-d;allexons<-e;
                      return(getMotifLocation(y,x,annotatedBSJs,allexons,'RBP'))},
                    y=id,z=getMotifLocation,d=annotatedBSJs,e=allexons))
  motif_row<-na.omit(motif_row)
  chrom<-strsplit(id,':')[[1]][3]
  strand<-strsplit(id,':')[[1]][2]
  if (length(motif_row)==0){
    print('next')
    next
  }
  # loc_data[[i]]<-data.frame('chrom'=id,'start'=motif_row,'end'=motif_row+6,'strand'=strand,'circRNA'=id,'motif'=names(motif_row))
  loc_data[[i]]<-data.frame('chrom'=id,'start'=motif_row,'end'=motif_row+6,'strand'=strand,'motif'=names(motif_row))
  
  if (i ==1){
    motif_loc_data<-loc_data[[i]]
  }else{
    motif_loc_data<-rbind(motif_loc_data,loc_data[[i]])
  }
}
write.table(motif_loc_data,'./output/RBPpred/predCirc/motif_loc_data.txt',row.names = F,quote = F,sep='\t')
save(i,file='i.Rdata')
motif_loc_data$chrom<-unlist(parLapply(cl,motif_loc_data$chrom,function(x){strsplit(x,":")[[1]][3]}))
ss<-motif_loc_data

# Total circQTLs or Filt the splice sites
## covert to GRanges object
index<-str_detect(ss$Chrom,'NKLS')
ss_5p<-ss
ss_5p$center<-rowSums(ss_5p[,c('start','end')])/2
ss_5p$start<-ss_5p$center-flanking
ss_5p$end<-ss_5p$center+flanking
ss_5p<-GRanges(Rle(ss_5p[,1]),IRanges(start=ss_5p$start,end=ss_5p$end),Rle(strand(ss_5p$strand)))

snps<-list(all_snps,circQTLs,gwas_circQTL)
labels<-c("Total snps","circQTLs","GWAS-related circQTLs")

## get related exon info of each circRNA
unique_circRNA<-unique(annotatedBSJs$id)
circRNA_exons<-read.table('./output/RBPpred/predCirc/circRNA_exons.txt',header = T,sep='\t')
circRNA_exons<-data.frame()
for (i in 1:length(unique_circRNA)){
  # i=1
  if (i %%100==0){
    print(i)
  }
  id=unique_circRNA[i]
  chr<-str_sub(annotatedBSJs$chrom[annotatedBSJs$id==id],4,str_count(annotatedBSJs$chrom[annotatedBSJs$id==id]))
  startUpBSE<-annotatedBSJs$startUpBSE[annotatedBSJs$id==id]
  endDownBSE<-annotatedBSJs$endDownBSE[annotatedBSJs$id==id]
  strand<-annotatedBSJs$strand[annotatedBSJs$id==id]
  if (strand=="-"){
    id_exon_loc<-allexons[seqnames(allexons)==chr&start(allexons)<=startUpBSE&end(allexons)>=endDownBSE]
    id_exon_loc$width_sum<-rev(cumsum(rev(width(id_exon_loc))))
  } else {
    id_exon_loc<-allexons[seqnames(allexons)==chr&start(allexons)>=startUpBSE&end(allexons)<=endDownBSE]
    id_exon_loc$width_sum<-cumsum(width(id_exon_loc))
  }
    temp<-as.data.frame(id_exon_loc)
    temp$id<-id
  if (i==1){
    circRNA_exons<-temp
  }
    circRNA_exons<-rbind(circRNA_exons,temp)
}
write.table(circRNA_exons,'./output/RBPpred/predCirc/circRNA_exons.txt',quote = F,row.names = F,sep='\t')

## get the all circQTLs within circRNA
getQTLinRNA<-function(snps,RNAs,RNA_index,window_size){
  snp<-snps[[1]]
  snp<-merge(snp,all_snps_ref_alt,by=c("Chrom","ARS_UCD1_2_Pos","VariantID"))
  snp<-snp[,c("Chrom","ARS_UCD1_2_Pos","VariantID","Ref","Alt")]
  # snp$Chrom<-paste0("chr",snp$Chrom)
  # snp$Pos<-lapply(snp$ARS_UCD1_2_Pos,function(x){as.numeric(strsplit(x,':',fixed=T)[[1]][2])})
  ##
  #(Genomic location): chrom   start     end strand  a_Gene_ID motif/miRNA
  RNAs<-miRNA_loc_data[,c("chrom","start","end","strand","a_Gene_ID","miRNA_name")]
  na_index<-!is.na(RNAs$start);RNAs<-RNAs[na_index,]
  
  RNAs<-RNAs[RNA_index,]
  window_size=10
  RNAs$RNAs_loc<-paste0(RNAs$chrom,":",(RNAs$start+RNAs$end)/2-window_size,":",(RNAs$start+RNAs$end)/2+window_size)
  # RNAs$RNAs_loc<-paste0(RNAs$chrom,":",RNAs$start,":",RNAs$end)
  isInRNAs<-function(RNA_loc,snp_loc){
    snp_chr<-paste0("chr",strsplit(snp_loc,":")[[1]][1])
    snp_pos<-as.numeric(strsplit(snp_loc,":")[[1]][2])
    RNA_chr<-strsplit(RNA_loc,":")[[1]][1]
    temp<-c(as.numeric(strsplit(RNA_loc,":")[[1]][2]),as.numeric(strsplit(RNA_loc,":")[[1]][3]))
    RNA_start<-min(temp)
    RNA_end<-max(temp)
    if (snp_chr==RNA_chr){
        res<-snp_pos>=RNA_start&snp_pos<=RNA_end
      if (res){
        loc<-snp_pos-RNA_start+1
        # print('ok')
        return(loc)
      }else{
        return(NA)
      }
      
    } else{
      return(NA)
    }
  
  }
  ## isInExon
  
  ## isInCircRNAs
  isIncircRNAs<-function(circRNA_id,snp_loc,circRNA_exons){
    ## return location in transcript of circRNA (not genomic location)
    isInExon<-function(snp_pos,exon_info){
      # snp_pos<-2881168
      exon_info<-exon_info[!duplicated(exon_info[,c("start","end")]),]
      res<-which(exon_info$start<=snp_pos&exon_info$end>=snp_pos)
      if (length(res)!=0){
        if (exon_info$strand[1]=="+"){
          if (res==1){former=0}else{former=sum(exon_info$width[1:(res-1)])}
          loc<-former+snp_pos-exon_info$start[res]+1
          return(paste0(loc,collapse = ','))
        } else{
          if (res==length(exon_info$seqnames)){latter=0}else{latter=sum(exon_info$width[(res+1):3])}
          loc<-latter+exon_info$end[res]-snp_pos+1
          return(paste0(loc,collapse = ','))
        }
      } else{
        return(NA)
      }
      
    }
    
    exon_info<-circRNA_exons[circRNA_exons$id==circRNA_id,]
    exon_info<-exon_info[!duplicated(exon_info[,c("start","end")]),]
    
    # RNA_loc<-strsplit(circRNA_id,'')[[1]][3]
    snp_chr<-paste0("chr",strsplit(snp_loc,":")[[1]][1])
    snp_pos<-as.numeric(strsplit(snp_loc,":")[[1]][2])
    # RNA_loc<-RNAs$RNAs_loc[1]
    RNA_chr<-strsplit(circRNA_id,":")[[1]][3]
    temp<-c(as.numeric(strsplit(circRNA_id,":")[[1]][4]),as.numeric(strsplit(circRNA_id,":")[[1]][5]))
    RNA_start<-min(temp)
    RNA_end<-max(temp)
    RNA_strand<-strsplit(circRNA_id,":")[[1]][2]
    if (snp_chr==RNA_chr){
      res<-snp_pos>=RNA_start&snp_pos<=RNA_end
      if (res){
        
        loc<-isInExon(snp_pos,exon_info)
        # print('ok')
        return(loc)
      }else{
        return(NA)
      }
      
    } else{
      return(NA)
    }
    
  }
  
  ## snps;for loop of snps
  isIn_foreach<-function(x){
    a=unlist(parLapply(cl,RNAs$RNAs_loc,isInRNAs,x))
    return(a)
  }
  
  # res = lapply(snp$ARS_UCD1_2_Pos,isIn_foreach)
  unique_circRNA<-annotatedBSJs$id
  res<-read.table('./output/RBPpred/predCirc/all_snps_in_circRNA.txt',sep='\t',header = T)
  # res<-data.frame('id'=unique_circRNA)
  t1=Sys.time()
  for (i in 9501:length(snp$ARS_UCD1_2_Pos)) {
    if (i %% 100==0){
      print(i)
      if (i %% 1000==0){
        write.table(res,'./output/RBPpred/predCirc/all_snps_in_circRNA.txt',quote = F,row.names = F,sep='\t')
        save(i,file='i.Rdata')
        print("writted to file!")
      }
    }
    
    find_snp=unlist(parLapply(cl,unique_circRNA,function(x,y,z,de){isIncircRNAs<-de;snp_loc<-y;circRNA_exons<-z
                        return(isIncircRNAs(x,snp_loc,circRNA_exons))},
                       y=snp$ARS_UCD1_2_Pos[i],z=circRNA_exons[,-c(6:7)],de=isIncircRNAs))
    
    
    b=sum(!is.na(find_snp))
    if (b!=0){
      print("Find it!")
      res[,snp$VariantID[i]]=find_snp
    }
    
  }
  t2=Sys.time()
  print(t2-t1)
  write.table(res,'./output/RBPpred/predCirc/all_snps_in_circRNA.txt',quote = F,row.names = F,sep='\t')
}
## Calculate the tagMatrix
distri_bins_res<-list()
distri_denstiy_res<-list()
for (i in 1:length(snps)){
  circqtl_filt<-snps[[i]]
  # circqtl_filt<-all_snps
  circqtl_filt<-subset(circqtl_filt,select = c(VariantID,ARS_UCD1_2_Pos))
  label<-labels[i]
  circqtl_filt$Chrom<-apply(X=data.frame(circqtl_filt$ARS_UCD1_2_Pos),1,FUN=function(x){paste0('chr',strsplit(x,':',fixed=T)[[1]][1])})
  circqtl_filt$Pos<-apply(X=data.frame(circqtl_filt$ARS_UCD1_2_Pos),1,FUN=function(x){strsplit(x,':',fixed=T)[[1]][2]})
  circqtl_filt<-circqtl_filt[,-2]
  peak<-circqtl_filt[,-1]
  peak<-peak[!is.na(peak$Pos),]
  ## config the flanking
  flanking<-100
  
  peak<-GRanges(Rle(peak[,1]),IRanges(peak[,2]))
  tagmatrix<-getTagMatrix(peak, weightCol = NULL, windows=ss_5p, flip_minor_strand = TRUE)
  # RNA_index<-as.numeric(rownames(tagmatrix))
  # ss_5p[RNA_index,]
  distri_ss5p<-colSums(tagmatrix)
  
  # Distribution of snps in bins
  bins=10
  xaxis<-seq(-flanking,flanking,bins)
  g<-rep(1:((length(distri_ss5p)-1)/bins),each=bins)
  distri_ss5p_sum<-tapply(distri_ss5p[-1],g,sum)
  # distri_ss3p_sum<-tapply(distri_ss3p[-1],g,sum)
  ### bins!!!
  distri_data_bins<-data.frame(t(rbind(c(distri_ss5p_sum),xaxis)))
  ## rollsum!!
  # distri_data_bins<-data.frame(t(rbind(c(rollsum(distri_ss5p,window_size)),
  # seq(-flanking+window_size/2,flanking-window_size/2+1,1))))
  
  colnames(distri_data_bins)<-c('Number','Distance')
  distri_data_bins$QTL<-labels[i]
  distri_data_bins$Splice_type<-'Upstream'
  distri_data_bins$Splice_type[distri_data_bins$Distance<0]<-'Downstream'
  
  distri_bins_res[[i]]<-distri_data_bins
  
  # Distribution of snps in sliding windowsize
  window_size=10
  distri_data<-data.frame(t(rbind(c(rollmean(distri_ss5p,window_size)
  ),
  seq(-flanking+window_size/2,flanking-window_size/2+1,1))))
  colnames(distri_data)<-c('value','Distance')
 
  distri_data$value<-scale(distri_data$value,center = T,scale = T)
  distri_data$QTL<-labels[i]
  
  distri_denstiy_res[[i]]<-distri_data
  
  if (i == 1){
    df_distri_bins<-distri_bins_res[[i]]
    df_distri_density<-distri_denstiy_res[[i]]
  }else{
    df_distri_bins<-rbind(df_distri_bins,distri_bins_res[[i]])
    df_distri_density<-rbind(df_distri_density,distri_denstiy_res[[i]])
  }
  
}

#####################################
## Figure S17A/B
#####################################
ggplot(df_distri_bins)+geom_histogram(aes(x=Distance,y=Number,fill=QTL),stat = 'identity',position = "stack")+
  theme_bw()+xlab('Distance to RBP binding sites (bp)')+ylab('Number per window size')+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))+
  # annotate("text", x = -600, y = -0.5, label = "Exons",size=5)+
  # annotate("text", x = 0, y = -0.5, label = "Introns",size=5)+
  # annotate("text", x = 600, y = -0.5, label = "Exons",size=5)+
  geom_vline(xintercept=c(-flanking,flanking), linetype="dotted")
# density distribution of different circQTLs to the RBP binding sites
#####################################
## Figure S17C
#####################################
ggplot(df_distri_density,aes(x=Distance,y=value,color=QTL))+
  # geom_point()+
  geom_line(size=1)+ylim(-2.5,2.5)+
  scale_y_continuous(expand = c(0,0))+
  # geom_smooth(method = 'loess',span=0.5)+
  theme_bw()+xlab('Distance to RBP binding sites (bp)')+ylab('Z-score of Mean per window size (20bp)')+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))+
  # annotate("text", x = -600, y = -0.5, label = "Exons",size=5)+
  # annotate("text", x = 0, y = -0.5, label = "Introns",size=5)+
  # annotate("text", x = 600, y = -0.5, label = "Exons",size=5)+
  geom_vline(xintercept=c(-flanking,flanking), linetype="dotted")+
# annotate("text", x = 350, y = -1.5, label = "Exons")+
# annotate("text", x = 820, y = -1.5, label = "Introns")+
annotate("rect", xmin = -6-window_size/2 , xmax = 6+window_size/2, ymin = -2.5, ymax = 2.5, alpha = .5, fill = "lightblue")
# annotate("rect", xmin = -250, xmax = 0, ymin = -1.0, ymax = -1.3, alpha = .2, colour = "black")+
# annotate("rect", xmin = 0, xmax = 700, ymin = -1.1, ymax = -1.2, alpha = .99, colour = "black")+
# annotate("rect", xmin = 700, xmax = 950, ymin = -1.0, ymax = -1.3, alpha = .2, colour = "black")

## MiRNA bindings sites
# miRbind_data_filt$
miRbind_data_filt$binding_gstar=unlist(parLapply(cl,paste0(miRbind_data_filt$a_Gene_ID,"|",miRbind_data_filt$UTR_start), 
                   function(x,y,z,d,e){library(stringr);getMotifLocation<-z;annotatedBSJs<-d;allexons<-e;
                   return(getMotifLocation(x,y,annotatedBSJs,allexons,'miRNA'))},
                   y="",z=getMotifLocation,d=annotatedBSJs,e=allexons))

miRbind_data_filt$chrom<-unlist(lapply(miRbind_data_filt$a_Gene_ID,function(x){strsplit(x,":",fixed = T)[[1]][3]}))
miRbind_data_filt$start<-miRbind_data_filt$binding_gstar
miRbind_data_filt$end<-miRbind_data_filt$binding_gstar+(miRbind_data_filt$UTR_end-miRbind_data_filt$UTR_start)
miRbind_data_filt$strand<-unlist(lapply(miRbind_data_filt$a_Gene_ID,function(x){strsplit(x,":",fixed = T)[[1]][2]}))

miRbind_data_filt$length<-miRbind_data_filt$UTR_end-miRbind_data_filt$UTR_start
miRbind_data_filt$end[miRbind_data_filt$strand=="-"]<-miRbind_data_filt$start[miRbind_data_filt$strand=="-"]
miRbind_data_filt$start[miRbind_data_filt$strand=="-"]<-miRbind_data_filt$end[miRbind_data_filt$strand=="-"]-
  miRbind_data_filt$length[miRbind_data_filt$strand=="-"]
write.table(miRbind_data_filt,'./output/RBPpred/predCirc/miRbind_data_filt.txt',row.names = F,quote = F,sep="\t")
miRNA_loc_data<-miRbind_data_filt[,c("chrom","start","end","strand","a_Gene_ID","miRNA_name","Site_type")]

write.table(miRNA_loc_data,'./output/RBPpred/predCirc/miRNA_loc_data.txt',row.names = F,quote = F,sep='\t')

miRNA_loc_data<-read.table('./output/RBPpred/predCirc/miRNA_loc_data.txt',header = T)
# miRNA_genome_loc<-read.table('./reference/miRNA/bta_miRNA.bed',sep = '\t')
# colnames(miRNA_genome_loc)<-c("miR_chrom","miR_start",'miR_end','miR_strand',"miRNA_name")
# miRNA_loc_data<-merge(miRNA_loc_data,miRNA_genome_loc,by='miRNA_name')
ss<-miRNA_loc_data

# Total circQTLs or Filt the splice sites
## covert to GRanges object
index<-str_detect(ss$Chrom,'NKLS')
na_index<-is.na(ss$start)
ss_5p<-ss[!na_index,]
ss_5p$center<-rowSums(ss_5p[,c('start','end')])/2
ss_5p$start<-ss_5p$center-flanking
ss_5p$end<-ss_5p$center+flanking
ss_5p<-GRanges(Rle(ss_5p[,1]),IRanges(start=ss_5p$start,end=ss_5p$end),Rle(strand(ss_5p$strand)))

#####################################
## Figure S17D/E
#####################################
# distribution of miRNA binding sites 
ggplot(df_distri_bins)+geom_histogram(aes(x=Distance,y=Number,fill=QTL),stat = 'identity',position = "stack")+
  theme_bw()+xlab('Distance to miRNA binding sites (bp)')+ylab('Number per window size')+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))+
  # annotate("text", x = -600, y = -0.5, label = "Exons",size=5)+
  # annotate("text", x = 0, y = -0.5, label = "Introns",size=5)+
  # annotate("text", x = 600, y = -0.5, label = "Exons",size=5)+
  geom_vline(xintercept=c(-flanking,flanking), linetype="dotted")

#####################################
## Figure S17F
#####################################
# density distribution of different circQTLs to the RBP binding sites
ggplot(df_distri_density,aes(x=Distance,y=value,color=QTL))+
  # geom_point()+
  geom_line(size=1)+
  scale_y_continuous(expand = c(0,0))+
  # geom_smooth(method = 'loess',span=0.5)+
  theme_bw()+xlab('Distance to RBP binding sites (bp)')+ylab('Z-score of Mean per window size (20bp)')+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))+
  # annotate("text", x = -600, y = -0.5, label = "Exons",size=5)+
  # annotate("text", x = 0, y = -0.5, label = "Introns",size=5)+
  # annotate("text", x = 600, y = -0.5, label = "Exons",size=5)+
  geom_vline(xintercept=c(-flanking,flanking), linetype="dotted")+
  # annotate("text", x = 350, y = -1.5, label = "Exons")+
  # annotate("text", x = 820, y = -1.5, label = "Introns")+
  annotate("rect", xmin = -6-window_size/2 , xmax = 6+window_size/2, ymin = -2.5, ymax = 5, alpha = .5, fill = "lightblue")
# annotate("rect", xmin = -250, xmax = 0, ymin = -1.0, ymax = -1.3, alpha = .2, colour = "black")+
# annotate("rect", xmin = 0, xmax = 700, ymin = -1.1, ymax = -1.2, alpha = .99, colour = "black")+
# annotate("rect", xmin = 700, xmax = 950, ymin = -1.0, ymax = -1.3, alpha = .2, colour = "black")

############################################################################################################
## Venn plot of RBPs and miRNA binding sites among breeds
# import breed info
# breed info
breed_info = read.csv('./SRR_breeds.csv')
# breed_info<-breed_info[unlist(lapply(breed_info$SRR_list,function(x){x%in%all_gtf})),]
breeds = breed_info$Breeds[!duplicated(breed_info$Breeds)]
col=breed_info$SRR_list
## miRNA binding sites
split_loc_circprofile<-function(x){
  chr<-strsplit(x,':',fixed = T)[[1]][3]
  start<-strsplit(x,':',fixed = T)[[1]][4]
  end<-strsplit(x,':',fixed = T)[[1]][5]
  loc<-as.numeric(c(start,end))
  return(paste0(chr,':',min(loc),'-',max(loc)))
}
miRNA_count$circ<-unlist(lapply(miRNA_count$a_Gene_ID,split_loc_circprofile))
miRNA_count$a_Gene_ID<-unlist(lapply(df_show$a_Gene_ID,function(x){strsplit(x,':',fixed = T)[[1]][1]}))

miRNA_count<-dcast(miRNA_count, circ~Site_type, fun.aggregate = sum,	
               subset = NULL, drop = TRUE, value.var = 'count')
idx=unlist(lapply(miRNA_count$circ,function(x){which(circRNA_expr$circ==x)}))
circRNA_expr_miR<-circRNA_expr[idx,]
aidx=unlist(lapply(circRNA_expr_miR$circ,function(x){which(miRNA_count$circ==x)}))
miRNA_count<-miRNA_count[aidx,]
miRNA_count_matrix<-crossprod(as.matrix(miRNA_count[,2:5]),as.matrix(circRNA_expr_miR[,breed_info$SRR_list]!=0))
# merge breed by mean
miRNA_count_convert<-data.frame('Site_type'=rownames(miRNA_count_matrix))
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  miRNA_count_breed = subset(miRNA_count_matrix,select = srr)
  miRNA_count_convert[,as.character(breeds[i])]=apply(miRNA_count_breed,1,mean)/sum(apply(miRNA_count_breed,1,mean))
  
}
df_show<-melt(miRNA_count_convert,id.vars = 'Site_type',variable.name = 'breeds',value.name = 'count')
#####################################
## Figure S16C
#####################################
ggplot(df_show)+geom_bar(aes(x=breeds,y=count,fill=Site_type),stat = 'identity',position='stack')+
  scale_fill_manual(values=c(pal_npg("nrc")(10),pal_aaas("default")(12)))+
  ylab('Frequency')+theme_bw()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=15),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))
# merge miRNA and circexpr by circ
head(circRNA_expr)
head(miRbind_data_filt)

## Distribution of RBPs and miRNA binding sites among different breeds
circRNA_expr_convert<-data.frame('circ'=circRNA_expr[,'circ'])
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  circRNA_expr_breed = subset(circRNA_expr,select = srr)
  circRNA_expr_convert[,as.character(breeds[i])]=apply(circRNA_expr_breed,1,mean)
}
## diff in RBP
rbps<-motifsFTS_circ$circ$motifs
rbp_count<-motifsFTS_circ$circ$counts
colnames(rbp_count)<-c('id',c(rbps$id))
rbp_count$circ<-unlist(lapply(rbp_count$id,function(x){chr<-strsplit(x,':',fixed=T)[[1]][3];
  start<-strsplit(x,':',fixed=T)[[1]][4];end<-strsplit(x,':',fixed=T)[[1]][5];
  return(paste0(chr,':',start,'-',end))
  }))


## Meanlength for each motif 
rbp_count$Motifsum<-rowSums(rbp_count[,c(rbps$id)])
rbp_count_circ<-merge(rbp_count[,c('circ','Motifsum')],circRNA_expr,by='circ')
# rbp_count_circ[,breed_info$SRR_list]<-rbp_count_circ$Motifsum*(rbp_count_circ[,breed_info$SRR_list])
rbp_count_circ[,breed_info$SRR_list]<-rbp_count_circ$Motifsum*(rbp_count_circ[,breed_info$SRR_list])

for (i in 1:length(breed_info$SRR_list)) {
  rbp_count_circ[,breed_info$SRR_list[i]]<-rbp_count_circ[,breed_info$SRR_list[i]]/(
  meanlength*(sum(rbp_count_circ[,breed_info$SRR_lisy[i]]!=0)+1))
}
annotation_col<-breed_info
rownames(annotation_col)<-annotation_col$SRR_list
annotation_col<-subset(annotation_col,select = 'Breeds')
group=unique(annotation_col$Breeds)
names(group)<-c('#ff83ff','#ff9289','#d3ba00','#00d65c','#00dae0','#82b7ff')
annotation_color=list(
  group=group
)
#####################################
## Figure S15B-RBP(right)
#####################################
# annotation_col<-subset(annotation_col,select = 'Breeds');annotation_col['SRR87031987',]<-'Hereford'
pheatmap(rbp_count_circ[,breed_info$SRR_list],scale='column',cluster_rows = T,cluster_cols = T,show_rownames = F,
         show_colnames = F,annotation_col = annotation_col,annotation_colors = annotation_color)
rbp_count_circ[,breed_info$SRR_list]<-apply(rbp_count_circ[,breed_info$SRR_list],2,scale)

rbp_count_circ_convert<-data.frame('circ'=rbp_count_circ[,'circ'])
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  rbp_count_circ_breed = subset(rbp_count_circ,select = srr)
  rbp_count_circ_convert[,as.character(breeds[i])]=apply(rbp_count_circ_breed,1,mean)
}
rbp_count_circ_convert<-melt(rbp_count_circ_convert,id.vars = 'circ',variable.name = 'breeds',value.name = 'Mean')
#####################################
## Figure S15C-motif (left)
#####################################
ggplot(rbp_count_circ_convert)+geom_violin(aes(x=breeds,y=Mean,fill=breeds))+
  geom_boxplot(aes(x=breeds,y=Mean,fill=breeds),width=0.2)+
  stat_compare_means(aes(x=breeds,y=Mean),method = 'anova',size=5)+
  theme_bw()+
  scale_fill_manual(values = names(annotation_color$group))+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=15),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)
  )
ggstatsplot::ggbetweenstats(
  data = rbp_count_circ_convert,
  x=breeds,y=Mean,fill=breeds,
  type='p',
  centrality.plotting=F,
  messages = F
)
#####################################
## Figure S15B-motif (left)
#####################################
## Get Matrix for each motif (dim of motif X dim of circRNA)
# order the column of circRNA expression matrix
rbp_count[,c('circ','Motifsum')]
rbp_count$circ_convert<-unlist(lapply(rbp_count$circ,function(x){
  chr<-strsplit(x,':',fixed=T)[[1]][1];
  loc<-as.numeric(strsplit(strsplit(x,':',fixed=T)[[1]][2],'-',fixed=T)[[1]]);
  start<-min(loc);end<-max(loc);
  return(paste0(chr,':',start,'-',end))
}))
idx=unlist(lapply(rbp_count$circ_convert,function(x){which(circRNA_expr$circ==x)}))
circRNA_expr<-circRNA_expr[idx,]
aidx=unlist(lapply(circRNA_expr$circ,function(x){which(rbp_count$circ_convert==x)}))
rbp_count<-rbp_count[aidx,]
## Mean length for each RBP
allRBPlist<-rbps$id
uniqueRBPlist<-unique(unlist(lapply(allRBPlist,function(x){strsplit(x,',',fixed=T)[[1]]})))
rbp_loc_count<-function(x,allRBPlist){
  loc<-!is.na(str_match(allRBPlist,x))
}
rbp_count_rbp<-data.frame('circ'=rbp_count$circ)
for (i in 1:length(uniqueRBPlist)) {
  loc<-rbp_loc_count(uniqueRBPlist[i],allRBPlist)
  if (sum(loc)==1){
    rbp_count_rbp[,uniqueRBPlist[i]]<-rbp_count[,allRBPlist[loc]]
  }else{
    rbp_count_rbp[,uniqueRBPlist[i]]<-rowSums(rbp_count[,allRBPlist[loc]])
  }
  
}
## inner factor
rbp_breed_matrix<-crossprod(as.matrix(rbp_count_rbp[,uniqueRBPlist]),as.matrix(circRNA_expr[,breed_info$SRR_list]!=0))

annotation_col<-breed_info
rownames(annotation_col)<-annotation_col$SRR_list
annotation_col<-subset(annotation_col,select = 'Breeds')
# annotation_col<-subset(annotation_col,select = 'Breeds');annotation_col['SRR87031987',]<-'Hereford'
pheatmap(rbp_breed_matrix,scale='column',border_color=NA,cluster_rows = T,cluster_cols = T,show_rownames = F,
         show_colnames = F,annotation_col = annotation_col,annotation_colors = annotation_color)
# Statistics of (RBPs count) heatmap among different breeds
colnames(rbp_breed_matrix)<-colnames(circRNA_expr[,breed_info$SRR_list])
rbp_breed_matrix<-apply(rbp_breed_matrix, 2,scale)
rbp_count_circ_convert<-data.frame('RBP'=uniqueRBPlist)
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  rbp_count_circ_breed = subset(rbp_breed_matrix,select = srr)
  rbp_count_circ_convert[,as.character(breeds[i])]=apply(rbp_count_circ_breed,1,mean)
}
rbp_count_circ_convert<-melt(rbp_count_circ_convert,id.vars = 'RBP',variable.name = 'breeds',value.name = 'Mean')
#####################################
## Figure S15C-RBP (right)
#####################################
ggplot(rbp_count_circ_convert)+
  geom_violin(aes(x=breeds,y=Mean,fill=breeds))+
  geom_boxplot(aes(x=breeds,y=Mean,fill=breeds),width=0.2)+
  # geom_jitter(aes(x=breeds,y=Mean),
  #             shape=16,size=1,position=position_jitter(0.2))+
  theme_bw()+
  stat_compare_means(aes(x=breeds,y=Mean),method = "anova",size=5)+
  # scale_color_manual(values = names(annotation_color$group))+
  scale_fill_manual(values = names(annotation_color$group))+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=15),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)
  )
ggstatsplot::ggbetweenstats(
  data = rbp_count_circ_convert,
  x=breeds,y=Mean,fill=breeds,
  type='p',
  centrality.plotting=F,
  messages = F
)+
  theme_bw()+
  scale_fill_manual(values = names(annotation_color$group))+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1,size=15),
        axis.text.y = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)
  )



