#!/bin/Rscript

library(ggplot2)
library(reshape2)
library(stringr)
library(igraph)
library(dplyr)
library(circRNAprofiler)
library(wordcloud2)

############################################################################################################
working_dir <- "path_to_your_working_dir"
circRNAprofiler_dir<-paste(working_dir,'output/RBPpred/predCirc',sep='/')
all_RNA_seq_file='./SRR_list_CIRI.txt'
output_dir <- paste(working_dir,'./output/ceRNA_network',sep='/')
setwd(circRNAprofiler_dir)
dir.create(output_dir)
# Detect cores and initialize that in local computer
no_cores <- 3
cl <- makeCluster(no_cores)


## import motif data
# setwd()
load('./motifsFTS_circ.Rdata')
load('./motifsFTS_circ_alt.Rdata')
load('./motifsRFTS_circ_random.Rdata')
load('./motifsRFTS_circ.Rdata')

mergedmotifsRFTS_circ<-mergeMotifs(motifsRFTS_circ)
mergedmotifsRFTS_circ$count<-mergedmotifsRFTS_circ$count/dim(motifsRFTS_circ$circ$targets)[1]/650.2962 # average lenght of circRNA is 650.2962

RBP_data<-motifsFTS_circ$circ$targets[,c("id","gene","length")]
RBP_data_count<-motifsFTS_circ$circ$counts
motif_data<-motifsFTS_circ$circ$motifs
# rm(motifsFTS_circ)
# rm(mergedmotifsRFTS_circ)
colnames(motif_data)<-c("motif","rbp")
RBP_data<-merge(RBP_data,RBP_data_count,by="id")
head(RBP_data)
all_rbp<-unique(unlist(strsplit(motif_data$rbp,",")))
RBP_stat_coount<-RBP_data[,c("id","gene","length")]
for (i in 1:length(all_rbp)){
   # i=1
   rbp<-all_rbp[i]
   rbp_random_count<-mergedmotifsRFTS_circ$count[mergedmotifsRFTS_circ$id==rbp]
   eachrbp_motifs<-motif_data$motif[!is.na(str_match(motif_data$rbp,rbp))]
   if (length(eachrbp_motifs)==1){
      RBP_stat_coount[,rbp]<-RBP_data[,eachrbp_motifs]/RBP_stat_coount[,"length"]/rbp_random_count
   } else{
      RBP_stat_coount[,rbp]<-rowSums(RBP_data[,eachrbp_motifs])/RBP_stat_coount[,"length"]/rbp_random_count
   }
}
RBP_stat_coount<-RBP_stat_coount[rowSums(RBP_stat_coount[,all_rbp])>0,]
## alt vs ref
RBP_stat_coount_ref<-RBP_stat_coount[unlist(lapply(alt_genotype_circseq$id[alt_genotype_circseq$isconvert=="alt"],function(x){which(RBP_stat_coount$id==x)})),]
RBP_stat_coount_ref<-subset(RBP_stat_coount_ref,select=-length)
raw_data_ref<-melt(RBP_stat_coount_ref,id.vars = c("id","gene"),
               variable.name = "RBP",value.name = "count")
raw_data_ref<-raw_data_ref[,c("id","gene","RBP","count")]
## wordCloud
raw_data_ref<-raw_data_ref[raw_data_ref$count>0,]
rbp_freq_ref<-data.frame("RBP"=raw_data_ref[,c("RBP")])
rbp_freq_ref$count<-1
rbp_freq_ref<-data.frame(rbp_freq_ref%>%group_by(RBP)%>%summarize(count=n()))
colnames(rbp_freq_ref)<-c('id','count')
##
RBP_stat_coount_alt<-RBP_stat_coount[unlist(lapply(alt_genotype_circseq$id[alt_genotype_circseq$isconvert=="alt"],function(x){which(RBP_stat_coount$id==x)})),]
RBP_stat_coount_alt<-subset(RBP_stat_coount_alt,select=-length)
raw_data_alt<-melt(RBP_stat_coount_alt,id.vars = c("id","gene"),
                   variable.name = "RBP",value.name = "count")
raw_data_alt<-raw_data_alt[,c("id","gene","RBP","count")]
## wordCloud
raw_data_alt<-raw_data_alt[raw_data_alt$count>0,]
rbp_freq_alt<-data.frame("RBP"=raw_data_alt[,c("RBP")])
rbp_freq_alt$count<-1
rbp_freq_alt<-data.frame(rbp_freq_alt%>%group_by(RBP)%>%summarize(count=n()))
colnames(rbp_freq_alt)<-c('id','count')
rbp_freq_alt$motif<-"AAAACC"
rbp_freq_ref$motif<-"AAAACC"
#####################################
## Figure S18A
#####################################
## logFC
p <-
   plotMotifs(
      
      rbp_freq_alt,
      rbp_freq_ref,
      nf1 = 1, 
      nf2 = 1,
      log2FC = 0.05,
      df1Name = "alt",
      df2Name = "ref",
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

library(ggsci)
raw_data_alt$geno<-"alt"
raw_data_ref$geno<-"ref"
rbp<-rbp_freq_alt$id[which(abs(log2(rbp_freq_alt$count/rbp_freq_ref$count))>0)]

rbp_ref_alt_data<-rbind(raw_data_ref,raw_data_alt)
idx<-unlist(lapply(rbp,function(x){which(rbp_ref_alt_data$RBP==x)}))
# 
rbp_ref_alt_data$geno<-factor(rbp_ref_alt_data$geno,levels = c('ref','alt'))
# shapiro.test(rbp_ref_alt_data$count[idx])
library(rstatix)
df_wilcox <- rbp_ref_alt_data[idx,] %>%group_by(RBP)%>%
   pairwise_wilcox_test(count ~ geno) %>%
   add_y_position(step.increase = 0.02)
df_wilcox_filt<-df_wilcox[df_wilcox$p.adj<=0.05,]
#####################################
## Figure S18B
#####################################
ggplot(rbp_ref_alt_data[idx,])+geom_boxplot(aes(x=RBP,y=count,fill=geno),width=0.5,vjust=1)+
   # geom_violin(aes(x=RBP,y=count,fill=geno))+
   theme_bw()+ylab("Normalized count")+
   scale_fill_manual(values = c(pal_npg("nrc")(10)[6:5]))+
   
   # stat_compare_means(aes(x=RBP,y=count),bracket.size='20',method='t.test',size = 5,hide.ns = F)+
   theme(axis.title.x = element_text(size=25),
         axis.title.y = element_text(size=25),
         axis.text.x = element_text(angle=30,hjust=1,size=18),
         axis.text.y = element_text(size=18),
         legend.title = element_text(size=18),
         legend.text = element_text(size=16))

#############################################################################################
RBP_stat_coount<-RBP_stat_coount[!is.na(RBP_stat_coount$id),]
RBP_stat_coount<-subset(RBP_stat_coount,select=-length)
raw_data<-melt(RBP_stat_coount,id.vars = c("id","gene"),
                variable.name = "RBP",value.name = "count")

raw_data<-raw_data[,c("id","gene","RBP","count")]
## wordCloud
raw_data<-raw_data[raw_data$count>0,]
rbp_freq<-data.frame("RBP"=raw_data[,c("RBP")])
rbp_freq$count<-1

rbp_freq<-data.frame(rbp_freq%>%group_by(RBP)%>%summarize(count=n()))
colnames(rbp_freq)<-c('word','freq')
rownames(rbp_freq)<-rbp_freq$word
rbp_freq<-rbp_freq[order(rbp_freq$freq,decreasing = T),]
#####################################
## Figure 7B
#####################################
wordcloud2(rbp_freq)

raw_data_filt<-raw_data
raw_data<-raw_data_filt[raw_data_filt$count>1,]
# dim(raw_data)
edge_data<-raw_data[,c("id","RBP","count")]
edge_data$type<-"RBP binding"
# dim(raw_data)

## create nodes data
node_circRNA<-raw_data[,c("id","gene")]
node_circRNA<-node_circRNA[!duplicated(node_circRNA$id),]
node_circRNA$type<-"circRNA"
node_rbp<-data.frame("id"=raw_data[,"RBP"],"gene"=raw_data[,"RBP"])
node_rbp<-node_rbp[!duplicated(node_rbp),]
node_rbp$type<-"RBP"
# colnames(node_rbp)<-c("gene")
nodes<-rbind(node_circRNA,node_rbp)
net <- graph_from_data_frame(d=edge_data, vertices=nodes, directed=F)

## Set the attr to color
# Generate colors based on media type:
colrs <- c("#80b1d3", "#fb8072")
names(colrs)<-unique(V(net)$type)
V(net)$color<-NA
V(net)$color <- colrs[unlist(lapply(V(net)$type,function(x){
   which(names(colrs)==x)}))]
# Compute node degrees (#links) and use that to set node size:
deg <- degree(net, mode="all")
summary(deg)
# V(net)$size <- log10(deg+1)*10
# We could also use the audience size value:
# V(net)$size <- V(net)$audience.size*0.6
# The labels are currently node IDs.
# Setting them to NA will render no labels:
# node_type<-V(net)$type
# node_type[node_type=="circRNA"]<-NA
# V(net)$label<-V(net)
V(net)$label[V(net)$type=="RBP"]<-names(V(net)[V(net)$type=="RBP"])
V(net)$label[V(net)$type=="circRNA"] <- NA
V(net)$group<-"RBP"
V(net)$group[V(net)$type=="circRNA"&deg<=2]<-"lower degree"
V(net)$group[V(net)$type=="circRNA"&deg>2]<-"higher degree"
unique(V(net)$group)
# Set edge width based on weight:
E(net)$width <- 0.01
#change arrow size and edge color:
E(net)$arrow.size <- .02
E(net)$edge.color <- "gray80"
# # We can even set the network layout:
# graph_attr(net, "layout") <- layout_with_lgl
# plot(net)
# legend(x=-1.5, y=-0.5, c("circRNA","RBP"), pch=21,
       # col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)


edge.end <- ends(net, es=E(net), names=F)[,2]
edge.col <- V(net)$color[edge.end]
# plot(net, edge.color=edge.col, edge.curved=.1)

## whole layout
V(net)$label.cex=log10(deg+1)*0.5
V(net)$label.color="black"
V(net)$frame.color <- "white"
E(net)$arrow.mode <- 0
# plot(net,edge.curved=.1)
# legend(x=-1.5, y=-0.5, c("circRNA","RBP"), pch=21,
       # col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

# ## layout_nicely
# l<-layout_nicely(net,weights=log2(E(net)$count*1))
# plot(net,edge.curved=.1,layout=l,arrow.mode=0,edge.width=.01)
# ##
# tkid <- tkplot(net) #tkid is the id of the tkplot that will open
# l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
# plot(net, layout=l)

# ## Fruchterman-Reingold layout
# l <- layout_with_fr(net,niter=50,weights=log2(E(net)$count*1))
# plot(net,layout=l,margin=c(0,0,0,0),arrow.mode=0,edge.width=.01
     # )


# ## layout_with_graphopt
# l <- layout_with_graphopt(net,charge=0.02)
# plot(net, layout=l,edge.color=edge.col,)

##
library(qgraph)
e <- get.edgelist(net)
l <- qgraph.layout.fruchtermanreingold(seq(1,length(e[,1])),vcount=vcount(net),
         area=8*(vcount(net)^2),repulse.rad=(vcount(net)^3.1))
# plot(net, layout=l)

# Community detection (by optimizing modularity over partitions):
clp <- cluster_optimal(net)
clp_1<-cluster_edge_betweenness(net,directed = F)
clp<-cluster_infomap(net,modularity = T,nb.trials = 10)
class(clp)
# dendPlot(clp, mode="hclust")
modularity(clp) 
compare(clp_1, clp, method="rand")
compare(membership(sg), membership(le))
# Community detection returns an object of class "communities"
# which igraph knows how to plot:
# plot(clp, net,layout=l)

# We can also plot the communities without relying on their built-in plot:
V(net)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
# plot(net, vertex.color=colrs[V(net)$community])
##
# library(ggraph)
V(net)$community<-unlist(lapply(membership(clp),as.character))
# install.packages('pillar')V(net)$group
#####################################
## Figure 7A-degree/Figure S20
#####################################
ggraph(net,layout = l) +
   geom_edge_fan(color="gray50",alpha=0.5) +
   geom_node_point(aes(color=type,size=size)) +
   geom_node_text(label=V(net)$label,size=V(net)$size/4,color="gray10", repel=T)+
   scale_color_manual(values = c("#80b1d3","#fb8072","#FFA500"))+
   theme_void()+theme(
      legend.title = element_text(size=15),
      legend.text = element_text(size=12),
   )+
   guides(color = guide_legend(override.aes = list(size = 5)))

#####################################
## Figure 7A-community/Figure S21
#####################################
library(RColorBrewer)
community_colobar<-colorRampPalette(brewer.pal(9, "GnBu"))(100)
ggraph(net,layout = l) +
   geom_edge_fan(color="gray50",alpha=0.5) +
   geom_node_point(aes(color=community,size=size)) +
   geom_node_text(label=V(net)$label,size=V(net)$size/4,color="gray10", repel=T)+
   scale_color_gradientn(colors = community_colobar)+
   theme_void()+theme(
      legend.title = element_text(size=15),
      legend.text = element_text(size=12),
   )+
   guides(color = guide_legend(override.aes = list(size = 5)))
## RBP GO and KEGG enrichment
library(clusterProfiler)
library(org.Bt.eg.db)
#####################################
## Figure 7C
#####################################
##  GO and KEGG enrichment of different degrees or community circRNA host genes
degree_circRNA_list<-names(V(net)[V(net)$group=="higher degree"])
community_circRNA_list<-names(V(net)[V(net)$community>22&V(net)$type=="circRNA"])

# raw_data[raw_data$count>1]
hostgene_list<-unlist(lapply(degree_circRNA_list,function(x){strsplit(x,":")[[1]][1]}))
hostgene_list<-unlist(lapply(community_circRNA_list,function(x){strsplit(x,":")[[1]][1]}))

hostgene_list
# rbp_freq_golist<-rbp_freq$word[rbp_freq$freq>=7000]
# enrichGO(rbp_freq_golist,OrgDb = "bta")
eg_down = bitr(hostgene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Bt.eg.db")
id=eg_down$ENTREZID

ego_BP <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis

BP_res <- subset(ego_BP@result,select=c('Description','geneID','qvalue'))
BP <- BP_res[1:10,]
# BP_dot=dotplot(ego_BP,showCategory=10,title="Enrichment GO Top10") #dotplot
# BP_bar=barplot(ego_BP,showCategory=10,title="Enrichment GO Top10") #barplot
# BP_dot
# BP_bar
ego_MF <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "MF",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis
# MF_dot=dotplot(ego_MF,showCategory=10,title="Enrichment GO Top10") #dotplot
# MF_bar=barplot(ego_MF,showCategory=10,title="Enrichment GO Top10") #barplot
MF_res <- subset(ego_MF@result,select=c('Description','geneID','qvalue'))
MF <- MF_res[1:10,]
ego_CC <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "CC",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis
# CC=dotplot(ego_CC,showCategory=10,title="Enrichment GO Top10") #dotplot
# CC_bar=barplot(ego_CC,showCategory=10,title="Enrichment GO Top10") #barplot
CC_res <- subset(ego_CC@result,select=c('Description','geneID','qvalue'))
CC <- CC_res[1:10,]

#
GO_result_total <- rbind(BP_res,MF_res,CC_res)
GO_result_total$GO_term<-rep(c('BP','MF','CC'),c(length(BP_res$Description),length(MF_res$Description),length(CC_res$Description)))
# write.csv(GO_result_total,paste(output,'/','cis-genes for lncRNA',' GO enrichment result.csv',sep=''))
# GO enrichment result of Top10
GO_result<-rbind(BP,MF,CC)

GO_result$GO_term<-rep(c('BP','MF','CC'),rep(10,3))
GO_result$GO_term_1<-rep(c('BP','MF','CC'),rep(10,3))
GO_result<-GO_result[!is.na(GO_result$qvalue),]
library(stringr)
GO_result$GeneNumber<-str_count(GO_result$geneID,'/')
labels<-GO_result[GO_result$`-log10(qvalue)`>=-log10(0.05),]
GO_result$'-log10(qvalue)'= (-1)*log10(GO_result$qvalue)
png(paste(output, '/','cis-genes for lncRNA'," Go enrichment analysis-dotplot.png", sep=""),width=4000,height=3000,res=300)
p<-ggplot(GO_result,aes(x=GO_term,y=Description,color=-log10(qvalue),size=GeneNumber))+
   geom_point()+theme_bw()+
   theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),
         axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),
         legend.title = element_text(size = 22),
         legend.text = element_text(size = 18))+
   guides(color = guide_legend(override.aes = list(size = 5)))+
   scale_color_gradient(limits = c(min(GO_result$`-log10(qvalue)`), max(GO_result$`-log10(qvalue)`)), 
                        low = "green",
                        high = "red")+  
   annotate('text',x = labels$GeneRatio,y=labels$`-log10(qvalue)`,label=rownames(labels))

p

##############################################################################################################################
### miRNA binding site network
miRbind_data_filt<-read.table('./miRbind_data_filt.txt',header = T,sep='\t')
head(miRbind_data_filt)
raw_data<-miRbind_data_filt[,c(1:2,5:7)]
rm(miRbind_data_filt)
edge_data<-raw_data[raw_data$Site_type=="8mer-1a",]
mirna_node<-data.frame("id"=unique(edge_data[,"miRNA_name"]),"type"="miRNA")
circ_node<-data.frame("id"=unique(edge_data[,"a_Gene_ID"]),"type"="circRNA")
node_data<-rbind(mirna_node,circ_node)
dim(edge_data)
#####################################
## Figure 7F
#####################################
## wordcloud
rbp_freq<-data.frame("miRNA"=raw_data[,c("miRNA_name")])
rbp_freq$count<-1

rbp_freq<-data.frame(rbp_freq%>%group_by(miRNA)%>%summarize(count=n()))
colnames(rbp_freq)<-c('word','freq')
rownames(rbp_freq)<-rbp_freq$word
rbp_freq<-rbp_freq[order(rbp_freq$freq,decreasing = T),]
# pdf('./a.pdf')
wordcloud2(rbp_freq)
## miRNA target gene prediction
top10_miRNA<-rbp_freq$word[1:10]
top10_miRNA_family<-unlist(lapply(top10_miRNA,function(x){paste0(strsplit(x,'-',fixed = T)[[1]][2:3],collapse = '-')}))
# import prediction results from targetScan database
predict_res<-read.table('./prediction/Predicted_Targets_Info.default_predictions.txt',sep='\t',header = T)
head(predict_res)
predict_res$miRNA_family<-unlist(lapply(predict_res$miR.Family,function(x){len<-length(strsplit(x,'/')[[1]]);
   if (len == 1){
      return(x)
   } else {
      res<-paste0(c(strsplit(x,'/')[[1]][1],paste0("miR-",strsplit(x,'/')[[1]][2:len],collapse = ';')),collapse = ';')
      return(res)
   }
   }))
#
idx<-unlist(lapply(top10_miRNA_family, function(x){grep(x,predict_res$miRNA_family,ignore.case = T)}))
predict_res_filt<-predict_res[idx,]
predict_res_filt$PCT<-as.numeric(predict_res_filt$PCT,)
predict_res_final<-predict_res_filt[predict_res_filt$Species.ID==9913&!is.na(predict_res_filt$PCT)&
                                       predict_res_filt$Seed.match!="6-mer"&predict_res_filt$PCT>0.5,]
target_mRNA_symbol<-unique(predict_res_final$Gene.Symbol)
# Go enrichment analysis
target_mRNA_gene <- bitr(target_mRNA_symbol, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Bt.eg.db)
id<-target_mRNA_gene$ENTREZID
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

BP_res <- subset(ego_BP@result,select=c('Description','geneID','qvalue'))
BP <- BP_res[1:10,]

ego_MF <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "MF",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis

MF_res <- subset(ego_MF@result,select=c('Description','geneID','qvalue'))
MF <- MF_res[1:10,]
ego_CC <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "CC",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis

CC_res <- subset(ego_CC@result,select=c('Description','geneID','qvalue'))
CC <- CC_res[1:10,]


# write.csv(GO_result_total,paste(output,'/','cis-genes for lncRNA',' GO enrichment result.csv',sep=''))
# GO enrichment result of Top10
GO_result<-rbind(BP,MF,CC)

GO_result$GO_term<-rep(c('BP','MF','CC'),rep(10,3))
GO_result$GO_term_1<-rep(c('BP','MF','CC'),rep(10,3))
GO_result<-GO_result[!is.na(GO_result$qvalue),]
library(stringr)
GO_result$GeneNumber<-str_count(GO_result$geneID,'/')
labels<-GO_result[GO_result$`-log10(qvalue)`>=-log10(0.05),]
GO_result$'-log10(qvalue)'= (-1)*log10(GO_result$qvalue)
# png(paste(output, '/','cis-genes for lncRNA'," Go enrichment analysis-dotplot.png", sep=""),width=4000,height=3000,res=300)
#####################################
## Figure 7G
#####################################
ggplot(GO_result,aes(x=GO_term,y=Description,color=-log10(qvalue),size=GeneNumber))+
   geom_point()+theme_bw()+
   theme(axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
         axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),
         legend.title = element_text(size = 18),
         legend.text = element_text(size = 15))+
   guides(color = guide_legend(override.aes = list(size = 4)))+
   scale_color_gradient(limits = c(min(GO_result$`-log10(qvalue)`), max(GO_result$`-log10(qvalue)`)), 
                        low = "green", 
                        high = "red")+  
   annotate('text',x = labels$GeneRatio,y=labels$`-log10(qvalue)`,label=rownames(labels))


## construct the circRNA-miRNA network
net <- graph_from_data_frame(d=edge_data, vertices=node_data, directed=F)

# plot the network
# plot(net)
net <- simplify(net, remove.multiple = F, remove.loops = T)
# plot(net, edge.arrow.size=.04,vertex.size=5,vertex.label=NA,
#      edge.curved=.1)
## Set the attr to color
# Generate colors based on media type:
colrs <- c("#80b1d3", "#fb8072")
names(colrs)<-unique(V(net)$type)
V(net)$color<-NA
V(net)$color <- colrs[unlist(lapply(V(net)$type,function(x){
   which(names(colrs)==x)}))]
# Compute node degrees (#links) and use that to set node size:
deg <- degree(net, mode="all")
V(net)$size <- log10(deg+1)*10
# We could also use the audience size value:
# V(net)$size <- V(net)$audience.size*0.6
# The labels are currently node IDs.
# Setting them to NA will render no labels:
# node_type<-V(net)$type
# node_type[node_type=="circRNA"]<-NA
# V(net)$label<-V(net)
V(net)$label[V(net)$type=="miRNA"]<-names(V(net)[V(net)$type=="miRNA"])
V(net)$label[V(net)$type=="circRNA"] <- NA
V(net)$group<-"miRNA"
V(net)$group[V(net)$type=="circRNA"&deg<=2]<-"lower degree"
V(net)$group[V(net)$type=="circRNA"&deg>2]<-"higher degree"
# Set edge width based on weight:
E(net)$width <- 0.01
#change arrow size and edge color:
E(net)$arrow.size <- .02
E(net)$edge.color <- "gray80"
# We can even set the network layout:
graph_attr(net, "layout") <- layout_with_lgl

##
edge.end <- ends(net, es=E(net), names=F)[,2]
edge.col <- V(net)$color[edge.end]

## whole layout
V(net)$label.cex=log10(deg+1)*0.5
V(net)$label.color="black"
V(net)$frame.color <- "white"
E(net)$arrow.mode <- 0

##
library(qgraph)
install.packages('qgraph')
e <- get.edgelist(net)
l <- qgraph.layout.fruchtermanreingold(seq(1,length(e[,1])),vcount=vcount(net),
                                       area=8*(vcount(net)^2),repulse.rad=(vcount(net)^3.1))
# plot(net, layout=l)

# Community detection (by optimizing modularity over partitions):
clp <- cluster_optimal(net)
clp_1<-cluster_edge_betweenness(net,directed = F)
clp<-cluster_infomap(net,modularity = T,nb.trials = 10)
class(clp)
# dendPlot(clp, mode="hclust")
modularity(clp) 
compare(clp_1, clp, method="rand")
compare(membership(sg), membership(le))
# Community detection returns an object of class "communities"
# which igraph knows how to plot:
# plot(clp, net,layout=l)

##

# We can also plot the communities without relying on their built-in plot:
V(net)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net, vertex.color=colrs[V(net)$community])
##
library(ggraph)
V(net)$comm<-membership(clp)
# install.packages('pillar')V(net)$group
#####################################
## Figure 7D
#####################################
ggraph(net,layout = l) +
   geom_edge_fan(color="gray50",alpha=0.5) +
   geom_node_point(aes(color=type,size=size)) +
   geom_node_text(label=V(net)$label,size=V(net)$size/4,color="gray10", repel=T)+
   scale_color_manual(values = c("#80b1d3","#fb8072","#FFA500"))+
   theme_void()+theme(
      legend.title = element_text(size=15),
      legend.text = element_text(size=12),
   )+
   guides(color = guide_legend(override.aes = list(size = 5)))

#####################################
## Figure 7E
#####################################
library(RColorBrewer)
community_colobar<-colorRampPalette(brewer.pal(9, "GnBu"))(100)
ggraph(net,layout = l) +
   geom_edge_fan(color="gray50",alpha=0.5) +
   geom_node_point(aes(color=comm,size=size)) +
   geom_node_text(label=V(net)$label,size=V(net)$size/4,color="gray10", repel=T)+
   scale_color_gradientn(colors = community_colobar)+
   theme_void()+theme(
      legend.title = element_text(size=15),
      legend.text = element_text(size=12),
   )+
   guides(color = guide_legend(override.aes = list(size = 5)))


#####################################
## Figure S22
#####################################
##  GO and KEGG enrichment of different degrees or community circRNA host genes
degree_circRNA_list<-names(V(net)[V(net)$group=="lower degree"])
community_circRNA_list<-names(V(net)[V(net)$community>44&V(net)$type=="circRNA"])

# raw_data[raw_data$count>1]
hostgene_list<-unlist(lapply(degree_circRNA_list,function(x){strsplit(x,":")[[1]][1]}))
hostgene_list<-unlist(lapply(community_circRNA_list,function(x){strsplit(x,":")[[1]][1]}))

hostgene_list
# rbp_freq_golist<-rbp_freq$word[rbp_freq$freq>=7000]
# enrichGO(rbp_freq_golist,OrgDb = "bta")
eg_down = bitr(hostgene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Bt.eg.db")
id=eg_down$ENTREZID

ego_BP <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis

BP_res <- subset(ego_BP@result,select=c('Description','geneID','qvalue'))
BP <- BP_res[1:10,]
# BP_dot=dotplot(ego_BP,showCategory=10,title="Enrichment GO Top10") #dotplot
# BP_bar=barplot(ego_BP,showCategory=10,title="Enrichment GO Top10") #barplot
# BP_dot
# BP_bar
ego_MF <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "MF",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis
# MF_dot=dotplot(ego_MF,showCategory=10,title="Enrichment GO Top10") #dotplot
# MF_bar=barplot(ego_MF,showCategory=10,title="Enrichment GO Top10") #barplot
MF_res <- subset(ego_MF@result,select=c('Description','geneID','qvalue'))
MF <- MF_res[1:10,]
ego_CC <- enrichGO(
   gene  = id,
   keyType = "ENTREZID",
   OrgDb   = org.Bt.eg.db,
   ont     = "CC",
   pAdjustMethod = "BH",
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.05,
   readable      = TRUE) #GO enrichment analysis
# CC=dotplot(ego_CC,showCategory=10,title="Enrichment GO Top10") #dotplot
# CC_bar=barplot(ego_CC,showCategory=10,title="Enrichment GO Top10") #barplot
CC_res <- subset(ego_CC@result,select=c('Description','geneID','qvalue'))
CC <- CC_res[1:10,]

#
GO_result_total <- rbind(BP_res,MF_res,CC_res)
GO_result_total$GO_term<-rep(c('BP','MF','CC'),c(length(BP_res$Description),length(MF_res$Description),length(CC_res$Description)))
# write.csv(GO_result_total,paste(output,'/','cis-genes for lncRNA',' GO enrichment result.csv',sep=''))
# GO enrichment result of Top10
GO_result<-rbind(BP,MF,CC)

GO_result$GO_term<-rep(c('BP','MF','CC'),rep(10,3))
GO_result$GO_term_1<-rep(c('BP','MF','CC'),rep(10,3))
GO_result<-GO_result[!is.na(GO_result$qvalue),]
library(stringr)
GO_result$GeneNumber<-str_count(GO_result$geneID,'/')
labels<-GO_result[GO_result$`-log10(qvalue)`>=-log10(0.05),]
GO_result$'-log10(qvalue)'= (-1)*log10(GO_result$qvalue)
png(paste(output, '/','cis-genes for lncRNA'," Go enrichment analysis-dotplot.png", sep=""),width=4000,height=3000,res=300)
p<-ggplot(GO_result,aes(x=GO_term,y=Description,color=-log10(qvalue),size=GeneNumber))+
   geom_point()+theme_bw()+
   theme(axis.title.x = element_text(size=25),axis.title.y = element_text(size=25),
         axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),
         legend.title = element_text(size = 22),
         legend.text = element_text(size = 18))+
   guides(color = guide_legend(override.aes = list(size = 5)))+
   scale_color_gradient(limits = c(min(GO_result$`-log10(qvalue)`), max(GO_result$`-log10(qvalue)`)), # 数据上下限
                        low = "green", # 下限颜色
                        high = "red")+  # 中点值
   annotate('text',x = labels$GeneRatio,y=labels$`-log10(qvalue)`,label=rownames(labels))

p

