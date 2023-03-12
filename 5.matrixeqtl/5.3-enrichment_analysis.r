#!/bin/Rscript
library(ggplot2)

############################################################################################################
working_dir <- "path_to_your_working_dir"
base.dir='./output/circQTLs'
setwd(working_dir)

############################################################################################################
#####################################
## Figure S3-4
#####################################
#### Enrichment of functional regions for circQTLs
## MonteCarlo function
monteCarlo<-function(circqtl_res_cons,region_type,iteration){
   # ''' 
   # circqtl_res_cons: data.frame marked dataset
   # region_type: character functional region in genome
   # iteration: int  number of iteration for sampling
   # return: pvalue
   # '''
   # iteration=1000
   
   background_total_number<-length(circqtl_res_cons$snps[circqtl_res_cons$isCircQTL==F])
   circqtl_total_number<-length(circqtl_res_cons$isCircQTL[circqtl_res_cons$isCircQTL==T])
   enrichment_freq<-0
   odds_ratio<-0
   odds_ratio_list<-list()
   for (i in 1:iteration){
      # sampling without duplicates
      random_res<-sample(circqtl_res_cons$snps,size = 5000,replace = FALSE)
      # length(random_res)
      # length(unique(random_res))  
      region_type="intergenic_variant"
      circ_enrich_number<-length(intersect(random_res,
                                           circqtl_res_cons$snps[circqtl_res_cons$isCircQTL==T&
                                                                    circqtl_res_cons$ConsequenceType==region_type]))
      back_enrich_number<-length(intersect(random_res,
                                           circqtl_res_cons$snps[circqtl_res_cons$isCircQTL==F&
                                                                    circqtl_res_cons$ConsequenceType==region_type]))
      circqtl_freq<-circ_enrich_number/circqtl_total_number;back_freq<-back_enrich_number/background_total_number
      odds_ratio<-odds_ratio+circqtl_freq/(back_freq)
      odds_ratio_list<-append(odds_ratio_list,circqtl_freq/(back_freq))
      # if (back_freq=)
      if (back_freq>=circqtl_freq){
         enrichment_freq<-enrichment_freq+1
      }
      
   }
   print(enrichment_freq/iteration)
   return(list(c(enrichment_freq/iteration,odds_ratio/iteration),odds_ratio_list))
}
## import circQTL result
circqtl<-read.table('./output/circQTLs/significant_circQTL_log_quant_hfilt.txt',header = T)
circqtl_filt<-subset(circqtl,select=c(1,8))
length(circqtl_filt$VariantID)
colnames(circqtl_filt)[1]<-'snps'

## import high-confidence circQTL result
circqtl_res<-read.table('./output/circQTLs/hvaild_circQTLs.txt',header = T)
snp_info<-snp_info[,c("VariantID","ConsequenceType")]
colnames(snp_info)[1]<-"snps"
circqtl_res<-merge(circqtl_res,snp_info,by='snps')
circqtl_res$isCircQTL<-TRUE
circqtl_res_rbind<-circqtl_res[,c("snps","ConsequenceType","isCircQTL")]

## The enrichment of circQTLs in distribution of consequenceType (background: all SNPs(circQTL+eQTLs+others) in BGVD)
# all_snp_info<-read.table('../bovine_vcf (BGVD)/Btau_5.0.1_SNPs.anno.info',header = T,sep = '\t')
# snp_info<-read.table('../bovine_vcf (BGVD)/Btau_5.0.1_SNPs.anno.info',header = T,sep = '\t')

## The enrichment of circQTLs in distribution of consequenceType (background: all SNPs in this study)
all_snps<-snp_info[,c("snps","ConsequenceType")];colnames(all_snps)<-c("snps","ConsequenceType")
colnames(snp_info)
idx<-unlist(parLapply(cl,unique(circqtl_res$snps),function(x,all_snps){which(all_snps$snps==x)},all_snps))
# rm(snp_info);gc()
all_snps_filt<-all_snps
all_snps_filt$isCircQTL<-FALSE
all_snps_filt$isCircQTL[idx]<-TRUE


## The enrichment of circQTLs in distribution of consequenceType (background: circQTL+eQTLs in this study)
all_snps_filt<-eqtl_filt_res[,c("snps","ConsequenceType")]

# idx<-unique(round(runif(1e5,min=1,max=dim(all_snps)[1])))
# all_snps_filt<-all_snps[idx,]

# idx<-unlist(parLapply(cl,unique(circqtl_res$snps),function(x,all_snps){which(all_snps$snps==x)},all_snps))
# length(idx)
# get the intersection of eQTL and circQTL
intersects<-intersect(all_snps_filt$snps,circqtl_res_rbind$snps)

all_snps_filt$isCircQTL<-FALSE
all_snps_filt<-all_snps_filt[!duplicated(all_snps_filt$snps),]

all_snps_filt<-rbind(all_snps_filt,circqtl_res_rbind[!duplicated(circqtl_res_rbind$snps),])
all_snps_filt<-all_snps_filt[!duplicated(all_snps_filt$snps),]
idx<-unlist(lapply(intersects,function(x){which(all_snps_filt$snps==x)}))
all_snps_filt$isCircQTL[idx]<-F

sum(all_snps_filt$isCircQTL==T)
# dim(circqtl_res_rbind[!duplicated(circqtl_res_rbind$snps),])
# length(circqtl_res_rbind$isCircQTL)
## Fisher's exact Test or MonteCarlo method
circqtl_res_cons<-all_snps_filt
type<-unique(circqtl_res_cons$ConsequenceType)
pvalue<-list();OR<-list()
# construct the 2X2 table for fisher's exact test
for (i in 1:length(type)) {
   #beta>0 & isin the type
   a<-nrow(circqtl_res_cons[circqtl_res_cons$isCircQTL==TRUE&circqtl_res_cons$ConsequenceType==type[i],])
   b<-nrow(circqtl_res_cons[circqtl_res_cons$isCircQTL==TRUE&circqtl_res_cons$ConsequenceType!=type[i],])
   c<-nrow(circqtl_res_cons[circqtl_res_cons$isCircQTL==FALSE&circqtl_res_cons$ConsequenceType==type[i],])
   d<-nrow(circqtl_res_cons[circqtl_res_cons$isCircQTL==FALSE&circqtl_res_cons$ConsequenceType!=type[i],])
   table<-data.frame('Is'=c(a,c),'Not'=c(b,d))
   print(table)
   res<-fisher.test(table,alternative = 'greater')
   pvalue[i]<-res$p.value
   OR[i]<-res$estimate
   
}
# MonteCarlo method
pvalue<-list();OR<-list()
monteCarloDf<-data.frame('oddsRatio'='','pvalue'='')
for (i in 1:length(type)) {
   # type[15]
   res<-monteCarlo(circqtl_res_cons,type[i],1e4)
   
   pvalue[i]<-res[[1]][1]
   OR[i]<-res[[1]][2]
   # print(unlist(res[[2]]))
   odds_ratio_list<-res[[2]]
   print(length(unlist(odds_ratio_list)))
   idx<-(i-1)*1e4+1
   monteCarloDf[idx:(idx+1e4-1),"oddsRatio"]<-unlist(odds_ratio_list)
   monteCarloDf[idx:(idx+1e4-1),"pvalue"]<-res[[1]][1]
   monteCarloDf[idx:(idx+1e4-1),"type"]<-type[i]
}
## 
library(plyr)
library(Rmisc)

# The following auxiliary function will be used to calculate the mean and standard deviation
data_summary <- function(data, varname, groupnames){
   require(plyr)
   summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
   }
   data_sum<-ddply(data, groupnames, .fun=summary_func,
                   varname)
   data_sum <- rename(data_sum, c("mean" = varname))
   return(data_sum)
}
monteCarloDf$oddsRatio<-as.numeric(monteCarloDf$oddsRatio)
monteCarloDf<-monteCarloDf[monteCarloDf$oddsRatio!=Inf,]
monteCarloDf_temp<-monteCarloDf[!duplicated(monteCarloDf$type),]

# plot_df<-data_summary(monteCarloDf,varname = "oddsRatio",groupnames = "type")
plot_df <- summarySE(monteCarloDf, measurevar = "oddsRatio",
                     groupvars = c("type"))
idx<-unlist(lapply(plot_df$type,function(x){which(monteCarloDf_temp$type==x)}))
plot_df$pvalue<-as.numeric(monteCarloDf_temp$pvalue[idx])
plot_df<-plot_df[!is.nan(plot_df$oddsRatio)&plot_df$oddsRatio!=0&!is.na(plot_df$type),]
plot_df$CI95<-plot_df$se*1.96

ggplot(plot_df)+geom_point(aes(x=oddsRatio,y=type))+
   geom_errorbar(aes(xmin = oddsRatio - CI95, xmax = oddsRatio + CI95, y = type),size=1, 
                 width = 0.2, position = position_dodge(0.9))+
   geom_vline(xintercept = 1,color='black',linetype='dashed',size=1)+
   xlab('Odds ratio (compared to non-circQTL)')+
   ggrepel::geom_text_repel(data = plot_df,
                            aes(x=oddsRatio,y=type,label=c(paste0("~italic(P) ==",round(pvalue,3)))),size=5,,parse = T)+
   theme_bw()+theme(axis.title.x = element_text(size=25),
                    axis.title.y = element_blank(),
                    axis.text = element_text(size=18),
                    legend.title = element_text(size=15),
                    legend.text = element_text(size=15))
## bubble plot for Fisher's test
# combine
enrich_res<-data.frame('Type'=type,'pvalue'=unlist(pvalue),'OR'=unlist(OR))

highlight<-enrich_res[(enrich_res$pvalue<=10**(-0.25))&(enrich_res$OR>=1),]
rand<-runif(length(highlight[,1]),0,0.2)
rand1<-runif(length(highlight[,1]),-0.1,0.1)
highlight$ORend<-highlight$OR-rand;highlight$pvalueend<-highlight$pvalue+rand1
n=length(enrich_res[,1])-length(highlight[,1])
blank<-data.frame('Type'=rep('',n),'pvalue'=rep(0,n),'OR'=rep(0,n),'pvalueend'=0,'ORend'=0)
highlight<-rbind(highlight,blank)
highlight$ORend[highlight$Type=='non_coding_transcript_exon_variant']=1.28
# bubble plot for Fisher's test
enrich_res<-enrich_res[enrich_res$OR!=Inf,]
ggplot(enrich_res)+geom_point(aes(x=OR,y=-log10(pvalue)),color='lightblue',size=10,alpha=0.8)+
   xlab('Odds ratio #(circQTL)/#(non-circQTL)')+
   ylab(expression(paste(-log[10],'(',italic("P")," value)",sep='')))+
   geom_vline(xintercept = 1,color='red',linetype='dashed')+
   geom_hline(yintercept = -log10(0.05),color='red',linetype='dashed')+
   # ylim(0,1000)+
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


### The enrichment of effects of circQTLs in distribution of consequenceType
circqtl_res<-read.table('./output/circQTLs/hvaild_circQTLs.txt',header = T)
circqtl_res_cons<-merge(circqtl_res,circqtl_filt,by='snps')
# circqtl_res_cons<-circqtl_res_cons[circqtl_res_cons$ConsequenceType]
type<-circqtl_res_cons$ConsequenceType[!duplicated(circqtl_res_cons$ConsequenceType)]
pvalue<-list();OR<-list()
for (i in 1:length(type)) {
   #beta>0 & isin the type
   a<-nrow(circqtl_res_cons[circqtl_res_cons$beta>0&circqtl_res_cons$ConsequenceType==type[i],])
   b<-nrow(circqtl_res_cons[circqtl_res_cons$beta>0&circqtl_res_cons$ConsequenceType!=type[i],])
   c<-nrow(circqtl_res_cons[circqtl_res_cons$beta<0&circqtl_res_cons$ConsequenceType==type[i],])
   d<-nrow(circqtl_res_cons[circqtl_res_cons$beta<0&circqtl_res_cons$ConsequenceType!=type[i],])
   table<-data.frame('Is'=c(a,c),'Not'=c(b,d))
   print(table)
   res<-fisher.test(table,alternative = 'greater')
   pvalue[i]<-res$p.value
   OR[i]<-res$estimate
   
}
##
circqtl_res_cons$isCircQTL<-FALSE
circqtl_res_cons$isCircQTL[circqtl_res_cons$beta>0]<-T
## MonteCarlo method
monteCarloDf<-data.frame('oddsRatio'='','pvalue'='')
for (i in 1:length(type)) {
   # type[15]
   res<-monteCarlo(circqtl_res_cons,type[i],1e4)
   
   pvalue[i]<-res[[1]][1]
   OR[i]<-res[[1]][2]
   # print(unlist(res[[2]]))
   odds_ratio_list<-res[[2]]
   print(length(unlist(odds_ratio_list)))
   idx<-(i-1)*1e4+1
   monteCarloDf[idx:(idx+1e4-1),"oddsRatio"]<-unlist(odds_ratio_list)
   monteCarloDf[idx:(idx+1e4-1),"pvalue"]<-res[[1]][1]
   monteCarloDf[idx:(idx+1e4-1),"type"]<-type[i]
}

# combine
enrich_res<-data.frame('Type'=type,'pvalue'=unlist(pvalue),'OR'=unlist(OR))

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
   xlab('Odds ratio #(Gamma>0)/#(Gamma)')+
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

##
monteCarloDf$oddsRatio<-as.numeric(monteCarloDf$oddsRatio)
monteCarloDf<-monteCarloDf[monteCarloDf$oddsRatio!=Inf,]
monteCarloDf_temp<-monteCarloDf[!duplicated(monteCarloDf$type),]

# plot_df<-data_summary(monteCarloDf,varname = "oddsRatio",groupnames = "type")
plot_df <- summarySE(monteCarloDf, measurevar = "oddsRatio",
                     groupvars = c("type"))
idx<-unlist(lapply(plot_df$type,function(x){which(monteCarloDf_temp$type==x)}))
plot_df$pvalue<-as.numeric(monteCarloDf_temp$pvalue[idx])
plot_df<-plot_df[!is.nan(plot_df$oddsRatio)&plot_df$oddsRatio!=0&!is.na(plot_df$type),]
plot_df$CI95<-plot_df$se*1.96

ggplot(plot_df)+geom_point(aes(x=oddsRatio,y=type))+
   geom_errorbar(aes(xmin = oddsRatio - CI95, xmax = oddsRatio + CI95, y = type),size=1, 
                 width = 0.2, position = position_dodge(0.9))+
   # geom_vline(xintercept = 1,color='black',linetype='dashed',size=1)+
   xlab('Odds ratio (compared to Gamma < 0)')+
   ggrepel::geom_text_repel(data = plot_df,
                            aes(x=oddsRatio,y=type,label=c(paste0("~italic(P) ==",round(pvalue,3)))),size=5,,parse = T)+
   theme_bw()+theme(axis.title.x = element_text(size=25),
                    axis.title.y = element_blank(),
                    axis.text = element_text(size=18),
                    legend.title = element_text(size=15),
                    legend.text = element_text(size=15))
