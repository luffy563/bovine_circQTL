#!/bin/Rscript

library(MatrixEQTL)
library(CMplot)

############################################################################################################
working_dir <- "path_to_your_working_dir"
base.dir='./output/eQTL'
setwd(working_dir)

##Set the analysis model
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 
##Set the SNP file
SNP_file_name = paste(base.dir, "/snp.txt", sep="");
snps_location_file_name = paste(base.dir, "/snp_pos.txt", sep="");

## adjust the expression data grouped by breeds
mRNA_expression_data <- read.table('./combined_mrna.txt',header = T,sep='\t')
header = colnames(mRNA_expression_data)
# read the breed information
breed_info = read.csv('./SRR_breeds.csv')
# breed_info<-breed_info[unlist(lapply(breed_info$SRR_list,function(x){x%in%all_gtf})),]
breeds = breed_info$Breeds[!duplicated(breed_info$Breeds)]
mRNA_expression_data_convert=data.frame(mRNA_expression_data[,header[1:6]])
for (i in 1:length(breeds)) {
  srr = breed_info$SRR_list[breed_info$Breeds==breeds[i]]
  mRNA_expression_data_breed = subset(mRNA_expression_data,select = srr)
  mRNA_expression_data_convert[,as.character(breeds[i])]=apply(mRNA_expression_data_breed,1,mean)
}
# output the related annotation info to file
write.table(mRNA_expression_data_convert[,c("Gene.ID","Reference","Start","End")],'./output/matrixeqtl/mRNA_info.txt',
            sep = '\t',row.names = F,quote = F)
# output the expression_data_convert to file
write.table(mRNA_expression_data_convert[,c("Gene.ID",breeds)],'./output/matrixeqtl/mRNA.txt',
            sep = '\t',row.names = F,quote = F)

# Set the expression file of mRNA
expression_file_name = paste(base.dir, "/mRNA.txt", sep="")
gene_location_file_name = paste(base.dir, "/mRNA_info.txt", sep="")

# set the location of output file
output_file_name_cis = './output/matrixeqtl/cis-eqtl.txt'
output_file_name_tra = './output/matrixeqtl/tra-eqtl.txt'

# set the covariate
covariates_file_name = paste(base.dir, "./output/cov.txt", sep="")

# pvOutputThreshold = 1e-4
pvOutputThreshold_tra = 1e-4
pvOutputThreshold_cis = 1e-4

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
  noFDRsaveMemory = FALSE)

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values
plot(me)

## 
trans_eqtl<-read.table('./output/matrixeqtl/tra-eqtl.txt',header = T)
cis_eqtl<-read.table('./output/matrixeqtl/cis-eqtl.txt',header = T)
eqtl<-rbind(trans_eqtl,cis_eqtl)
eqtl_filt<-eqtl[eqtl$FDR<=0.05,]
# eqtl_filt<-eqtl
eqtl_filt<-eqtl_filt[!duplicated(eqtl_filt$SNP),]
head(eqtl_filt)
colnames(eqtl_filt)[1]<-'snps'
# dim(eqtl_filt)

# import all_test_snps
data=read.table('./output/combined/combined.vcf',sep=' ',header = T)
snp_info=data[,1:15]
colnames(snp_info)[8]<-"snps"
eqtl_filt_res<-merge(eqtl_filt,snp_info,by='snps')
write.table(eqtl_filt_res,'./output/matrixeqtl/eqtl_filt_res.txt',sep='\t',quote = F,row.names = F)
