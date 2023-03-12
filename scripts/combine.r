# Combine all genotype_table files
rm(list = ls()) #绉婚櫎鎵€鏈夊彉閲?
if (!require(readr) == TRUE) {
  install.packages("readr")
}
if (!require(parallel) == TRUE) {
  install.packages("parallel")
}
if (!require(optparse) == TRUE) {
  install.packages("optparse")
}
library(optparse)
library(readr)
library(parallel)

option_list = list(
  make_option(c("-p", "--pthreads"), type = "integer", default = 1,
              action = "store", help = "The thread number of running this script"),
  make_option(c("-t", "--tholstein"), type = "character", default = FALSE,
              action = "store", help = "The input freq table of Holstein cattle population"),
  make_option(c("-q", "--qinchuan"), type = "character", default = FALSE,
              action = "store", help = "The input freq table of Qinchuan cattle population"),
  make_option(c("-b", "--bgvd"), type = "character", default = FALSE,
              action = "store", help = "The input tab file from BGVD database"),
  make_option(c("-o", "--out"), type = "character", default = FALSE,
              action = "store", help = "Output file path"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (opt$p == F) {
  print('Please input the options')
}
# Making cluster notes
cl <- makeCluster(3)

# GSE95358: Holstein
hols_path <- opt$t
qinc_path <- opt$q
bgvd_data <- opt$b
output_path <- opt$o

# hols_path <- './output/freq_table/GSE95358.frq'
# qinc_path <- './output/freq_table/GSE100038.frq'
# bgvd_data <- 'E:/Btau_5.0.1_SNPs.anno.tab'
# output_path <- './combine.vcf'

# reading datasets
hols <- read.table(hols_path,row.names=NULL,header = T, sep='\t')
colnames(hols)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'Ref', 'Alt')
n_hols <- hols$N_CHR[1]/2
hols_loc <- paste(hols$CHROM,':',hols$POS,sep='')
# GSE100038: Qinchuan
qinc <- read.table(qinc_path,row.names=NULL,header = T, sep='\t')
colnames(qinc)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'Ref', 'Alt')
qinc_loc <- paste(qinc$CHROM,':',qinc$POS,sep='')
n_qinc <- qinc$N_CHR[1]/2
inter <- intersect(qinc_loc,hols_loc)


# 鐗瑰畾閮ㄥ垎鐨勬暟鎹?
ptm <- proc.time() #鐢ㄤ簬璁＄畻璇诲彇鑰楁椂
con <- file(bgvd_data, "r")
line=readLines(con,n=1)
colname <- strsplit(line,'\t',fixed = T)
combined_line <- colname
init <- matrix(0,1,length(colname[[1]]))
data <- as.data.frame(init)
colnames(data) <- unlist(colname)
temp <- data
data$Qinchuan<-0

write.table(data[0,],output_path,
            quote = F ,sep='\t')
n<-100000
ptm <- proc.time() #鐢ㄤ簬璁＄畻璇诲彇鑰楁椂
i <- 1
while( length(line) != 0 ) {
  line=readLines(con,n=100000)
  term <- strsplit(line,'\t',fixed = T)

  temp[1:length(term),] <- data.frame(t(data.frame(term)))
  # 鏇存柊Holstein閮ㄥ垎鐨勬暟鎹?
  colnames(temp) <- unlist(colname)
  bgvd_loc <- temp$UMD3_1_1_Pos
  idx <- parLapply(cl,unlist(inter),function(c,b)which(b==c),bgvd_loc)
  hols_idx <- parLapply(cl,bgvd_loc[unlist(idx)],function(c,b)which(b==c),hols_loc)
  qinc_idx <- parLapply(cl,bgvd_loc[unlist(idx)],function(c,b)which(b==c),qinc_loc)
  bgvd_hols_f <- data.frame(t(as.data.frame(strsplit(temp$Holstein[unlist(idx)],':',fixed = T))))
  colnames(bgvd_hols_f) <- c('Number','Freq')
  n_bgvd <-  as.numeric(bgvd_hols_f$Number[1])
  hols_f <- data.frame(t(as.data.frame(strsplit(hols$Alt[unlist(hols_idx)],':',fixed = T))))
  colnames(hols_f) <- c('Alt','Freq')
  qinc_f <- data.frame(t(as.data.frame(strsplit(qinc$Alt[unlist(qinc_idx)],':',fixed = T))))
  colnames(qinc_f) <- c('Alt','Freq')
  up_hols_f <- (n_bgvd *  as.numeric(bgvd_hols_f$Freq) + n_hols *  as.numeric(hols_f$Freq))/(n_bgvd + n_hols)
  up_n_hols <- (n_bgvd + n_hols)
  temp$Holstein[unlist(idx)] <- paste(up_n_hols,':',up_hols_f,sep='')
  temp$Qinchuan<-0
  temp$Qinchuan[unlist(idx)] <- paste(n_qinc,':',qinc_f$Freq,sep='')
  sprintf('Complete iteration: %d',i * 100000)
  write.table(temp[unlist(idx),],output_path,
              quote=F,row.names=F,col.names=F,append = T)
  i <- i+1
  temp <- data.frame(matrix(0,1,length(colname[[1]])))

}
t <- proc.time() - ptm
close(con)
stopCluster(cl)
sprintf('Complete costs: %f min',t[3]/60)

