
#!/usr/lib/Rscript
# input files names
library('optparse')
library(stringr)

option_list = list(
  make_option(c("-s", "--single"), type = "integer", default = FALSE,
              action = "store", help = "The input bed file is output based on single-end seq or paired-end seq!"),
  make_option(c("-r", "--remap"), type = "integer", default = FALSE,
              action = "store", help = "The output file is the remapped and combined bed file! (only used for URBORUS)"),
  make_option(c("-q", "--quant"), type = "integer", default = FALSE,
              action = "store", help = "The output bed file is used for the quant by CIRIquant"),
  
  make_option(c("-a", "--algorithm"), type = "character", default = FALSE,
              action = "store", help = "The input bed file is based on which algorithm.
		  FC: find_circ
		  CIRC2: CIRCexplorer2
		  CIRI2: CIRI2
		  CDBG: CircDBG
		  CM: CircMarker
		  MP: Mapsplice
		  KF: KNIFE
		  NCL: NCLscan
		  UB: URBORUS
		  DCC:DCC
		  SE:Segemehl"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# print(opt)

# which circRNAs prediction software
if (opt$a=='FC') {
  path='output/find_circ/'} else if (opt$a=='CIRC2') {
    path='output/tophat_fusion/'
  } else if (opt$a=='CDBG') {
    path='output/CircDBG/'
  } else if (opt$a=='CM') {
    path='output/CircMarker/'
  } else if (opt$a=='MP') {
    path='output/Mapsplice/'
  } else if (opt$a=='KF') {
    path='output/KNIFE/'
  } else if (opt$a=='CF') {
    path='output/CircRNAfinder/'
  } else if (opt$a=='MP') {
    path='output/Mapsplice/'
  } else if (opt$a=='CIRI2') {
    path='output/CIRI2/'
  } else if (opt$a=='NCL') {
    path='output/NCLscan/'
  } else if (opt$a=='UB') {
    path='output/UROBORUS/'
  } else if (opt$a=='DCC') {
    path='output/DCC/'
  } else if (opt$a=='SE') {
    path='output/segemehl/'
  }
# single-end seq or paired-end seq
if (opt$s==1) {
  filenames=read.table('SRR_list.txt')}  else {
    filenames=read.table('sra_list.txt')
  }

for (i in 1:length(filenames[,1])) {
  a=as.character(filenames[,1][i])
  print(a)
  if (opt$a=='FC') {
    ipath=paste(path,a,'/circ_candidates.bed',sep='')} else if (opt$a=='CIRC2') {
    ipath=paste(path,a,'/circ_candidates.bed',sep='')} else if (opt$a=='CF'){
    ipath=paste(path,a,'/',a,'s_filteredJunctions.bed',sep='')} else if (opt$a=='MP'){
    ipath=paste(path,a,'/fusions_raw.txt',sep='')} else if (opt$a=='CDBG'){
    ipath=paste(path,a,'/Detection_Result/Brief_sum.txt',sep='')} else if (opt$a=='CM'){
    ipath=paste(path,a,'/Detection_Result/Brief_sum.txt',sep='')} else if (opt$a=='CIRI2'){
    ipath=paste(path,a,'/out.ciri',sep='')} else if (opt$a=='NCL'){
    ipath=paste(path,a,'/NCLscan_out.result.info',sep='')} else if (opt$a=='UB'){
    ipath=paste(path,a,'/circRNA_list.txt',sep='')} else if (opt$a=='DCC'){
    ipath=paste(path,a,'/Alignment/',a,'Chimeric.out.junction.circRNA',sep='')} else if (opt$a=='SE'){
    ipath=paste(path,a,'/',a,'.sum.bed',sep='')}
  print(ipath)


  if (opt$a=='CIRI2') {
    data=read.table(ipath, comment.char='') 
    data=data[-1,] } else if (opt$a=='CM') { 
      data=read.table(ipath,sep='')} else if (opt$a=='CDBG') { 
        data=read.table(ipath,sep='')} else {
          data=read.table(ipath,sep='\t')}
V1=data$V1
if (opt$a!='CIRI2') {
  if (opt$a=='CM') {
  V3=data$V2
  data$V2=data$V3
  data$V3=V3
  }
V4=paste(data$V1,':',data$V2,'|',data$V3,sep = '')
data$V4=V4
}
V5=rep('.',length(V1))
if (opt$a=='CM') {
  V6=data$V5
  data$V5=V5
  data$V6=V6
} else if (opt$a=='CIRI2') {
  data$V5=data$V1
  data$V6=V5
  data$V7=data$V11
  } else {
  data$V5=V5
  }


if (opt$a=='CIRI2') {
  data = data[,1:7]
  data=data[,-1]
} else {
  data = data[,1:6]
}
write.table(data,paste(path, a, '/circ_candidates_CIRIquant.bed',sep=''),row.names=FALSE, col.names=FALSE,sep='\t',quote=F)
}
