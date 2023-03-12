#!/bin/bash

##############################################################################
# Format genotype datasets (including BGVD database and two SNP-chip sets)
# Author: Hongfei Liu
# 
##############################################################################

### Global variables
WORKING_DIR=${YOUR_WORKING_DIRECTORY}
GENOTYPE_DATASETS_DIR=${WORKING_DIR}/genotype_datasets
BTAU6_GENOME=${WORKING_DIR}/genome/bosTau6.fa
OUTPUT_DIR=${WORKING_DIR}/output
SUB_OUTPUT_DIR=${OUTPUT_DIR}/imputation
TEMP_DIR=${OUTPUT_DIR}/imputation/temp
INPUT_DIR=${OUTPUT_DIR}/format_genotype_datasets

# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
java=${JAVA}
beagle=${path_to_beagle_jar}

### Basic config
THREADS=80

### Imputation by beagle
## GSE95358
mkdir -p ${SUB_OUTPUT_DIR}/GSE95358
$java -Xmx4096m -jar $beagle gt=${INPUT_DIR}/GSE95358_vcf_format.vcf \
nthreads=$THREADS ne=236 out=${SUB_OUTPUT_DIR}/GSE95358/GSE95358_imp
## GSE100038
mkdir -p ${SUB_OUTPUT_DIR}/GSE100038
$java -Xmx4096m -jar $beagle gt=${INPUT_DIR}/GSE100038_vcf_format.vcf \
nthreads==$THREADS ne=288 out=${SUB_OUTPUT_DIR}/GSE100038/GSE100038_imp

### Filtration by vcftools to discard SNPs with MAF (minor allele frequency) < 0.01 
vcftools --gzvcf ${SUB_OUTPUT_DIR}/GSE95358/GSE95358_imp.vcf.gz \
--maf 0.01 --hwe 0.05 --recode --recode-INFO-all --out ${SUB_OUTPUT_DIR}/GSE95358/GSE95358_filt
vcftools --gzvcf ${SUB_OUTPUT_DIR}/GSE100038/GSE100038_imp.vcf.gz \
--maf 0.01 --hwe 0.05 --recode --recode-INFO-all --out ${SUB_OUTPUT_DIR}/GSE100038/GSE100038_filt
