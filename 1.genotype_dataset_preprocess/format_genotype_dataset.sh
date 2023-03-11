#!/bin/bash

##############################################################################################
# Format genotype datasets (including BGVD database and two SNP-chip sets) into vcf format
# Author: Hongfei Liu
# 
##############################################################################################

### Global variables
WORKING_DIR=${YOUR_WORKING_DIRECTORY}
GENOTYPE_DATASETS_DIR=${WORKING_DIR}/genotype_datasets
BTAU6_GENOME=${WORKING_DIR}/genome/bosTau6.fa
OUTPUT_DIR=${WORKING_DIR}/output
SUB_OUTPUT_DIR=${OUTPUT_DIR}/format_genotype_datasets
TEMP_DIR=${OUTPUT_DIR}/format_genotype_datasets/temp
# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/format_genotype_datasets
mkdir -p ${OUTPUT_DIR}/format_genotype_datasets/temp

### Required software
iaap=$HOME/bin/iaap-cli/iaap-cli
bcftools=$BCF_TOOLS
export BCFTOOLS_PLUGINS="${BCF_TOOLS_PLUGINS_DIR}"

### Basic config
THREADS=80

###########################################################
### GSE95358 SNP-chip set
###########################################################
### 1. Convert illumina idat files to GTC files by iaap
bpm_manifest_file=${GENOTYPE_DATASETS_DIR}/GSE95358/bovinehd-manifest-b.bpm
egt_cluster_file=${GENOTYPE_DATASETS_DIR}/GSE95358/BovineHD_A.egt
sample_sheet_file=${GENOTYPE_DATASETS_DIR}/GSE95358/SampleSheet.csv
path_to_idat_folder=${GENOTYPE_DATASETS_DIR}/GSE95358/idat
GSE95358_output_folder=${SUB_OUTPUT_DIR}/GSE95358
mkdir -p ${GSE95358_output_folder}

## genotype calling from idat signal
$iaap gencall \
  $bpm_manifest_file \
  $egt_cluster_file \
  $GSE95358_output_folder \
  -s $sample_sheet_file -t $THREADS -g
  
### 2. Convert illumina GTC files to VCF files
csv_files==${GENOTYPE_DATASETS_DIR}/GSE95358/BovineHD_B1.csv
out_prefix=GSE95358
## gtc to vcf by bcftools gtc2vcf plugin
$bcftools +gtc2vcf \
  --no-version -Ou \
  -c $csv_files \
  --bpm $bpm_manifest_file \
  --egt $egt_cluster_file \
  --gtcs $GSE95358_output_folder \
  --fasta-ref $BTAU6_GENOME | \
  bcftools sort -Ou -T ${TEMP_DIR}/bcftools-sort.XXXXXX | \
  bcftools norm --no-version -Ob -o ${GSE95358_output_folder}/$out_prefix.bcf -c x -f $ref && \
  bcftools index -f ${GSE95358_output_folder}/$out_prefix.bcf
# convert bcf to vcf
bcftools view ${GSE95358_output_folder}/$out_prefix.bcf -O v -o ${GSE95358_output_folder}/$out_prefix.vcf

###########################################################
### GSE100038 SNP-chip set
###########################################################
GSE100038_output_folder=${SUB_OUTPUT_DIR}/GSE100038
### Use R script tab2vcf.R to vcf format
Rscript tab2vcf.R -o $GSE100038_output_folder

###########################################################
### BGVD genotype datasets
###########################################################
### The clean vcf was directly downloaded from BGVD database.

###########################################################
### SNP-chip sets
###########################################################
### drop duplications by plink
## GSE95358
plink --vcf ${GSE95358_output_folder}/GSE95358.vcf --list-duplicate-vars ids-only suppress-first \
--cow --const-fid GSE95358 --out ${TEMP_DIR}/GSE95358_tmp 
plink --vcf ${GSE95358_output_folder}/GSE95358.vcf --exclude ${TEMP_DIR}/GSE95358_tmp.dupvar \
--recode vcf-iid --cow --const-fid GSE95358 --out ${GSE95358_output_folder}/GSE95358_vcf_format
## GSE100038
plink --vcf ${GSE100038_output_folder}/GSE100038.vcf --make-bed --out ${GSE100038_output_folder}/GSE100038_mb
plink --bfile ${GSE100038_output_folder}/GSE100038_mb.vcf --list-duplicate-vars ids-only suppress-first \
--cow --const-fid GSE100038 --out ${TEMP_DIR}/GSE100038_tmp 
plink --bfile ${TEMP_DIR}/GSE100038_mb.vcf --exclude ${TEMP_DIR}/GSE100038_tmp.dupvar \
--recode vcf-iid --cow --const-fid GSE100038 --out ${GSE100038_output_folder}/GSE100038_vcf_format


