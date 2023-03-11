#!/bin/bash

##############################################################################
# Remap genotype datasets from 
# Author: Hongfei Liu
# 
##############################################################################

### Global variables
WORKING_DIR=${YOUR_WORKING_DIRECTORY}
GENOTYPE_DATASETS_DIR=${WORKING_DIR}/genotype_datasets
BTA_REFERENCE=${WORKING_DIR}/reference
OUTPUT_DIR=${WORKING_DIR}/output
SUB_OUTPUT_DIR=${OUTPUT_DIR}/remap
TEMP_DIR=${OUTPUT_DIR}/remap/temp
INPUT_DIR=${OUTPUT_DIR}/imputation

# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
plink=${plink}
liftOver=${path_to_liftOver}
liftOverPlink=${path_to_liftOverPlink_py}

### Basic config
THREADS=80
BTAU6_GENOME=${BTA_REFERENCE}/genome/bosTau6.fa
bta6_to_bta9_chain=${WORKING_DIR}/liftOver_chain/bosTau6ToBosTau9.over.chain

### convert vcf to plink format
$plink --vcf ${INPUT_DIR}/GSE95358/GSE95358_filt.vcf --recode --cow --out ${SUB_OUTPUT_DIR}/GSE95358 --const-fid Holstein

$plink --vcf ${INPUT_DIR}/GSE100038/GSE100038_filt.vcf --recode --cow --out ${SUB_OUTPUT_DIR}/GSE100038 --const-fid Qinchuan

### remap by liftOverPlink
$liftOverPlink -m ${SUB_OUTPUT_DIR}/GSE95358.map -p ${SUB_OUTPUT_DIR}/GSE95358.ped \
-o $SUB_OUTPUT_DIR/GSE95358 -c $bta6_to_bta9_chain -e $liftOver

$liftOverPlink -m ${SUB_OUTPUT_DIR}/GSE100038.map -p ${SUB_OUTPUT_DIR}/GSE100038.ped \
-o $SUB_OUTPUT_DIR/GSE100038 -c $bta6_to_bta9_chain -e $liftOver
