#!/bin/bash

##############################################################################
# Format genotype datasets (including BGVD database and two SNP-chip sets)
# Author: Hongfei Liu
# 
##############################################################################

### Global variables
WORKING_DIR=${YOUR_WORKING_DIRECTORY}
GENOTYPE_DATASETS_DIR=${WORKING_DIR}/genotype_datasets
BTA_REFERENCE=${WORKING_DIR}/reference
OUTPUT_DIR=${WORKING_DIR}/output
SUB_OUTPUT_DIR=${OUTPUT_DIR}/combine
TEMP_DIR=${OUTPUT_DIR}/combine/temp
INPUT_DIR=${OUTPUT_DIR}/summary_statistic
SCRIPTS_DIR=${WORKING_DIR}/scripts

# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
Rscript=${Rscript}

### Basic config
THREADS=80

### Combine the genotype datasets of two SNP-chip datasets with that from BGVD database
$Rscript $SCRIPTS_DIR/combine.r -p $THREADS \
-t $INPUT_DIR/freq_table/GSE95358.frq -q $INPUT_DIR/freq_table/GSE100038.frq \
-b $GENOTYPE_DATASETS_DIR/BGVD/Btau_5.0.1_SNPs.anno.tab -o $SUB_OUTPUT_DIR/combine.vcf

