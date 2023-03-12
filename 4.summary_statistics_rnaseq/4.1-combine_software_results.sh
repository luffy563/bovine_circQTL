#!/bin/bash

##############################################################################
# Remap genotype datasets from 
# Author: Hongfei Liu
# 
##############################################################################

### Global variables
WORKING_DIR=${YOUR_WORKING_DIRECTORY}
RNASEQ_DATASETS_DIR=${WORKING_DIR}/RNA_seq_datasets
BTA_REFERENCE=${WORKING_DIR}/reference
OUTPUT_DIR=${WORKING_DIR}/output
SUB_OUTPUT_DIR=${OUTPUT_DIR}/summary_statistic
TEMP_DIR=${OUTPUT_DIR}/summary_statistic/temp
SCRIPTS_DIR=${WORKING_DIR}/scripts


### Summary statistics of RNA-seq datasets and combine the detection results from different software
python $SCRIPTS_DIR/summary_rna_seq.py