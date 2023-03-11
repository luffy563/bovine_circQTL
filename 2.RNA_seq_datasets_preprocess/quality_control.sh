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
SUB_OUTPUT_DIR=${OUTPUT_DIR}/quality_control
TEMP_DIR=${OUTPUT_DIR}/quality_control/temp

# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
fastqc=${fastqc}
multiqc=${multiqc}
fastp=${path_to_fastp}

### Basic config
THREADS=80
BTAU6_GENOME=${BTA_REFERENCE}/genome/bosTau6.fa


### QC report
mkdir -p ${SUB_OUTPUT_DIR}/QC_report_former/
$fastqc -t $THREADS -o ${SUB_OUTPUT_DIR}/QC_report_former/ ${RNASEQ_DATASETS_DIR}/fastq/*.fastq.gz

## Combine QC_report by multiqc
$multiqc ${SUB_OUTPUT_DIR}/QC_report_former/*.zip

## Combine two sparate fastq.gz files (SRR8703198 SRR8703197) to SRR87031987
cat SRR8703198.fastq.gz SRR8703197.fastq.gz > SRR87031987.fastq.gz
## QC by fastp
# sparate PE/SE files
ls ${RNASEQ_DATASETS_DIR}/fastq > ${RNASEQ_DATASETS_DIR}/fastq_files.txt
mkdir -p ${RNASEQ_DATASETS_DIR}/SE ${RNASEQ_DATASETS_DIR}/PE
for i in $(cat ${RNASEQ_DATASETS_DIR}/fastq_files.txt)
do
	if [[ $i == *_* ]];then
		mv ${RNASEQ_DATASETS_DIR}/fastq/$i ${RNASEQ_DATASETS_DIR}/fastq/PE
	else
		mv ${RNASEQ_DATASETS_DIR}/fastq/$i ${RNASEQ_DATASETS_DIR}/fastq/SE
	fi
done

# get fastq meta info
ls ${RNASEQ_DATASETS_DIR}/fastq/PE >${RNASEQ_DATASETS_DIR}/fastq/PE.txt
ls ${RNASEQ_DATASETS_DIR}/fastq/SE >${RNASEQ_DATASETS_DIR}/fastq/SE.txt

mkdir -p ${SUB_OUTPUT_DIR}/clean_fastq ${SUB_OUTPUT_DIR}/QC_report_after
# QC for PE
n=1
for fastq_file in $(cat ${RNASEQ_DATASETS_DIR}/fastq/PE.txt)
do
	total=`expr $n % 2`
	if [ $total == 1 ]
	then
		R1=$fastq_file
	else
		R2=$fastq_file
		echo $R1
		echo $R2
		$fastp -i ${RNASEQ_DATASETS_DIR}/fastq/PE/$R1 -I ${RNASEQ_DATASETS_DIR}/fastq/PE/$R2 \
		-o ${SUB_OUTPUT_DIR}/clean_fastq/$R1 -O ${SUB_OUTPUT_DIR}/clean_fastq/$R2 \
		-q 20 -u 40 -l 15 --detect_adapter_for_pe -g -x -c
	fi
	n=`expr $n + 1`
done

# QC for SE
for fastq_file in $(cat ${RNASEQ_DATASETS_DIR}/fastq/SE.txt)
do
	R1=$fastq_file
	echo $R1
	$fastp -i ${RNASEQ_DATASETS_DIR}/fastq/SE/$R1 -o ${RNASEQ_DATASETS_DIR}/clean_fastq/$R1 \
	-q 20 -u 40 -l 15 -g -x
done

### QC report
$fastqc -t $THREADS -o ${SUB_OUTPUT_DIR}/QC_report_after ${RNASEQ_DATASETS_DIR}/clean_fastq/*.fastq.gz
$multiqc  ${SUB_OUTPUT_DIR}/QC_report_after/*.zip
