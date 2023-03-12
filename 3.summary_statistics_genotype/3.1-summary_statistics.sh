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
INPUT_DIR=${OUTPUT_DIR}/imputation

# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
PopLDdecay_DIR=$WORKING_DIR/software/PopLDdecay
PopLDdecay=$PopLDdecay_DIR/PopLDdecay
multiqc=${multiqc}
fastp=${path_to_fastp}

### Basic config
THREADS=80
BTAU6_GENOME=${BTA_REFERENCE}/genome/bosTau6.fa


###################################################
#   Summary statistics of genotype data           #
###################################################

# # # Compute the genptype frequency of snp chip datasets
mkdir -p $SUB_OUTPUT_DIR/freq_table
vcftools --gzvcf $INPUT_DIR/GSE95358/GSE95358_imp.vcf.gz --freq --out $SUB_OUTPUT_DIR/freq_table/GSE95358
vcftools --gzvcf $INPUT_DIR/GSE100038/GSE100038_imp.vcf.gz --freq --out $SUB_OUTPUT_DIR/freq_table/GSE100038

# # # Summary statistics of missing genotype/individuals of snp chip datasets
mkdir -p missing_summary
vcftools --gzvcf $INPUT_DIR/GSE95358/GSE95358_imp.vcf.gz --missing-indv --out $SUB_OUTPUT_DIR/missing_summary/GSE95358
vcftools --gzvcf $INPUT_DIR/GSE95358/GSE95358_imp.vcf.gz --missing-site --out $SUB_OUTPUT_DIR/missing_summary/GSE95358
vcftools --gzvcf $INPUT_DIR/GSE100038/GSE100038_imp.vcf.gz --missing-indv --out $SUB_OUTPUT_DIR/missing_summary/GSE100038
vcftools --gzvcf $INPUT_DIR/GSE100038/GSE100038_imp.vcf.gz --missing-site --out $SUB_OUTPUT_DIR/missing_summary/GSE100038

# # # LD decay for GSE95358
mkdir -p $SUB_OUTPUT_DIR/LD_decay

GSE95358_LDdecay_list=$SUB_OUTPUT_DIR/GSE95358_LDdecay.list

$PopLDdecay -InVCF $OUTPUT_DIR/GSE100038/GSE95358.vcf -OutStat $SUB_OUTPUT_DIR/GSE95358.stat
$PopLDdecay -InVCF $INPUT_DIR/GSE95358/GSE95358_imp.vcf -OutStat $SUB_OUTPUT_DIR/GSE95358_imp.stat
perl $PopLDdecay_DIR/bin/Plot_OnePop.pl -inFile $SUB_OUTPUT_DIR/GSE95358.stat.gz -output $SUB_OUTPUT_DIR/LD_decay/GSE95358.grpah
perl $PopLDdecay_DIR/bin/Plot_OnePop.pl -inFile $SUB_OUTPUT_DIR/GSE95358_imp.stat.gz -output $SUB_OUTPUT_DIR/LD_decay/GSE95358_imp.grpah
perl $PopLDdecay_DIR/bin/Plot_MultiPop.pl -inList $GSE95358_LDdecay_list -output $SUB_OUTPUT_DIR/LD_decay/GSE95358_compare

# # # LD decay for GSE100038
# gunzip -k ./output/imputation/GSE100038_imp.vcf.gz
GSE100038_LDdecay_list=$SUB_OUTPUT_DIR/GSE100038_LDdecay.list

$PopLDdecay -InVCF $OUTPUT_DIR/GSE100038/GSE100038.vcf -OutStat $SUB_OUTPUT_DIR/GSE100038.stat
$PopLDdecay -InVCF $INPUT_DIR/GSE100038/GSE100038_imp.vcf.gz -OutStat $SUB_OUTPUT_DIR/LD_decay/GSE100038_imp.stat
perl $PopLDdecay_DIR/bin/Plot_OnePop.pl -inFile $SUB_OUTPUT_DIR/GSE100038.stat.gz -output $SUB_OUTPUT_DIR/LD_decay/GSE100038.grpah
perl $PopLDdecay_DIR/bin/Plot_OnePop.pl -inFile $SUB_OUTPUT_DIR/LD_decay/GSE100038_imp.stat.gz -output $SUB_OUTPUT_DIR/LD decay/GSE100038_imp.grpah
perl $PopLDdecay_DIR/bin/Plot_MultiPop.pl -inList $GSE100038_LDdecay_list -output $SUB_OUTPUT_DIR/LD_decay/GSE100038_compare