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
SUB_OUTPUT_DIR=${OUTPUT_DIR}/Alu_iden
TEMP_DIR=${OUTPUT_DIR}/Alu_iden/temp
INPUT_DIR=${OUTPUT_DIR}/summary_statistic
SCRIPTS_DIR=${WORKING_DIR}/scripts

# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
convert2bed=${convert2bed}
bedtools=${bedtools}
blastn=${blastn}

### Configs
GENOME_FASTA_FILE=$BTA_REFERENCE/genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
GENE_ANNOT_FILE=$BTA_REFERENCE/genome/Bos_taurus.ARS-UCD1.2.101.gtf

### Identification of Alu-like elements
# Extract introns bed annotation
# Extract genes
awk '$3 == "gene"' $GENE_ANNOT_FILE > $SUB_OUTPUT_DIR/genes_temp.gtf

# Use BEDops convert2bed
# convert2bed won't work if the "transcript_id" field is not there. This makes an empty "transcript_id" field so convert2bed is happy.
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' $SUB_OUTPUT_DIR/genes_temp.gtf\
 > $SUB_OUTPUT_DIR/genes.gtf

#Extract exons (exons always have transcript_id field, no need to repeat above)
awk '$3 == "exon"' $GENE_ANNOT_FILE > $SUB_OUTPUT_DIR/exons.gtf

#Convert both to bed format so we can use bedtools
$convert2bed -i gtf < $SUB_OUTPUT_DIR/exons.gtf > $SUB_OUTPUT_DIR/exons.bed
$convert2bed -i gtf < $SUB_OUTPUT_DIR/genes.gtf > $SUB_OUTPUT_DIR/genes.bed

# Remove exons from the genes.bed file, enforce strand specificy (-s)
$bedtools subtract -a $SUB_OUTPUT_DIR/genes.bed -b $SUB_OUTPUT_DIR/exons.bed -s > $SUB_OUTPUT_DIR/introns.bed

# Convert intron bed file to GTF file
awk '{print $1"\t"$7"\t""intron""\t"($2+1)"\t"$3"\t"$5"\t"$6"\t"$9"\t"(substr($0, index($0,$10)))}' $SUB_OUTPUT_DIR/introns.bed \
> $SUB_OUTPUT_DIR/introns.gtf

# Get the fasta file of introns
$bedtools getfasta -fi $GENOME_FASTA_FILE -bed $SUB_OUTPUT_DIR/introns.bed -fo $SUB_OUTPUT_DIR/introns.fasta

## filt the SINE from rmsk.txt
awk '$12 == "SINE"' $SUB_OUTPUT_DIR/rmsk.txt > $SUB_OUTPUT_DIR/bt9_SINE.txt
awk '{print $6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13}' $SUB_OUTPUT_DIR/bt9_SINE.txt > $SUB_OUTPUT_DIR/bt9_SINE.bed

# get the flanking (upstream or downstream) alu seqs of each circRNA (for preaparing)
$bedtools getfasta -name -fi $GENOME_FASTA_FILE -bed $SUB_OUTPUT_DIR/bt9_SINE.bed -fo $SUB_OUTPUT_DIR/alu_like.fasta

$bedtools intersect -a $SUB_OUTPUT_DIR/bt9_SINE.bed -b $SUB_OUTPUT_DIR/introns.gtf -wb > $SUB_OUTPUT_DIR/bt9_SINE.gtf


#### Blast for calculating CSI (Conplementary Stability Index) of Alu-like elements in flanking introns of circRNAs with ABS
# extract flanking Alu like elements to circRNAs
Rscript extract_flanking_Alu.r

# extract query and subject seqs of Alu-like elements to files
cat $SUB_OUTPUT_DIR/alu_flank.txt | while read circ query subject
do
	echo $query
	echo $subject
	var=${query//,/ }
	var1=${subject//,/ }
	rm $SUB_OUTPUT_DIR/query_1.list && rm $SUB_OUTPUT_DIR/subject_1.list
	for element in $var   
	do  
		echo $element>>$SUB_OUTPUT_DIR/query_1.list
	done  
	for element in $var1  
	do  
		echo $element>>$SUB_OUTPUT_DIR/subject_1.list
	done
	python $SCRIPTS_DIR/filt.py -i $SUB_OUTPUT_DIR/query_1.list -f $SUB_OUTPUT_DIR/alu_like.fasta -o $SUB_OUTPUT_DIR/query_1.fasta
	python $SCRIPTS_DIR/filt.py -i $SUB_OUTPUT_DIR/subject_1.list -f $SUB_OUTPUT_DIR/alu_like.fasta -o $SUB_OUTPUT_DIR/subject_1.fasta

	$blastn -query $SUB_OUTPUT_DIR/query_1.fasta -subject $SUB_OUTPUT_DIR/subject_1.fasta -evalue 1e-5 \
	-outfmt 6 -out $TEMP_DIR/blastn_output_total.txt
	size=$(ls -l $TEMP_DIR/blastn_output_total.txt)
	score=0
	if [ ${size:34:2} == 0 ]
	then
		echo $circ "NA" >>$SUB_OUTPUT_DIR/blast_output_total.txt
	else
		awk '{print "'"$circ"'\t" $0}' $TEMP_DIR/blastn_output_total.txt > $TEMP_DIR/blastn_output_tmp.txt

		awk '{print $1"\t"$2"\t"$3"\t"$13}' $TEMP_DIR/blastn_output_tmp.txt >>$SUB_OUTPUT_DIR/blast_output_total.txt

	fi
done
