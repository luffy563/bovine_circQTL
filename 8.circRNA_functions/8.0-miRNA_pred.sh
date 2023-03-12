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
SUB_OUTPUT_DIR=${OUTPUT_DIR}/miRNA_pred
TEMP_DIR=${OUTPUT_DIR}/miRNA_pred/temp
SCRIPTS_DIR=${WORKING_DIR}/scripts

# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
Rscript=${Rscript}
vienna_home=/home/liuhongfei/miniconda3/pkgs/viennarna-2.4.18-py27h8eb80aa_0/share/ViennaRNA/bin
ssw_test=$WORKING_DIR/software/Complete-Striped-Smith-Waterman-Library/src/ssw_test
IntaRNA=$IntaRNA
seqkit=$seqkit
RNAfold=$RNAfold

### Basic config
THREADS=80
genome=$BTA_REFERENCE/genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa

###########################################################################################
#   miRNA binding sites predicted by Targetscan                                           #
###########################################################################################
## miRNA sites predict

# get miR_family_info & miRNA_info  8094
# get all miRNA family info from fasta file of miRBase
python $SCRIPTS_DIR/get_seed_region.py

awk '{if ($3=="9913") print $1"\t"$2"\t"$3}' $BTA_REFERENCE/miRNA/miR_Family_Info_convert.txt>.$BTA_REFERENCE/miRNA/bta_miR_family_info.txt
cat $BTA_REFERENCE/miRNA/bta_miR_family_info.txt | while read miRNAFamily seed taxID
do
# awk '{if (($1 ~ /^'"$miRNAFamily"'$/) || ($1 ~ /^'"$miRNAFamily"'\//)) print $1"\t"$3"\t"$4"\t"$5}' $BTA_REFERENCE/miRNA/miR_Family_Info.txt>>miRNA_for_context.txt
awk '{if (($1 ~ /^'"$miRNAFamily"'$/)) print $1"\t"$3"\t"$4"\t"$5}' $BTA_REFERENCE/miRNA/miR_Family_Info.txt>>$BTA_REFERENCE/miRNA/miRNA_for_context.txt

done

# classify all miRNA family site type in cow
cd $WORKING_DIR/software/targetscan
awk '{if($1 != "id") print $1"\t"9913"\t"$9}' $SUB_OUTPUT_DIR/targetsFTS_circ_seq.txt
awk '{if($1 != "id") print $1"\t"9913"\t"$9}' $SUB_OUTPUT_DIR/targetsFTS_circ_seq.txt >$SUB_OUTPUT_DIR/circRNA_internal_seq.txt
awk '{if(($1 != "id") && ($12 ~ "alt")) print $1"\t"9913"\t"$11}' $SUB_OUTPUT_DIR/alt_genotype_circseq.txt >$SUB_OUTPUT_DIR/alt_circRNA_internal_seq.txt

perl targetscan_70.pl $BTA_REFERENCE/miRNA/all_bta_miR_family_info.txt $SUB_OUTPUT_DIR/circRNA_internal_seq.txt \
$SUB_OUTPUT_DIR/all_miRNA_binding_sites_internal_by_targetscan.txt

perl targetscan_70.pl $BTA_REFERENCE/miRNA/all_bta_miR_family_info.txt $SUB_OUTPUT_DIR/alt_circRNA_internal_seq.txt \
$SUB_OUTPUT_DIR/all_miRNA_binding_sites_alt_internal_by_targetscan.txt

# classify conserved miRNA family site type
perl targetscan_70.pl $BTA_REFERENCE/miRNA/bta_miR_family_info.txt \
$SUB_OUTPUT_DIR/circRNA_internal_seq.txt .$SUB_OUTPUT_DIR/conserved_miRNA_binding_sites_internal_by_targetscan.txt

## Calculate free Energy by RNAfold and aligned Score by SW algorithm by EMBOSS
# get all circRNA and miRNA fasta files
# miRNA fasta
# $BTA_REFERENCE/miRNA/bos_mature.fa
# circRNA fasta
cat $SUB_OUTPUT_DIR/circRNA_internal_seq.txt | while read circ taxID seq
do
if [[ "${seq}" =~ "circ" ||  "${seq}" =~ "NA" ]]
then
	echo "OK"
	continue

else

	echo \>$circ>>$SUB_OUTPUT_DIR/circRNA_internal_seq.fa
	echo $seq>>$SUB_OUTPUT_DIR/circRNA_internal_seq.fa
fi
done

## Alt circRNA internal 
cat $SUB_OUTPUT_DIR/alt_circRNA_internal_seq.txt | while read circ taxID seq
do
if [[ "${seq}" =~ "circ" ||  "${seq}" =~ "NA" ]]
then
	echo "OK"
	continue

else

	echo \>$circ>>$SUB_OUTPUT_DIR/alt_circRNA_internal_seq.fa
	echo $seq>>$SUB_OUTPUT_DIR/alt_circRNA_internal_seq.fa
fi
done


# calculate mfw by RNAflod through InaRNA
$IntaRNA --threads 20 -m H --outMode C \
-q $SUB_OUTPUT_DIR/circRNA_internal_seq.fa -t $BTA_REFERENCE/miRNA/bos_mature.fa \
--out $SUB_OUTPUT_DIR/test_byIntaRNA.csv

# ref circRNA
awk '{print $1}' $SUB_OUTPUT_DIR/alt_circRNA_internal_seq.txt > $SUB_OUTPUT_DIR/alt_circRNA.list
python $SCRIPTS_DIR/filt.py -i $SUB_OUTPUT_DIR/alt_circRNA.list -f $SUB_OUTPUT_DIR/circRNA_internal_seq.fa -o $SUB_OUTPUT_DIR/ref_circRNA_internal_seq.fa

python $SCRIPTS_DIR/filt.py -i $SUB_OUTPUT_DIR/alt_circRNA.list -f $SUB_OUTPUT_DIR/circRNA_internal_seq.fa -o $SUB_OUTPUT_DIR/ref_circRNA_internal_seq.fa
$IntaRNA --threads 20 -m H --outMode C \
-q $SUB_OUTPUT_DIR/ref_circRNA_internal_seq.fa -t $BTA_REFERENCE/miRNA/bos_mature.fa \
--out $SUB_OUTPUT_DIR/ref_byIntaRNA.csv
# alt circRNA
$IntaRNA --threads 20 -m H --outMode C \
-q $SUB_OUTPUT_DIR/alt_circRNA_internal_seq.fa -t $BTA_REFERENCE/miRNA/bos_mature.fa \
--out $SUB_OUTPUT_DIR/alt_byIntaRNA.csv

# aligned Score by SW algorithm through ssw_test
awk -F ";" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"}' $SUB_OUTPUT_DIR/test_byIntaRNA.csv>$SUB_OUTPUT_DIR/test_byIntaRNA.txt
cat $SUB_OUTPUT_DIR/test_byIntaRNA.txt | while read miR a b c d miRstart miRend circ start2 end2 subseqDP hybridDP E
do
echo $miR >$TEMP_DIR/tempmiR.list;echo $circ >$TEMP_DIR/tempcirc.list
python $SCRIPTS_DIR/filt.py -i $TEMP_DIR/tempcirc.list -f $SUB_OUTPUT_DIR/circRNA_internal_seq.fa -o $SUB_OUTPUT_DIR/mRNA_raw.fa
# start=`$start2 - 10`;end=`$end2 + 10`
$seqkit subseq -r ${start2}:${end2} $TEMP_DIR/mRNA_raw.fa >$TEMP_DIR/mRNA-1.fa 
python $SCRIPTS_DIR/filt.py -i $TEMP_DIR/tempmiR.list -f $BTA_REFERENCE/miRNA/bos_mature.fa -o $TEMP_DIR/miR.fa
$seqkit seq -r -p -t RNA $TEMP_DIR/mRNA-1.fa >$TEMP_DIR/mRNA_com.fa
# ./software/Complete-Striped-Smith-Waterman-Library/src/ssw_test \
# -m2 -x2 -o2 -e2 -c $SUB_OUTPUT_DIR/miR.fa $SUB_OUTPUT_DIR/mRNA_com.fa
# get the aligned score of each miRNA
alignedScore=$($ssw_test \
-m2 -x2 -o2 -e1 -r $TEMP_DIR/miR.fa $TEMP_DIR/mRNA_com.fa | grep 'optimal_alignment_score: ' \
| grep -oP 'optimal_alignment_score: \d+' | grep -oP '\d+')
echo ${miR}$'\t'${miRstart}$'\t'${miRend}$'\t'${circ}$'\t'${E}$'\t'${alignedScore} >>$SUB_OUTPUT_DIR/miRNA_bySW.txt
done
done

### ref by SW
awk -F ";" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"}' $SUB_OUTPUT_DIR/ref_byIntaRNA.csv>$SUB_OUTPUT_DIR/ref_byIntaRNA.txt
cat $SUB_OUTPUT_DIR/ref_byIntaRNA.txt | while read miR a b c d miRstart miRend circ start2 end2 subseqDP hybridDP E
do
echo $miR >$TEMP_DIR/tempmiR.list;echo $circ >$TEMP_DIR/tempcirc.list
python $SCRIPTS_DIR/filt.py -i $TEMP_DIR/tempcirc.list -f $SUB_OUTPUT_DIR/ref_circRNA_internal_seq.fa -o ./mRNA_raw.fa
# start=`$start2 - 10`;end=`$end2 + 10`
$seqkit subseq -r ${start2}:${end2} $TEMP_DIR/mRNA_raw.fa >$TEMP_DIR/mRNA-1.fa 
python $SCRIPTS_DIR/filt.py -i $TEMP_DIR/tempmiR.list -f $BTA_REFERENCE/miRNA/bos_mature.fa -o $TEMP_DIR/miR.fa
$seqkit seq -r -p -t RNA ./mRNA-1.fa >./mRNA_com.fa
# ./software/Complete-Striped-Smith-Waterman-Library/src/ssw_test \
# -m2 -x2 -o2 -e2 -c ./miR.fa ./mRNA_com.fa
# get the aligned score of each miRNA
alignedScore=$(./software/Complete-Striped-Smith-Waterman-Library/src/ssw_test \
-m2 -x2 -o2 -e1 -r $TEMP_DIR/miR.fa $TEMP_DIR/mRNA_com.fa | grep 'optimal_alignment_score: ' \
| grep -oP 'optimal_alignment_score: \d+' | grep -oP '\d+')
echo ${miR}$'\t'${miRstart}$'\t'${miRend}$'\t'${circ}$'\t'${E}$'\t'${alignedScore} >>$SUB_OUTPUT_DIR/ref_test_bySW.txt
done
done


### alt genotype by SW
awk -F ";" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"}' $SUB_OUTPUT_DIR/alt_byIntaRNA.csv>$SUB_OUTPUT_DIR/alt_byIntaRNA.txt
cat $SUB_OUTPUT_DIR/alt_byIntaRNA.txt | while read miR a b c d miRstart miRend circ start2 end2 subseqDP hybridDP E
do
echo $miR >$TEMP_DIR/tempmiR.list;echo $circ >$TEMP_DIR/tempcirc.list
python $SCRIPTS_DIR/filt.py -i $TEMP_DIR/tempcirc.list -f $SUB_OUTPUT_DIR/alt_circRNA_internal_seq.fa -o $TEMP_DIR/mRNA_raw.fa
# start=`$start2 - 10`;end=`$end2 + 10`
$seqkit subseq -r ${start2}:${end2} $TEMP_DIR/mRNA_raw.fa >$TEMP_DIR/mRNA-1.fa 
python $SCRIPTS_DIR/filt.py -i $TEMP_DIR/tempmiR.list -f $BTA_REFERENCE/miRNA/bos_mature.fa -o $TEMP_DIR/miR.fa
$seqkit seq -r -p -t RNA ./mRNA-1.fa >./mRNA_com.fa
# ./software/Complete-Striped-Smith-Waterman-Library/src/ssw_test \
# -m2 -x2 -o2 -e2 -c ./miR.fa ./mRNA_com.fa
# get the aligned score of each miRNA
alignedScore=$($ssw_test \
-m2 -x2 -o2 -e1 -r $TEMP_DIR/miR.fa $TEMP_DIR/mRNA_com.fa | grep 'optimal_alignment_score: ' \
| grep -oP 'optimal_alignment_score: \d+' | grep -oP '\d+')
echo ${miR}$'\t'${miRstart}$'\t'${miRend}$'\t'${circ}$'\t'${E}$'\t'${alignedScore} >>$SUB_OUTPUT_DIR/alt_test_bySW.txt
done

### RBP motif search
# generate random seqs
python extract_random_seq.py
$bedtools getfasta -fi $genome -bed $BTA_REFERENCE/miRNA/random_seq.bed -fo $BTA_REFERENCE/miRNA/random_seq.fasta

#####################################################################################################################
###########################################################################################
#   circRNA secondary structure predicted by RNAfold                                      #
###########################################################################################
## To search the RNA secondary structure between ref and alt
## Secondary structure prediction by RNAfold
mkdir -p $OUTPUT_DIR/RNAstructure/ref
c -v -j20 -p -c -d2 --noLP --noClosingGU <$SUB_OUTPUT_DIR/ref_circRNA_internal_seq.fa >>$OUTPUT_DIR/RNAstructure/ref/ref_rna_structure.txt
cp *.ps *.fold $OUTPUT_DIR/RNAstructure/ref
mkdir -p $OUTPUT_DIR/RNAstructure/alt
$RNAfold -v -j20 -p -c -d2 --noLP --noClosingGU <$SUB_OUTPUT_DIR/alt_circRNA_internal_seq.fa >>$OUTPUT_DIR/RNAstructure/alt/alt_rna_structure.txt
cp *.ps *.fold $OUTPUT_DIR/RNAstructure/alt
## output the mountain plot
ls $OUTPUT_DIR/RNAstructure/ref | grep ss.ps >$OUTPUT_DIR/RNAstructure/all_circRNA.list
cat $OUTPUT_DIR/RNAstructure/all_circRNA.list | while read circ_ss
do
circ=${circ_ss:0:-6}
# add the colored-annotation of base-pair probability 
# perl $vienna_home/relplot.pl ZNF536_+_chr18_40972964_40975136_ss.ps ZNF536_+_chr18_40972964_40975136_dp.ps >ZNF536_+_chr18_40972964_40975136_colored.ps
perl $vienna_home/relplot.pl $OUTPUT_DIR/RNAstructure/ref/${circ}_ss.ps $OUTPUT_DIR/RNAstructure/ref/${circ}_dp.ps >$OUTPUT_DIR/RNAstructure/ref/${circ}_colored.ps
perl $vienna_home/relplot.pl $OUTPUT_DIR/RNAstructure/alt/${circ}_ss.ps $OUTPUT_DIR/RNAstructure/alt/${circ}_dp.ps >$OUTPUT_DIR/RNAstructure/alt/${circ}_colored.ps
# get the mfe, pf-score, and entropy score of each position 
perl $vienna_home/cmount.pl
perl $vienna_home/mountain.pl $OUTPUT_DIR/RNAstructure/ref/${circ}_dp.ps >$OUTPUT_DIR/RNAstructure/ref/${circ}_score.txt
perl $vienna_home/mountain.pl $OUTPUT_DIR/RNAstructure/alt/${circ}_dp.ps >$OUTPUT_DIR/RNAstructure/alt/${circ}_score.txt

done
## split two different structures (mfe and centroid) 
# ref
rm $OUTPUT_DIR/RNAstructure/ref/mfe_structure.txt
rm $OUTPUT_DIR/RNAstructure/ref/therm_structure.txt
rm $OUTPUT_DIR/RNAstructure/ref/cent_structure.txt
cat $OUTPUT_DIR/RNAstructure/ref/ref_rna_structure.txt | while read line
do
	if [[ $line =~ ">" ]]
	then
	n=0
	echo $line>>$OUTPUT_DIR/RNAstructure/ref/mfe_structure.txt
	echo $line>>$OUTPUT_DIR/RNAstructure/ref/therm_structure.txt
	echo $line>>$OUTPUT_DIR/RNAstructure/ref/cent_structure.txt
	elif [[ $line =~ ".(" ]]
	then
	n=`expr $n + 1`
	if [[ $n == 1 ]]
	then
	echo $line>>$OUTPUT_DIR/RNAstructure/ref/mfe_structure.txt
	elif [[ $n == 2 ]]
	then
	echo $line>>$OUTPUT_DIR/RNAstructure/ref/therm_structure.txt
	elif [[ $n == 3 ]]
	then
	echo $line>>$OUTPUT_DIR/RNAstructure/ref/cent_structure.txt
	fi
	fi
done
# alt
rm $OUTPUT_DIR/RNAstructure/alt/mfe_structure.txt
rm $OUTPUT_DIR/RNAstructure/alt/therm_structure.txt
rm $OUTPUT_DIR/RNAstructure/alt/cent_structure.txt
cat $OUTPUT_DIR/RNAstructure/alt/alt_rna_structure.txt | while read line
do
	if [[ $line =~ ">" ]]
	then
	n=0
	echo $line>>$OUTPUT_DIR/RNAstructure/alt/mfe_structure.txt
	echo $line>>$OUTPUT_DIR/RNAstructure/alt/therm_structure.txt
	echo $line>>$OUTPUT_DIR/RNAstructure/alt/cent_structure.txt
	elif [[ $line =~ ".(" ]]
	then
	n=`expr $n + 1`
	if [[ $n == 1 ]]
	then
	echo $line>>$OUTPUT_DIR/RNAstructure/alt/mfe_structure.txt
	elif [[ $n == 2 ]]
	then
	echo $line>>$OUTPUT_DIR/RNAstructure/alt/therm_structure.txt
	elif [[ $n == 3 ]]
	then
	echo $line>>$OUTPUT_DIR/RNAstructure/alt/cent_structure.txt
	fi
	fi
done