#!/binbwa/bash

##############################################################################
# Remap genotype datasets from 
# Author: Hongfei Liu
# 
##############################################################################

### Global variables
WORKING_DIR=${YOUR_WORKING_DIRECTORY}
RNASEQ_DATASETS_DIR=${WORKING_DIR}/RNA_seq_datasets
BTA_REFERENCE=${WORKING_DIR}/reference
CONFIG_DIR=${WORKING_DIR}/config
OUTPUT_DIR=${WORKING_DIR}/output
SUB_OUTPUT_DIR=${OUTPUT_DIR}/circRNA_identification
TEMP_DIR=${OUTPUT_DIR}/circRNA_identification/temp
INPUT_DIR=${SUB_OUTPUT_DIR}/clean_fastq
SCRIPTS_DIR=${WORKING_DIR}/scripts
# create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p ${SUB_OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

### Required software
bwa=${bwa}
STAR=${STAR}
hisat2_build=${hisat2-build}
CIRCexplorer2=${CIRCexplorer2}
CIRI2_path=$WORKING_DIR/software/CIRI_v2.0.6
CircMarker_path=${WORKING_DIR}/software/CircMarker-master/CircRnaDetectDraft/CircRNA_Detect_Draft/MakeFile
circRNA_finder_path=$WORKING_DIR/software/circRNA_finder-master
CIRIquant=${CIRIquant}


### Basic config
THREADS=80
SE_RNA_SEQ_META=${RNASEQ_DATASETS_DIR}/fastq/SE.txt
PE_RNA_SEQ_META=${RNASEQ_DATASETS_DIR}/fastq/PE.txt
BTAU6_GENOME=${BTA_REFERENCE}/genome/bosTau6.fa
BOS_REF_BED=${BTA_REFERENCE}/genome/bos_ref.txt
ARS_UCD_GENOME=${BTA_REFERENCE}/genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
ARS_UCD_GTF=${BTA_REFERENCE}/annotation/Bos_taurus.ARS-UCD1.2.101.gtf

### Building HISAT2 index for bovine
mkdir -p ${BTA_REFERENCE}/ht2_index
$SCRIPTS_DIR/extract_splice_sites.py $ARS_UCD_GTF >${BTA_REFERENCE}/ht2_index/btau9.ss 
$SCRIPTS_DIR/extract_exons.py $ARS_UCD_GTF >${BTA_REFERENCE}/ht2_index/btau9.exon 
$hisat2_build -p $THREADS --ss ${BTA_REFERENCE}/ht2_index/btau9.ss --exon ${BTA_REFERENCE}/ht2_index/btau9.exon \
$ARS_UCD_GENOME ${BTA_REFERENCE}/ht2_index/btau9_tran 

######################################
# Identification by CircMarker
######################################
### CircMarker
# For PE
cat $PE_RNA_SEQ_META | while read name len
do
	mkdir -p $SUB_OUTPUT_DIR/CircMarker/$name
	loc=${INPUT_DIR}\/$name
	echo $loc
	cp ${CONFIG_DIR}/CM_config.ini $path/config.ini
	sed -i "/Reads1=/s/$/${INPUT_DIR}\/${name}_1.fastq.gz/" $path/config.ini
	sed -i "/Reads2=/s/$/${INPUT_DIR}\/${name}_2.fastq.gz/" $path/config.ini
	sed -i "/ReadsLen=/s/$/${len}/" $path/config.ini
	$path/CircRnaDetectDraft $path/config.ini
	cp -r Detection_Result $SUB_OUTPUT_DIR/CircMarker/$name
done
# For SE
cat $SE_RNA_SEQ_META | while read name len
do
	mkdir -p $SUB_OUTPUT_DIR/CircMarker/$name
	loc=${INPUT_DIR}\/$name
	echo $loc
	cp ${CONFIG_DIR}/CM_config.ini $path/config.ini
	sed -i "/Reads1=/s/$/${INPUT_DIR}\/${name}.fastq.gz/" $path/config.ini
	sed -i "/ReadsLen=/s/$/${len}/" $path/config.ini
	$path/CircRnaDetectDraft $path/config.ini
	cp -r Detection_Result $SUB_OUTPUT_DIR/CircMarker/$name
done
## Conversion
Rscript $SCRIPTS_DIR/conversion.r -a CM -s 1
Rscript $SCRIPTS_DIR/conversion.r -a CM -s 0
## Annotate
cat $PE_RNA_SEQ_META | while read i len
do
$CIRCexplorer2 annotate -r $BOS_REF_BED -g $ARS_UCD_GENOME \
-b $SUB_OUTPUT_DIR/CircMarker/${i}/circ_candidates_convert.bed -o $SUB_OUTPUT_DIR/CircMarker/${i}/circularRNA_known.txt
done
cat $SE_RNA_SEQ_META | while read i
do
$CIRCexplorer2 annotate -r $BOS_REF_BED -g $ARS_UCD_GENOME \
-b $SUB_OUTPUT_DIR/CircMarker/${i}/circ_candidates_convert.bed -o $SUB_OUTPUT_DIR/CircMarker/${i}/circularRNA_known.txt
done


### CircMarker quantification results calibrated by CIRIquant quantification
# Rscript $SCRIPTS_DIR/prep_quant.r -a CM -s 1

## For PE only
# cat $PE_RNA_SEQ_META | while read name len
# do
# $CIRIquant -t $THREADS \
          # -1 ${INPUT_DIR}/${name}_1.fastq.gz \
          # -2 ${INPUT_DIR}/${name}_2.fastq.gz \
          # --config ${CONFIG_DIR}/btau9.yml \
          # -o $SUB_OUTPUT_DIR/CircMarker/$name \
          # -p ${name}_quant \
          # --bed $SUB_OUTPUT_DIR/CircMarker/$name/circ_candidates_CIRIquant.bed
# done


######################################
# Identification by CIRI2
######################################	  
# CIRI2 prediction

# BWA indexing
$bwa index $ARS_UCD_GENOME
# Alignment and identification for PE
cat $PE_RNA_SEQ_META | while read name len
do
mkdir -p $SUB_OUTPUT_DIR/CIRI2/${name}/Alignment
$bwa mem -t $THREADS -T 19 $ARS_UCD_GENOME ${INPUT_DIR}/${name}_1.fastq.gz ${INPUT_DIR}/${name}_2.fastq.gz \
> $SUB_OUTPUT_DIR/CIRI2/${name}/Alignment/aln-pe.sam
perl $CIRI2_path/CIRI2.pl -I $SUB_OUTPUT_DIR/CIRI2/${name}/Alignment/aln-pe.sam -O $SUB_OUTPUT_DIR/CIRI2/${name}/out.ciri \
-F $ARS_UCD_GENOME -A $ARS_UCD_GTF
done
# Alignment and identification for SE
cat $SE_RNA_SEQ_META | while read name
do
mkdir -p output/CIRI2/${name}/Alignment
$bwa mem -t $THREADS -T 19 $ARS_UCD_GENOME ${INPUT_DIR}/${name}.fastq.gz > $SUB_OUTPUT_DIR/CIRI2/${name}/Alignment/aln-se.sam
perl $CIRI2_path/CIRI2.pl -I $SUB_OUTPUT_DIR/CIRI2/${name}/Alignment/aln-se.sam -O $SUB_OUTPUT_DIR/CIRI2/${name}/out.ciri \
-F $ARS_UCD_GENOME -A $ARS_UCD_GTF
done

### Delete Alignment temp files of CIRI2 (Optional)
# for i in $(cat $SE_RNA_SEQ_META)
# do
	# mkdir -p temp
	# mkdir -p $SUB_OUTPUT_DIR/CIRI2/$i/Alignment
	# rm -r ./output/CIRI2/$i/Alignment 
# done

# for i in $(cat $PE_RNA_SEQ_META)
# do
	# mkdir -p temp
	# mkdir -p $SUB_OUTPUT_DIR/CIRI2/$i/Alignment
	# rm -r $SUB_OUTPUT_DIR/CIRI2/$i/Alignment
# done

# # Conversion
Rscript $SCRIPTS_DIR/conversion.r -a CIRC2 -s 1
# # Annotate
cat $PE_RNA_SEQ_META | while read i len
do
$CIRCexplorer2 annotate -r $BOS_REF_BED -g $ARS_UCD_GENOME \
-b $SUB_OUTPUT_DIR/CIRI2/${i}/circ_candidates_convert.bed -o $SUB_OUTPUT_DIR/CIRI2/${i}/circularRNA_known.txt
done
cat $SE_RNA_SEQ_META | while read i
do
$CIRCexplorer2 annotate -r $BOS_REF_BED -g $ARS_UCD_GENOME \
-b $SUB_OUTPUT_DIR/CIRI2/${i}/circ_candidates_convert.bed -o $SUB_OUTPUT_DIR/CIRI2/${i}/circularRNA_known.txt
done

### CIRI2 CIRIquant quantification (Required for Paired-end mRNA identification and quantity)
Rscript $SCRIPTS_DIR/prep_quant.r -a CIRI2 -s 1

cat $PE_RNA_SEQ_META | while read i len
do
$CIRIquant -t $THREADS \
          -1 $INPUT_DIR/${i}_1.fastq.gz \
          -2 $INPUT_DIR/${i}_2.fastq.gz \
          --config $CONFIG_DIR/btau9.yml \
          -o $SUB_OUTPUT_DIR/CIRI2/$i \
          -p ${i}_quant \
		  --tool CIRI2 \
		  --circ $SUB_OUTPUT_DIR/CIRI2/$i/out.ciri
done

### Single-end mRNA-seq
for fastq_file in $(cat $SE_RNA_SEQ_META)
do
mkdir -p $SUB_OUTPUT_DIR/CIRI2/${fastq_file}/Alignment/
mkdir -p $SUB_OUTPUT_DIR/CIRI2/${fastq_file}/gene/
R1=${fastq_file}
echo $R1
$hisat2 -p $THREADS --dta -x $BTA_REFERENCE/ht2_index/btau9_tran -U $INPUT_DIR/$R1.fastq.gz | $samtools view -hbuS - | \
$samtools sort -o $SUB_OUTPUT_DIR/CIRI2/${R1}/Alignment/${R1}.bam
# Assemble
$stringtie $SUB_OUTPUT_DIR/CIRI2/${R1}/Alignment/${R1}.bam -p $THREADS \
-G $ARS_UCD_GTF -o $SUB_OUTPUT_DIR/CIRI2/${R1}/gene/${R1}.gtf \
-l ${R1} $SUB_OUTPUT_DIR/CIRI2/${R1}/Alignment/${R1}.bam
done

for bam_file in $(cat $SE_RNA_SEQ_META)
do
echo $bam_file
mkdir -p $SUB_OUTPUT_DIR/CIRI2/ballgown/$bam_file
stringtie -e -B -p $THREADS -G $ARS_UCD_GTF \
-o $SUB_OUTPUT_DIR/CIRI2/ballgown/$bam_file/${bam_file}.gtf $SUB_OUTPUT_DIR/CIRI2/${bam_file}/Alignment/$bam_file.bam
done
python $SCRIPTS_DIR/getTPM.py -i $SUB_OUTPUT_DIR/CIRI2/ballgown/ \
-g $SUB_OUTPUT_DIR/CIRI2/gene_tpm_matrix.csv -t $SUB_OUTPUT_DIR/CIRI2/transcript_tpm_matrix.csv


######################################
# Identification by CircRNAfinder
######################################
# CircRNAfinder
# indexing
mkdir -p $BTA_REFERENCE/star_index
$STAR --runThreadN 30 --runMode genomeGenerate --genomeDir $BTA_REFERENCE/star_index \
--genomeFastaFiles $ARS_UCD_GENOME --sjdbGTFfile $ARS_UCD_GTF --sjdbOverhang 99

# Alignment FOR PE
cat $PE_RNA_SEQ_META | while read i len
do

mkdir -p $SUB_OUTPUT_DIR/CircRNAfinder/${i}/Alignment
$STAR --genomeDir $BTA_REFERENCE/star_index \
--readFilesCommand zcat \
--readFilesIn $INPUT_DIR/${i}_1.fastq.gz $INPUT_DIR/${i}_2.fastq.gz \
--readMatesLengthsIn NotEqual \
--runThreadN 15 \
--chimSegmentMin 20 \
--chimScoreMin 1 \
--alignIntronMax 500000 \
--outFilterMismatchNmax 4 \
--alignTranscriptsPerReadNmax 100000 \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 2 \
--outFileNamePrefix $SUB_OUTPUT_DIR/CircRNAfinder/${i}/Alignment/${i}

perl $circRNA_finder_path/postProcessStarAlignment.pl --starDir $SUB_OUTPUT_DIR/CircRNAfinder/${i}/Alignment/${i} \
--minLen 100 --outDir $SUB_OUTPUT_DIR/CircRNAfinder/${i}/
done
## Accurate quantity by CIRIquant to CircRNAfinder
# Rscript $SCRIPTS_DIR/prep_quant.r -a CF -s 1

# cat $PE_RNA_SEQ_META | while read name len
# do  
# $CIRIquant -t $THREADS \
          # -1 $INPUT_DIR/${name}_1.fastq.gz \
          # -2 $INPUT_DIR/${name}_2.fastq.gz \
          # --config $CONFIG_DIR/btau9.yml \
          # -o $SUB_OUTPUT_DIR/CircRNAfinder/$name \
          # -p ${name}_quant \
		  # --tool circRNA_finder \
		  # --circ $SUB_OUTPUT_DIR/CircRNAfinder/${name}/${name}s_filteredJunctions.bed
# done


# Alignment FOR SE
cat $SE_RNA_SEQ_META | while read i
do
mkdir -p $SUB_OUTPUT_DIR/CircRNAfinder/${i}/Alignment
$STAR --genomeDir $BTA_REFERENCE/star_index/ \
--readFilesCommand zcat \
--readFilesIn $INPUT_DIR/${i}.fastq.gz \
--readMatesLengthsIn NotEqual \
--runThreadN 15 \
--chimSegmentMin 20 \
--chimScoreMin 1 \
--alignIntronMax 500000 \
--outFilterMismatchNmax 4 \
--alignTranscriptsPerReadNmax 100000 \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 2 \
--outFileNamePrefix $SUB_OUTPUT_DIR/CircRNAfinder/${i}/Alignment/${i}

perl $path/postProcessStarAlignment.pl --starDir $SUB_OUTPUT_DIR/CircRNAfinder/${i}/Alignment/${i} \
--minLen 100 --outDir $SUB_OUTPUT_DIR/CircRNAfinder/${i}/
done
