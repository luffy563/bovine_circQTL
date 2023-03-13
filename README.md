# BOVINE-circQTL (Breed-specific circQTL regulatory networks in bovine muscle tissue: Insights into circRNA biogenesis and multifunctional mechanisms)

***
[![](https://img.shields.io/badge/License-GPL3.0-green)](https://github.com/luffy563/bovine_circQTL/blob/main/LICENSE)
[![](https://img.shields.io/badge/Python-3.5.2-brightgreen)](https://www.python.org/downloads/release/python-352/)
![](https://img.shields.io/badge/matplotlib-3.3.3-blue)
[![](https://img.shields.io/badge/R-4.1.0-orange)](https://cloud.r-project.org/src/base/R-4/R-4.1.0.tar.gz)
![](https://img.shields.io/badge/ggplot2-3.3.5-red)

Copyright (C) 2020-2023
Northwest A&F University,
Mingzhi Liao, Xianyong Lan
  
Authors: Hongfei Liu

Contact: lhf563@nwafu.edu.cn

BOVINE-circQTL is a code repository for the pipeline of circQTL identification and downstream analysis in bovine muscle
tissue.

***
# Table of contents
***
<!--ts-->
* [BOVINE-circQTL (Breed-specific circQTL regulatory networks in bovine muscle tissue: Insights into circRNA biogenesis and multifunctional mechanisms)](#bovine-circqtl-breed-specific-circqtl-regulatory-networks-in-bovine-muscle-tissue-insights-into-circrna-biogenesis-and-multifunctional-mechanisms)
* [Table of contents](#table-of-contents)
   * [Overview of pipeline](#overview-of-pipeline)
   * [Config](#config)
   * [Required scripts](#required-scripts)
      * [Python scripts](#python-scripts)
      * [R scripts](#r-scripts)

<!-- Created by https://github.com/ekalinin/github-markdown-toc -->
<!-- Added by: luffy, at: Mon Mar 13 11:23:49 CST 2023 -->

<!--te-->
***

## Overview of pipeline
- 1.genotype_dataset_preprocess

    Raw genotype datasets contain vcf file from BGVD ([Bovine Genome Variant Database][BGVD]) 
and two SNP-chip datasets ([GSE95358](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95358) and [GSE100038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100038)).
This preprocess step includes format, imputation, and remap (from btau 6 to btau 9 genome assembly).
- 2.RNA_seq_datasets_preprocess
    
  RNA-seq datasets were retreived from from NCBI SRA database. This preprocess includes quality control, circRNA and mRNA identification.

- 3.summary_statistics_genotype
  
  This step includes summary statistics report of genotype datasets and combination of different genotype datasets from various sources.
- 4.summary_statistics_rnaseq
  
  CircRNA identification was implemented by different software, including [CircMarker][CircMarker], [CIRI2][CIRI2], and [circRNAfinder][circRNAfinder], and then the detection
results were merged in this step.
- 5.matrixeqtl
  
  CircQTL and eQTL were both identified by [MatrixEQTL][MatrixEQTL] in R basic environment. Furthermore, the functional genomic region and phenotype enrichment analysis was also carried out. 
- 6.trans-circQTL
  Substantial trans-circQTLs were found in this study, so the preliminary investigation of characteristic of tran-circQTL was involved in this step.
- 7.ABS_events
  
  First, the ABS (alternativ back-splicing) events profile was constructed to investigate the overall patterns and divergence among breeds. 
Then, to explore the association of ABS events (intron length and number/pairing ability of SINE (Short interspersed nuclear elements) /Alu-like elements) with flanking circQTL mediated by cis-elements (SINE), we extracted the flanking SINE elements closer to
ABS-circRNAs.
- 8.circRNA_functions
  
  CircRNA function mainly include miRNA sponge and RBP interaction at cytoplasm, so the potential miRNA and RBP interact with circRNAs
were predicted, by which we constructed the circRNA-miRNA/RBP interaction networks. To evaluate the distribution and effect of circQTL,
we also investigated the distribution of circQTLs within binding sites and then analyzed the change of binding ability, enrichment degree,
and secondary structure through altering the genotype of circQTLs witinin circRNAs.

## Config
- btau9.yml: a YAML-formated config file for [CIRIquant][CIRIquant] to find software and reference needed
- CM_config.ini: config file of [CircMarker][CircMarker], which contains path of reference genome fasta file, annotation file, reads1/2, and other required or optional parameters. Of which, Reference, GTF, Reads1/2, and options in Parameter section is important and required.
- software_list.txt: a list contains three circRNA detection software 
- SRR_list.txt: SRA accession ID list for RNA-seq dataset 
- SRR_breeds.csv: SRA accession ID and corresponding breed name

## Required scripts
### Python scripts
- extract_exons.py: This file is part of [HISAT 2][HISAT 2] for extracting exons from gtf annotation file.
- extract_flanking_introns.py: Extract flanking introns closer to ABS-circRNAs
- extract_splice_sites.py: This file is part of [HISAT 2][HISAT 2] for extracting splice sites from gtf annotation file.
- extract_random_seq_bed.py: Generate the custom number of random sequences from the bovine genome (ARS-UCD1.2 genome assembly) according
- filt.py: Extract the specific sequence set from a fasta file based on a sequence ID list.
- getTPM.py: Extract the TPM quantity matrix at genes and transcripts levels from the output generated by `stringtie -e`.
- get_seed_region.py: Get the seed region of miRNAs
- summary_rna_seq.py: To merge circRNA detection results and then report the summary statistics
to the specific normal distribution of sequence length.
### R scripts
- annotation.r: merge annotated circRNA list (circularRNA_known.txt) by the known circRNA
- combine.r: Combination of BGVD genotype and two SNP-chip sets
- conversion.r: convert raw genome coordinate (0-based or 1-based) to 0-based
- prep_quant.r: convert raw genome coordinate (0-based or 1-based) to 0-based for [CIRIquant][CIRIquant] calibration


[HISAT 2]:http://daehwankimlab.github.io/hisat2/
[CIRCexplorer2]: https://circexplorer2.readthedocs.io/en/latest/
[CIRI2]: https://sourceforge.net/projects/ciri/files/CIRI2/
[CIRIquant]: https://github.com/bioinfo-biols/CIRIquant
[CircMarker]:https://github.com/lxwgcool/CircMarker
[circRNAfinder]:https://github.com/bioxfu/circRNAFinder
[KNIFE]: https://github.com/blawney/knife_circ_rna
[segemehl]: https://www.bioinf.uni-leipzig.de/Software/segemehl/
[MatrixEQTL]:https://github.com/andreyshabalin/MatrixEQTL
[circBase]: http://www.circbase.org/
[circAtlas]: http://circatlas.biols.ac.cn/
[BGVD]: http://animal.omics.pro/code/index.php/BosVar