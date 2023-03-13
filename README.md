# BOVINE-circQTL (Breed-specific circQTL regulatory networks in bovine muscle tissue: Insights into circRNA biogenesis and multifunctional mechanisms)

***
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
   * [1.genotype_dataset_preprocess](#1genotype_dataset_preprocess)
      * [1.1 format genotype dataset](#11-format-genotype-dataset)
      * [1.1 impute genotype](#11-impute-genotype)
      * [1.1 remap genotype dataset](#11-remap-genotype-dataset)
   * [2.RNA_seq_datasets_preprocess](#2rna_seq_datasets_preprocess)
      * [2.1 quality control](#21-quality-control)
      * [2.2 circRNA and mRNA identification](#22-circrna-and-mrna-identification)
   * [3.summary_statistics_genotype](#3summary_statistics_genotype)
      * [3.1 summary statistics](#31-summary-statistics)
      * [3.2 combine genotype datasets](#32-combine-genotype-datasets)
   * [4.summary_statistics_rnaseq](#4summary_statistics_rnaseq)
      * [4.1 combine software results](#41-combine-software-results)
   * [5.matrixeqtl](#5matrixeqtl)
      * [5.1 identify circQTL](#51-identify-circqtl)
      * [5.2 identify eQTL](#52-identify-eqtl)
      * [5.3 enrichment analysis](#53-enrichment-analysis)
   * [6.trans-circQTL](#6trans-circqtl)
      * [6.1 circRNA-mRNA relationship](#61-circrna-mrna-relationship)
   * [7.ABS_events](#7abs_events)
      * [7.1 ABS events construction](#71-abs-events-construction)
      * [7.2 Alu-like elements identification](#72-alu-like-elements-identification)
      * [7.3 intron-ABS relationship](#73-intron-abs-relationship)
      * [7.4 Alu-like ABS relationship](#74-alu-like-abs-relationship)
      * [7.5 Alu-like effect](#75-alu-like-effect)
   * [8.circRNA_functions](#8circrna_functions)
      * [8.0 miRNA binding prediction](#80-mirna-binding-prediction)
      * [8.1 RBP binding prediction](#81-rbp-binding-prediction)
      * [8.2 circ-miRNA network](#82-circ-mirna-network)
      * [8.3 circQTL effect](#83-circqtl-effect)

<!-- Created by https://github.com/ekalinin/github-markdown-toc -->
<!-- Added by: luffy, at: Mon Mar 13 10:28:59 CST 2023 -->

<!--te-->
***

## Overview of pipeline
- 1.genotype_dataset_preprocess

    Raw genotype datasets contain vcf file from BGVD ([Bovine Genome Variant Database](http://animal.omics.pro/code/index.php/BosVar)) 
and two SNP-chip datasets ([GSE95358](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95358) and [GSE100038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100038)).
This preprocess step includes format, imputation, and remap (from btau 6 to btau 9 genome assembly).
- 2.RNA_seq_datasets_preprocess
    
  RNA-seq datasets were retreived from from NCBI SRA database. This preprocess includes quality control, circRNA and mRNA identification.

- 3.summary_statistics_genotype
  
  This step includes summary statistics report of genotype datasets and combination of different genotype datasets from various sources.
- 4.summary_statistics_rnaseq
  
  CircRNA identification was implemented by different software, including CircMarker, CIRI2, and circRNAfinder, and then the detection
results were merged in this step.
- 5.matrixeqtl
  
  CircQTL and eQTL were both identified by MatrixEQTL in R basic environment. Furthermore, the functional genomic region and phenotype enrichment analysis was also carried out. 
- 6.trans-circQTL
  Substantial trans-circQTLs were found in this study, so the preliminary investigation of characteristic of tran-circQTL was involved in this step.
- 7.ABS_events
  
  First, the ABS (alternativ back-splicing) events profile was constructed to investigate the overall patterns and divergence among breeds. 
Then, to explore the association of ABS events (intron length and number/pairing ability of SINE/Alu-like elements) with flanking circQTL mediated by cis-elements (SINE), we extracted the flanking SINE elements closer to
ABS-circRNAs.
- 8.circRNA_functions
  
  CircRNA function mainly include miRNA sponge and RBP interaction at cytoplasm, so the potential miRNA and RBP interact with circRNAs
were predicted, by which we constructed the circRNA-miRNA/RBP interaction networks. To evaluate the distribution and effect of circQTL,
we also investigated the distribution of circQTLs within binding sites and then analyzed the change of binding ability, enrichment degree,
and secondary structure through altering the genotype of circQTLs witinin circRNAs.

## 1.genotype_dataset_preprocess
### 1.1 format genotype dataset
### 1.1 impute genotype
### 1.1 remap genotype dataset

## 2.RNA_seq_datasets_preprocess
### 2.1 quality control
### 2.2 circRNA and mRNA identification

## 3.summary_statistics_genotype
### 3.1 summary statistics
### 3.2 combine genotype datasets

## 4.summary_statistics_rnaseq
### 4.1 combine software results

## 5.matrixeqtl
### 5.1 identify circQTL
### 5.2 identify eQTL
### 5.3 enrichment analysis

## 6.trans-circQTL
### 6.1 circRNA-mRNA relationship

## 7.ABS_events
### 7.1 ABS events construction
### 7.2 Alu-like elements identification
### 7.3 intron-ABS relationship
### 7.4 Alu-like ABS relationship
### 7.5 Alu-like effect

## 8.circRNA_functions
### 8.0 miRNA binding prediction
### 8.1 RBP binding prediction
### 8.2 circ-miRNA network
### 8.3 circQTL effect

