# RBP Modules: Integrative Analysis of RBP Similarities and Functional Annotation

This repository contains a set of scripts designed for (1) an integrative assessment of RNA-binding protein (RBP) regulatory interactions, and (2) functional annotation of RBPs using BioID data.

## System Requirements

### Software Dependencies and Versions
- R version 4.2.2
- R packages:
  - tidyr (v1.3.0)
  - gage (v2.48.0)
  - AnnotationDbi (v1.60.2)
  - org.Hs.eg.db (v3.16.0)
  - fgsea (v1.24.0)
  - gplots (v3.1.3.1)
  - colorspace (v2.1-0)
  - extrafont (v0.19)
  - extrafontdb (v1.0)
  - tidyverse (v2.0.0)
  - metap (v1.1)
  - stringr (v1.5.0)

For the installation of individual package dependencies, please follow the official installation guides provided by each package's documentation. Installation time might vary depending on your system.

### Operating Systems
- Tested on Ubuntu 20.04 LTS and Windows 10

### Non-Standard Hardware
No non-standard hardware is required.

## Installation Guide

1. Clone the repository:
```
git clone https://github.com/goodarzilab/RBP_modules.git
```
2. Navigate to the repository folder:
```
cd RBP_modules
```

## Demo

### Instructions

1. Unzip the archive with the example files into the data folder
```
unzip data.zip -d data
```
2. Run the dataset integration:
```
Rscript bin/Datasets_integration.R BioID=data/signed_log_pv_BioID.tsv eCLIP=data/signed_log_pv_eCLIP.tsv PerturbSeq=data/signed_log_pv_Perturbseq.tsv
```
3. Run the functional annotation:
```
Rscript bin/BioID_annotation.R nperm=1000 nproc=10 input=data/signed_log_pv_BioID.tsv
```

### Expected Output
- The scripts will generate PDF files `Logitp_heatmap.pdf` and `BioID_based_annotation_BP.pdf` containing the respective heatmaps.

### Expected Run Time for Demo
- The expected run time for the demo depends on the capabilities of your computer and the number of processes used (nproc parameter).

## Additional Information

### License
This software is distributed under the GNU GENERAL PUBLIC LICENSE. For more information, see the LICENSE file in the repository.


## Analysis Workflow Description

### 1. Dataset integration for estimating functional RBP similarities
The data integration procedure is handled by **Datasets_integration.R**. The script normalizes signed p-values derived from the analysis of the BioID, eCLIP, and Perturb-seq experiments using Z-scoring followed by cosine distances estimation. Next, the calculation of the empirical p-values using sorted RBP-RBP distances and p-values aggregation are performed to produce the integrated pairwise distances.

In more detail, the functional similarity of RBPs is estimated by joint analysis of eCLIP, BioID, and Perturb-seq data. First, using signed p-values, Z-scores are calculated for every row across 50 and 68 RBPs separately for BioID (for 2936 preys) and Perturb-seq (for 10 latent vectors), respectively. For eCLIP, the Z-scoring step is skipped. 
Next, cosine distance is computed for all RBP pairs (1225 for BioID, 7140 for eCLIP and 2278 for Perturb-seq) followed by ranking and calculation of empirical p-values defined as a fraction of RBP pairs with the cosine distance less than the score of the tested pair. The resulting empirical p-values are aggregated with the logitp function of metap. Raw (non-aggregated) p-values are used for the RBP pairs assessed in a single type of experiment. Finally, global integration heatmap is generated using cosine distance (1 - cosine similarity) and Ward’s clusterization.

The three source files containing signed p-values obtained from the experimental data (signed_log_pv_BioID.tsv, signed_log_pv_eCLIP.tsv, signed_log_pv_Perturbseq.tsv) are available in the 'data' subfolder.
The file contents are as follows:
1. **signed_log_pv_BioID.tsv** contains -log<sub>10</sub>(p-values) * sing(log<sub>2</sub>FC) estimates for each bait-prey pair (50 RBPs x 2936 preys)
2. **signed_log_pv_eCLIP.tsv** is gene-level (120 RBPs x 22471 genes)
3. **signed_log_pv_Perturbseq.tsv** includes the gene expression alteration estimates obtained from the autoencoder (first 10 latent vectors, 68 RBPs x 10 LVs).

### 2. Functional annotation of RBPs using BioID data
**BioID_annotation.R** annotates 50 RBPs with GO terms based on their protein neighborhoods identified by BioID using signed_log_pv_BioID.tsv source file. For this, lists of protein partners for each RBP are ranked by z-scored signed p-values obtained from BioID and tested for the GO terms enrichment using gene set enrichment analysis (gsea). The resulting GO pathways are filtered by NES > 2 for at least one RBP followed by bi-clustering and heatmap visualization.

In detail, the script annotates the RBPs based upon preys identified in BioID experiments starting from the signed p-values estimated as -log<sub>10</sub>(p-values)\*sign(log<sub>2</sub>FC) for every bait-prey pair separately. For each prey, signed p-values are converted to Z-scores by estimating the mean and sd across baits. Z-scores rank the preys and gsea is applied using three GO terms annotation sets (biological processes (BP), molecular functions (MF), and cell components (CC), each taken separately). Lists of 2865 preys Entrez ids are used in gene set enrichment analysis for each RBP of the whole set of fifty. In the end, Ward’s clusterization with cosine distance is used to generate the GO BP annotation heatmap.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10498278.svg)](https://doi.org/10.5281/zenodo.10498278)
