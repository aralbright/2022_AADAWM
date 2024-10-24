# 2022_AADAWM

This repository accompanies version 2 of our manuscript, Genome-wide analysis of anterior-posterior mRNA localization in Stentor coeruleus reveals a role for the microtubule cytoskeleton

Manuscript Link: https://www.biorxiv.org/content/10.1101/2023.01.09.523364v2.full
Raw Data Link: https://doi.org/10.6078/D1WT6W

As of Jan. 8, 2023 the file paths in all of these scripts reflect the local environment on Ashley's computer and would have to be changed if you choose to utilize this code. We are continuing to work on making this analysis more reproducible as we continue to update the accompanying preprint-in-progress. 

### 1. Adapter trimming

Code in the file 1_flexbar.sh was adapted to run in a loop from NEB's "Adapter trimming directions for NEBNext Single Cell/Low Input RNA Library Prep Kit for Illumina" 
https://github.com/nebiolabs/nebnext-single-cell-rna-seq

### 2. Kallisto 

The file 2_kallisto_loop.sh contains 2 separate scripts that run locally on the command line. The first section generates the custom reference index used for kallisto, the second is a loop that runs kallisto for all samples. For your reference, all Stentor genome files can be found here: https://stentor.ciliate.org/downloads.php

### 3. Sleuth

3_sleuth.R takes the kallisto results as input, and creates a sleuth dataframe containing normalized estimated transcripts per million (TPM) that we filtered for anterior and posterior samples, log(1+x) transformed, and saves the file 'log_ap_df.csv' for funneling into the Python analysis. 

### 4. Skew Analysis

4_SkewAnalysis.ipynb contains all of the skew analysis. This was used to generate Figure 2 and 3, and contains all tables associated with the manuscript. 

### Extra

The bin folder contains a couple of files required or generated by this analysis, but may not include all necessary files. For example, the CDS fasta required to generate a kallisto index is too big but can be found at the link under the kallisto section above. 
