# Heldreichia Paper

## Overview
This repository contains R scripts used in the study on Heldreichia species, focusing on SNP data analysis, population structure exploration, and gene enrichment. These scripts form the core analytical framework for the publication. Below is an overview of the key functionalities and contributions of each script:

### Scripts and Their Functions

1. **DAPC Analysis (`dapcnet.R`)**
   - Performs Discriminant Analysis of Principal Components (DAPC) to examine population structure.
   - Identifies SNPs that contribute to group differentiation.
   - Visualizes membership probabilities and SNP contributions using scatter plots, bar plots, and heatmaps.

2. **GO Enrichment Analysis (`go.R`)**
   - Conducts Gene Ontology (GO) enrichment analysis on identified genes.
   - Focuses on Biological Processes (BP) and visualizes results through bar plots, dot plots, and network plots.
   - Outputs enriched pathways and associated genes in structured files for further interpretation.

3. **Random Forest Analysis (`randomf_snps.R`)**
   - Applies random forest classification to identify significant SNPs associated with phenotypes.
   - Includes permutation testing for significance thresholds.
   - Outputs lists of important SNPs, which are further used in other analyses.

4. **PCA Analysis (`snpR_pca.R`)**
   - Conducts Principal Component Analysis (PCA) to visualize population structure.
   - Integrates geographical and phenotypic metadata for enhanced scatter plots.
   - Provides insights into genetic differentiation and clustering.

## File Structure

```
heldreichia_paper/
|— dapcnet.R       # Script for DAPC analysis and SNP heatmap generation
|— go.R            # Script for GO enrichment analysis
|— randomf_snps.R  # Script for random forest analysis of SNPs
|— snpR_pca.R      # Script for PCA analysis of SNP data
```


## Usage Instructions

1. **Prepare Input Data**
   - Ensure VCF files are filtered and ready (e.g., `filtered_snps.vcf.gz`).
   - Provide metadata files for populations and phenotypes (e.g., `sample.list`).

2. **Run Scripts**
   - Execute scripts in the order required for your analysis.
     - For population structure: `dapcnet.R` and `snpR_pca.R`
     - For SNP significance with Random Forest: `randomf_snps.R`
     - For gene enrichment: `go.R`

3. **Output**
   - Outputs include heatmaps, PCA plots, bar plots, and structured data files for significant SNPs and enriched pathways.

