# Hybrid Autoencoder–Random Forest Approach for Cancer Subtype Classification from RNA-Seq Data

## Problem Statement
Cancer is a heterogeneous disease with several subtypes, each associated with a distinct gene expression profile. Subtype classification is essential for guiding personalized treatment. RNA-Seq provides a comprehensive view of gene expression across the transcriptome, making it a powerful tool for this task. However, its high-dimensionality poses a challenge to building computationally efficient classification models.

## Proposed Solution
To tackle the high-dimensionality of RNA-Seq data, this project explores a hybrid autoencoder-random (AE-RF) framework leveraging differentially expressed genes (DEGs). Autoencoder is a a neural network used for unsupervised learning tasks like dimensionality reduction which is a key part of this project. Here, we used an autoencoder to extract a compact set of relevant features from RNA-Seq data for downstream classification by a random forest. To build a simple and efficient autoencoder, we reduced the initial number of input genes by selecting only DEGs from the whole transcriptome.

The proposed approach is demonstrated on TCGA-BRCA breast cancer dataset as an example, but it can be applied to other cancer types with whole-transcriptome RNA-Seq data.

## Methodology

### Dataset
TCGA-BRCA dataset was retrieved from GDC Data Portal using `TCGAbiolinks` package, and includes whole transcriptome RNA-Seq raw counts and subtype labels of breast cancer patients. The dataset was cleaned, then split into 70% training and 30% test sets. Raw counts were opted instead of TPM- or FPKM-normalized counts for conducting differential gene expression analysis and implementing an appropriate normalization pipeline for downstream machine learning models.

### Differentially Expressed Gene (DEG) Identification
DEGs were identified from the raw counts of the training set using the `DESeq2` package. Breast cancer has four major molecular subtypes: Luminal A, Luminal B, HER2-enriched and Basal-like. Here are the steps performed for finding DEGs of the subtypes after running DESeq2 differential gene expression analysis:

1. For each reference subtype, perform pairwise comparison with each of the other subtypes.
2. Identify upregulated and downregulated genes for each pairwise comparison based on:
   * Log2 fold change > 1 and < −1 for biological significance
   * Adjusted p-value < 0.05 for statistical significance.
3. Retain only genes that are consistently up- or downregulated across all comparisons with the reference subtype.
4. Combine DEGs from all subtypes to obtain a final set of genes for downstream analysis.

### AE-RF Model Architecture
![Model Architecture](model_architecture.png)

## Code Files
- [`tcga_brca_gdc_retrieval.R`](R_scripts/tcga_brca_gdc_retrieval.R) – Retrieves TCGA-BRCA data from GDC Data Portal using `TCGAbiolinks`.  
- [`data_cleaning_splitting.ipynb`](jupyter_notebooks/data_cleaning_splitting.ipynb) – Cleans the TCGA-BRCA dataset and splits it into training and test sets.
- [`finding_degs_using_deseq2.R`](R_scripts/finding_degs_using_deseq2.R) – Finds DEGs from the whole transcriptome using `DESeq2`. 
- [`AE_RF_classification.ipynb`](jupyter_notebooks/AE_RF_classification.ipynb) – Normalizes RNA-Seq raw counts, trains AE-RF, evaluates on test set, and analyzes model performance

## Future Improvements
To be updated soon...

## References
To be updated soon...
