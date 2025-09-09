# Hybrid Autoencoder–Random Forest Approach for Cancer Subtype Classification from RNA-Seq Data

## Problem Statement
Cancer is a heterogeneous disease with several subtypes, each associated with a distinct gene expression profile. Subtype classification is essential for guiding personalized treatment. RNA-Seq provides a comprehensive view of gene expression across the transcriptome, making it a powerful tool for this task. However, its high-dimensionality poses a challenge to building computationally efficient classification models.

## Proposed Solution
To tackle the high-dimensionality of RNA-Seq data, this project explores a hybrid autoencoder-random (AE-RF) framework leveraging differentially expressed genes (DEGs). Autoencoder is a a neural network used for unsupervised learning tasks like dimensionality reduction which is a key part of this project. Here, we used an autoencoder to extract a compact set of relevant features from RNA-Seq data for downstream classification by a random forest. To build a simple and efficient autoencoder, we reduced the initial number of input genes by selecting only DEGs from the whole transcriptome.

The proposed approach is demonstrated on breast cancer subtypes of TCGA-BRCA cohort as a case study, but it can be applied to other cancer types with whole-transcriptome RNA-Seq data.

## Methodology

### Dataset
TCGA-BRCA dataset is used here to implement and evaluate the proposed AE-RF framework. It includes RNA-Seq raw counts and subtype labels of breast cancer patients. The dataset was cleaned, then split into 70% and 30% test set. Raw counts were used for conducting differential gene expression analysis and implementing appropriate normalization pipeline for downstream machine learning models.

#### Code Files and Description
- [`tcga_brca_gdc_retrieval.R`](R_scripts/tcga_brca_gdc_retrieval.R) – Retrieves TCGA-BRCA data from GDC Data Portal using `TCGAbiolinks`.  
- [`data_cleaning_splitting.ipynb`](jupyter_notebook/data_cleaning_splitting.ipynb) – Cleans the TCGA-BRCA dataset and splits it into training and test sets.

