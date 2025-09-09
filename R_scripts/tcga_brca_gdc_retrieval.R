### --------Important Notes--------
# Note 1:
# Create the project directory in C:\ or D:\ (eg: D:\dir_name)
# That is how GDCdownload() can download the query files. 
# The files could not be downloaded when the dir was nested in subfolders.

# Note 2:
# Update the libraries to avoid errors.


# Note 3:
# Load the exported csv files again to check if all the details are present.


### ---Retrieval of Gene Expression and Subtype Data of TCGA-BRCA from GDC---

# load libraries
library(TCGAbiolinks)
library(SummarizedExperiment)

# get info on TCGA-BRCA project
getProjectSummary('TCGA-BRCA')

# build a query to retrieve gene expression data
query_GeneExp = GDCquery(project = 'TCGA-BRCA',
                         data.category = 'Transcriptome Profiling',
                         data.type = 'Gene Expression Quantification',
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = 'STAR - Counts',
                         sample.type = 'Primary Tumor',
                         access = 'open')

# view the retrieved file details of the query
getResults(query = query_GeneExp)

# download gene expression data
# default dir name is used: "GDCdata"
GDCdownload(query_GeneExp)

# transform the downloaded gene expression files into a summarized count matrix
BRCA_GeneExp = GDCprepare(query_GeneExp, summarizedExperiment = TRUE)

# retrieve gene expression counts 
geneExp_unstranded = assay(BRCA_GeneExp, "unstranded") # for raw counts
# geneExp_tpm_unstrand = assay(BRCA_GeneExp, "tpm_unstrand") # for tpm counts
# geneExp_fpkm_unstrand = assay(BRCA_GeneExp, "fpkm_unstrand") # for fpkm counts


# retrieve metadata of TCGA-BRCA patients (which includes subtype labels)
sample_info = colData(BRCA_GeneExp)

# check column names in sample_info object
# column name containing subtype labels: "paper_BRCA_Subtype_PAM50"
colnames(sample_info) 

# convert sample info object into a data frame
sample_info_df = as.data.frame(sample_info)

# sample_info_df has list columns preventing from saving it as a csv
# check which columns are of list type
list_columns = names(sample_info_df[sapply(sample_info_df, is.list)])

# list columns are "treatments"   "primary_site" "disease_type"

# remove the list columns as they are not needed for subtype classification task
sample_info_df_new = subset(sample_info_df, select = -c(treatments,primary_site,disease_type))

# save the counts and edited sample info as csv files
write.csv(geneExp_unstranded,'geneExp_rawCounts.csv')
write.csv(sample_info_df_new,'sample_info_df.csv')

