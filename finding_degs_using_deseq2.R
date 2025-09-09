
# load libraries
library(DESeq2)
library(tidyverse)

#-----------------------Functions to streamline workflow------------------------

# Function: get_subtype_degs()
# Algorithm Overview:
#   1. For a given reference subtype, compare it with all other subtypes in a loop.
#   2. Identify upregulated and downregulated genes for each pairwise comparison.
#   3. Find the **common DEGs** that are consistently upregulated or downregulated 
#      across all comparisons with the reference subtype.
# Purpose:
#   - Identify subtype-specific DEGs (both up- and downregulated).
# Parameters:
#   - dds_obj: DESeq2 object containing normalized counts and design.
#   - ref_group: reference subtype to compare against others.
#   - all_groups: vector of all subtypes including the reference.
#   - lfc_thres: log2 fold change threshold for calling DEGs.
#   - adj_pval: adjusted p-value threshold for significance.
get_subtype_degs = function(dds_obj, ref_group, all_groups, lfc_thres,
                                 adj_pval)
{
  
  # define groups to be compared with reference subtype
  comparing_groups = all_groups[all_groups != ref_group]
  
  # initialize variables
  common_upreg_genes = NULL
  common_downreg_genes = NULL
  
  # loop through each comparison subtype
  for (x in comparing_groups) {
    
    # extract DESeq2 results for reference vs current subtype comparison 
    res = results(dds_obj, alpha = adj_pval,
                  contrast = c("breast_subtype", ref_group, x))
    
    # remove NA results and convert to dataframe
    res = na.omit(res) %>% as.data.frame()
    
    # ---upregulated genes ---
    upreg_genes = res[res$log2FoldChange > lfc_thres & res$padj < adj_pval, ]
    upreg_gene_names = row.names(upreg_genes)
    
    # keep only genes consistently upregulated across all comparisons
    if (is.null(common_upreg_genes)) {
      common_upreg_genes = upreg_gene_names
    } else {
      common_upreg_genes = intersect(common_upreg_genes, upreg_gene_names)
    }
   
    
    # ---downregulated genes ---
    downreg_genes = res[res$log2FoldChange < -lfc_thres & res$padj < adj_pval, ]
    downreg_gene_names = row.names(downreg_genes)
    
    # keep only genes consistently downregulated across all comparisons
    if (is.null(common_downreg_genes)) {
      common_downreg_genes = downreg_gene_names
    } else {
      common_downreg_genes = intersect(common_downreg_genes, downreg_gene_names)
    }
  }
  
  # convert to dataframes for both up- and downregulated genes
  df_upreg_genes = data.frame(gene_id = common_upreg_genes, type = "upregulated")
  df_downreg_genes = data.frame(gene_id = common_downreg_genes, type = "downregulated")
  
  # combine both up- and downregulated genes
  common_degs = rbind(df_upreg_genes, df_downreg_genes)
  
  # print messages summarizing the results
  message(paste("Common significant DEGs for LFC >", lfc_thres, ":", nrow(df_upreg_genes)))
  message(paste("Common significant DEGs for LFC <", lfc_thres, ":", nrow(df_downreg_genes)))
  message(paste("Common significant DEGs for LFC > and <", lfc_thres, ":", nrow(common_degs)))
  message(paste("Retrieval of DEGs specific to", ref_group, "done!"))
  
  return(common_degs)
}


#------------------------------1. Import files----------------------------------

# load sample metadata of the training set
col_data = read.csv("C:\\Users\\User\\Documents\\brca_subtype_ae\\data\\train_test_sets\\X_raw_counts_train.csv", header = TRUE, row.names = 1,
                    check.names = FALSE)

# load raw counts data of the training set
# Note: genes should be rows, samples should be columns
# The below file is not, so it is transposed
counts_data = t(read.csv("C:\\Users\\User\\Documents\\brca_subtype_ae\\data\\train_test_sets\\y_subtype_train.csv", header = TRUE, row.names = 1,
                         check.names = FALSE))


#--------------------2. Preparation prior to DESeq analysis---------------------

# check dimensions
dim(col_data) # (721, 1)
dim(counts_data) # (60660, 721)

# select subtype column in sample_info and omit the rest, if they are present
# drop = FALSE ensures it is a dataframe
col_data = col_data[, "paper_BRCA_Subtype_PAM50", drop = FALSE]

# rename the column for simplicity
colnames(col_data) = c("breast_subtype")

# convert the subtype feature to factor
col_data$breast_subtype = as.factor(col_data$breast_subtype)

# check whether all row names of col_data are in colnames of counts_data
# output should be true
all(row.names(col_data) %in% colnames(counts_data))

# check order row names of col_data match with colnames of counts_data
# output should be true
all(row.names(col_data) == colnames(counts_data))

# construct a DESeq dataset object
dds = DESeqDataSetFromMatrix(countData = counts_data,
                             colData = col_data,
                             design = ~breast_subtype)

# pre-filtering of low count genes to remove noise and speed up computation (based on DESeq2 vignette)
# keep only rows that have a count of at least 10 for a minimal number of samples.
# minimal number of samples is the smallest group size (Her2: 57)
smallestGroupsize = min(table(col_data$breast_subtype)) 
keep = rowSums(counts(dds) >= 10) >= smallestGroupsize
dds = dds[keep,]

# dimension after pre-filtering
dim(dds) # (26640, 721)

# extract the list of all groups/subtypes
subtype_list = levels(dds$breast_subtype)

#-------------------------------3. DESeq2 analysis-------------------------------

# run DESeq2 with default parameters
deseq_obj = DESeq(dds)

# plot dispersion estimate plot to check RNA-seq data quality
# expected pattern: 
 # Decreasing dispersion with increasing mean expression curve
 # Smooth and central red line (fitted dispersion trend)
 # Most points cluster around the red line
plotDispEsts(deseq_obj, ylim = c(1e-3, 1e3), legend = FALSE)


# extract consistent DEGs (up- and downregulated) for each subtype with:
#   log2 fold change > and < 1
#   adjusted p-value < 0.05
# use get_subtype_degs() that is defined at the top
basal_degs = get_subtype_degs(dds_obj = deseq_obj, ref_group = "Basal", 
                              all_groups = subtype_list, lfc_thres = 1, 
                              adj_pval = 0.05)
her2_degs = get_subtype_degs(dds_obj = deseq_obj, ref_group = "Her2", 
                             all_groups = subtype_list, lfc_thres = 1, 
                             adj_pval = 0.05)
lumA_degs = get_subtype_degs(dds_obj = deseq_obj, ref_group = "LumA", 
                             all_groups = subtype_list, lfc_thres = 1, 
                             adj_pval = 0.05)
lumB_degs = get_subtype_degs(dds_obj = deseq_obj, ref_group = "LumB", 
                             all_groups = subtype_list, lfc_thres = 1, 
                             adj_pval = 0.05)

# combine the DEGs of all subtypes
all_subtype_degs = Reduce(union,list(basal_degs$gene_id,her2_degs$gene_id,lumA_degs$gene_id,
                                lumB_degs$gene_id)) %>% as.data.frame(col.names = "gene_id")

write.csv(all_subtype_degs, "degs_training_set.csv")

