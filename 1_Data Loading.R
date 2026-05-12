library(vroom)
library(stringr)
library(Seurat)
library(patchwork)
library(CellChat)
library(harmony)
library(data.table)
library(ggplot2)
library(dplyr)
library(scales) 
library(cluster)
install.packages("ape")
library(ape)

options(Seurat.object.assay.version = "v3")

Sys.setenv("VROOM_CONNECTION_SIZE" = 500000000)
rm(list = ls()); gc() 

# ====================================================================
# 1. DATA LOADING & FORCING SPARSITY
# ====================================================================
print("Loading massive matrix safely with vroom...")
counts_df <- vroom("GSE120575_counts.txt.gz", col_names = TRUE, show_col_types = FALSE)
genes_clean <- make.unique(gsub("_", "-", counts_df[[1]]))

print("Converting to standard matrix (Bypassing the factor trap)...")
raw_matrix <- as.matrix(counts_df[, -1])
suppressWarnings(class(raw_matrix) <- "numeric")

print("Executing the Sparsity Fix...")
raw_matrix[is.na(raw_matrix)] <- 0
raw_matrix[raw_matrix < 0.01] <- 0 

print("Building the Sparse Matrix (Compression will actually work now!)...")
counts_sparse <- as(raw_matrix, "dgCMatrix")

rownames(counts_sparse) <- genes_clean
colnames(counts_sparse) <- colnames(counts_df)[-1]

rm(counts_df, raw_matrix); gc() 

print("Filtering NA genes and building QC-filtered Object...")
valid_rows <- which(!is.na(rownames(counts_sparse)) & rownames(counts_sparse) != "")
counts_sparse <- counts_sparse[valid_rows, ]

sade_seurat <- CreateSeuratObject(counts = counts_sparse, project = "SadeFeldman", min.cells = 3, min.features = 200) 

sade_seurat2 <- sade_seurat 