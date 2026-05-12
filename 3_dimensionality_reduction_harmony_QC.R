# ====================================================================
# 4. DIMENSIONALITY REDUCTION, HARMONY & CLUSTERING
# ====================================================================
print("Running PCA, Harmony, and UMAP...")
sade_seurat <- RunPCA(sade_seurat, features = VariableFeatures(sade_seurat), npcs = 50, verbose = FALSE)

# Note: Changed group.by.vars to "true_patient_id" to match the metadata explicitly
sade_seurat <- RunHarmony(sade_seurat, 
                          group.by.vars = "true_patient_id", 
                          theta = c(3), 
                          dims.use = 1:50, 
                          max_iter = 20)
ElbowPlot(sade_seurat)

sade_seurat <- RunUMAP(sade_seurat, reduction = "harmony", dims = 1:20, n.neighbors = 50, min.dist = 0.5)
sade_seurat <- FindNeighbors(sade_seurat, reduction = "harmony", dims = 1:20)
sade_seurat <- FindClusters(sade_seurat, resolution = 0.25)

sade_seurat$cell_type <- Idents(sade_seurat)

print("Dieting the Object & Clearing Memory...")
sade_seurat_slim <- DietSeurat(sade_seurat, assays = "RNA", layers = c("counts", "data"), dimreducs = c("pca", "harmony", "umap"))
sade_seurat <- sade_seurat_slim
rm(sade_seurat_slim, sade_seurat2); gc(); gc()

DimPlot(sade_seurat, group.by = "seurat_clusters", label = TRUE, label.size = 4, pt.size = 0.8) + 
  ggtitle("UMAP: Numeric Clusters")
# 1. Open the PDF device
# Using a square aspect ratio (7x7 or 8x8) usually looks best for UMAPs

# 3. Close the device
dev.off()

# ====================================================================
# 5. Silhouette scores
# ====================================================================

# 1. Extract PCA coordinates (usually top 20 or 30 PCs are used for clustering)
# We use the same number of PCs you used for FindNeighbors
pc_dims <- 1:20
data_dims <- Embeddings(sade_seurat, reduction = "pca")[, pc_dims]

# 2. Extract cluster assignments
clusters <- Idents(sade_seurat)

# 3. Calculate the distance matrix and silhouette scores
# Warning: For very large datasets (>50k cells), this distance matrix can eat RAM.
dist_matrix <- dist(data_dims)
sil_score <- silhouette(as.numeric(clusters), dist_matrix)

# 4. Get the Average Silhouette Width
avg_sil_width <- mean(sil_score[, 3])
print(paste("Average Silhouette Score:", round(avg_sil_width, 4)))

# Summary by cluster
summary(sil_score)

par(mar = c(5, 15, 5, 2) + 0.1)

plot(sil_score, 
     col = 1:length(unique(clusters)), 
     border = NA, 
     main = "Silhouette Analysis: Immune Cluster Separation Quality",
     cex.names = 0.8,   # Slightly smaller text to prevent overlapping
     nmax.lab = 50,     # Ensures labels are shown even for many clusters
     xlab = "Silhouette Width (s_i)")
