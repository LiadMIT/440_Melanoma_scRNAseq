# ====================================================================
# 6. Differential Gene Expression and Cell Types
# ====================================================================

print("Calculating marker genes for all clusters...")
cluster_markers <- FindAllMarkers(sade_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

print("Extracting the top 10 genes per cluster...")
top10_markers <- cluster_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC)

print(top10_markers, n = 90)

print("Applying new biological identities based on marker genes...")

# 1. Force Seurat to look at the raw numerical clusters

print("Applying biological labels to the downsampled object...")

# 1. Rename the numbers to text, ensuring we pull from and save to sade_seurat3
sade_seurat3 <- RenameIdents(
  sade_seurat,
  "0" = "CD8+ T Cells",
  "1" = "CD4+ T Cells",          
  "2" = "Macrophages",         # Kept separate!
  "3" = "B Cells",
  "4" = "NK & Innate-like T Cells",
  "5" = "pDCs",
  "6" = "cDCs"                 # Kept separate!
)

# 2. Lock it into the metadata (we only need to do this once!)
sade_seurat3$final_cell_type <- Idents(sade_seurat3)

print("Visualizing the mapped biological landscape...")
# 3. Print the UMAP
dim_plot <- DimPlot(sade_seurat3, group.by = "final_cell_type", label = TRUE, label.size = 4, pt.size = 0.1) + 
  ggtitle("UMAP: Initial Cell Types")
print(dim_plot)

# ====================================================================
# 4. Heatmap of initial clustering
# ====================================================================

print("Preparing data for Heatmap...")
# 4. Ensure the active identities are locked to your clean text labels
Idents(sade_seurat3) <- "final_cell_type"

# 5. Scale the marker genes specifically inside sade_seurat3
sade_seurat3 <- ScaleData(sade_seurat3, features = top10_markers$gene)

print("Generating the Single-Cell Clustering Heatmap...")
# 6. Run DoHeatmap strictly on sade_seurat3
cluster_heatmap <- DoHeatmap(
  sade_seurat3, 
  features = top10_markers$gene, 
  size = 4,                     # Size of the cell type labels at the top
  angle = 45                    # Tilts the top labels for readability
) + 
  # Classic, publication-ready color palette
  scale_fill_gradientn(colors = c("dodgerblue4", "white", "firebrick3"), name = "Expression") +
  
  labs(title = "Cell Cluster Heatmap") +
  
  # Clean up the theme so the 60 gene names don't overlap on the Y-axis
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.y = element_text(size = 6, face = "italic", color = "black") 
  )

# Display the plot
print(cluster_heatmap)

# ====================================================================
# 5. Validation of combining Myeloid populations 
# ====================================================================

print("Plotting Shared Pan-Myeloid Markers...")

# The genes that make a myeloid cell a myeloid cell
pan_myeloid_genes <- c(
  "SPI1",     # PU.1 (Master Myeloid Transcription Factor)
  "ITGAX",    # CD11c (Myeloid integrin)
  "ITGAM",    # CD11b (Myeloid integrin)
  "HLA-DRA",  # MHC-II (Shared by both)
  "CD68",     # General Myeloid/Macrophage marker
  "TYROBP"    # DAP12 (Critical myeloid signaling adaptor)
)

DotPlot(sade_seurat3, features = pan_myeloid_genes, group.by = "final_cell_type") +
  scale_color_gradient(low = "lightgrey", high = "darkgreen") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic")) +
  labs(title = "Shared Core Myeloid Phenotype")


print("Re-calculating tree without HLA interference...")

# 1. Get the current highly variable genes
var_genes <- VariableFeatures(sade_seurat3)

# 2. Identify and remove all HLA genes from that list
hla_genes <- grep("^HLA-", var_genes, value = TRUE)
var_genes_no_hla <- setdiff(var_genes, hla_genes)

# 3. NEW STEP: Scale the data specifically for our non-HLA gene list
sade_seurat3_temp <- ScaleData(sade_seurat3, features = var_genes_no_hla, verbose = FALSE)

# 4. Re-run the underlying PCA math using ONLY the non-HLA genes
sade_seurat3_temp <- RunPCA(sade_seurat3_temp, features = var_genes_no_hla, verbose = FALSE)

# 5. Build and plot the new tree
sade_seurat3_temp <- BuildClusterTree(sade_seurat3_temp, dims = 1:20)
PlotClusterTree(sade_seurat3_temp)

