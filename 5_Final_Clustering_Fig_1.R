# ====================================================================
# 6. Combining Myeloid populations
# ====================================================================

# 1. Rename and Merge Identities (Added the missing comma!)
sade_seurat <- RenameIdents(
  sade_seurat3,
  "CD8+ T Cells"= "CD8+ T Cells",
  "CD4+ T Cells"= "CD4+ T Cells",
  "Macrophages" = "Myeloid APCs", 
  "cDCs"        = "Myeloid APCs", 
  "B Cells"     = "B Cells",
  "NK & Innate-like T Cells" = "NK & Innate-like T Cells", # <-- Comma added here
  "pDCs"        = "pDCs"
)

# 2. Save the new merged names permanently to the metadata
sade_seurat$merged_cell_type <- Idents(sade_seurat)

# 3. Calculate markers for the new groupings
print("Calculating marker genes for all clusters...")
cluster_markers <- FindAllMarkers(sade_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 4. Extract and print the top 10 genes
print("Extracting the top 10 genes per cluster...")
top10_markers <- cluster_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC)

print(top10_markers, n = 90)


# Lock it into the metadata so it persists for all future plots
sade_seurat$final_cell_type <- Idents(sade_seurat)

print("Labels successfully applied to sade_seurat2!")

# Save the new names to your metadata for downstream plotting
sade_seurat$final_cell_type <- Idents(sade_seurat)

# Save the new names to the metadata
sade_seurat$cell_type <- Idents(sade_seurat)

print("Visualizing the mapped biological landscape...")
DimPlot(sade_seurat, group.by = "cell_type", label = TRUE, label.size = 5, pt.size = 0.4) + 
  ggtitle("UMAP: Final Cell Types")

print("1. Ensuring the active identities are your new biological labels...")
# This ensures the color bars at the top of the heatmap group by your final cell types
Idents(sade_seurat) <- "final_cell_type"

print("2. Scaling the data for your extracted markers...")
# Seurat's DoHeatmap requires genes to be in the 'scale.data' slot (Z-scored).
# Running ScaleData here ensures these specific 60 genes are perfectly scaled for visualization.
sade_seurat <- ScaleData(sade_seurat, features = top10_markers$gene)


# ====================================================================
# Cell populations Fig 1. A
# ====================================================================

print("Visualizing the mapped biological landscape...")
DimPlot(sade_seurat, group.by = "cell_type", label = TRUE, label.size = 5, pt.size = 0.4) + 
  ggtitle("UMAP: Final Cell Types")

# ====================================================================
# Cell populations Fig 1. B
# ====================================================================

print("1. Fetching data and calculating overall proportions for the bars...")
# Note: Ensure your patient ID column is correctly named below (e.g., "patient_id", "orig.ident")
prop_data <- FetchData(sade_seurat, vars = c("cell_type", "clinical_outcome", "true_patient_id"))

# Calculate the overall bar sizes
prop_summary <- prop_data %>%
  group_by(clinical_outcome, cell_type) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(
    Proportion = Count / sum(Count),
    Percentage = round(Proportion * 100, 1)
  ) %>%
  ungroup()

print("2. Calculating patient-level Wilcoxon statistics for the significance stars...")
# Calculate proportions per individual patient
patient_props <- prop_data %>%
  group_by(clinical_outcome, true_patient_id, cell_type) %>%
  summarise(Cell_Count = n(), .groups = "drop_last") %>%
  mutate(Patient_Proportion = Cell_Count / sum(Cell_Count)) %>%
  ungroup()

cell_types <- na.omit(unique(patient_props$cell_type))
patient_results_list <- list()

for (ct in cell_types) {
  ct_data <- patient_props %>% filter(cell_type == ct)
  # Run the test if data exists in both clinical groups
  if (length(unique(ct_data$clinical_outcome)) == 2) {
    test_res <- wilcox.test(Patient_Proportion ~ clinical_outcome, data = ct_data)
    patient_results_list[[as.character(ct)]] <- data.frame(
      cell_type = ct,
      P_Value = test_res$p.value
    )
  }
}

patient_sig_df <- bind_rows(patient_results_list)

# Apply FDR correction
patient_sig_df$FDR_Adjusted_P <- p.adjust(patient_sig_df$P_Value, method = "fdr")

# Assign stars (non-significant results get left blank so they don't clutter the plot)
patient_sig_df$Significance <- cut(patient_sig_df$FDR_Adjusted_P, 
                                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                                   labels = c("***", "**", "*", "")) 

print("3. Merging data and generating the final plot...")
plot_data <- merge(prop_summary, patient_sig_df, by = "cell_type")

# Create the final text label combining the overall Percentage with the patient-level Significance
# Hides labels for tiny populations (< 2%) to keep the plot clean
plot_data$Label_Text <- ifelse(
  plot_data$Percentage >= 2, 
  paste0(plot_data$Percentage, "% ", plot_data$Significance), 
  ""
)

# Plot it
patient_sig_barplot <- ggplot(plot_data, aes(x = clinical_outcome, y = Proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  
  geom_text(
    aes(label = Label_Text), 
    position = position_fill(vjust = 0.5), 
    size = 4, 
    fontface = "bold", 
    color = "black"
  ) +
  
  scale_y_continuous(labels = percent_format()) +
  
  labs(
    title = "Cell Type Proportions (Pre-R vs Pre-NR)",
    # Updated subtitle to explicitly state the rigorous methodology
    subtitle = "FDR-adjusted Patient-Level Wilcoxon Test (* p<0.05, ** p<0.01, *** p<0.001)",
    x = "Clinical Outcome",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# Display the plot
print(patient_sig_barplot)

print(as.data.frame(publication_table))

# ====================================================================
# Heatmap of new clustering Fig 1c
# ====================================================================

print("3. Generating the Single-Cell Clustering Heatmap...")
cluster_heatmap <- DoHeatmap(
  sade_seurat, 
  features = top10_markers$gene, 
  size = 4,                     # Size of the cell type labels at the top
  angle = 45                    # Tilts the top labels for readability
) + 
  # Classic, publication-ready color palette (Blue = Low, White = Average, Red = High)
  scale_fill_gradientn(colors = c("dodgerblue4", "white", "firebrick3"), name = "Expression") +
  
  labs(title = "Cell Cluster Heatmap") +
  
  # Clean up the theme so the 60 gene names don't overlap on the Y-axis
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.y = element_text(size = 7, face = "italic", color = "black") 
  )

# Display the plot
print(cluster_heatmap)

# Get total cells split by clinical outcome
table(sade_seurat$clinical_outcome)