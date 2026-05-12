# ====================================================================
# QC clustering supplementary Figure 2a
# ====================================================================

print("1. Isolating the Pre-Treatment Baseline (if not already done)...")
# Just in case the object isn't in your environment anymore
seurat_baseline <- sade_seurat

print("2. Generating the 3-Panel Baseline Landscape...")

# Plot A: The Biological Cell Types
p_base_cells <- DimPlot(seurat_baseline, group.by = "cell_type", label = TRUE, label.size = 3, pt.size = 0.2) + 
  ggtitle("Cell Types") + 
  theme(legend.position = "bottom")

# Plot B: Patient Distribution (Checking for batch effects)
p_base_patients <- DimPlot(seurat_baseline, group.by = "true_patient_id", pt.size = 0.2) + 
  ggtitle("Patient Integration") + 
  theme(legend.position = "none") # Hiding legend because there are too many patients

# Plot C: Global Outcome (Responders vs Non-Responders)
p_base_outcome <- DimPlot(seurat_baseline, group.by = "clinical_outcome", cols = c("R" = "firebrick", "NR" = "dodgerblue"), pt.size = 0.2) + 
  ggtitle("Clinical Outcome")

print("3. Stitching the plots together with Patchwork...")
# This syntax uses the patchwork library to arrange the plots in a grid
final_baseline_grid <- (p_base_cells | p_base_patients | p_base_outcome)

# Display the massive combined plot
print(final_baseline_grid)

# ====================================================================
# QC clustering supplementary Figure 2b
# ====================================================================

print("1. Safely determining the new order...")
# Use the correct column name: final_cell_type
current_types <- as.character(unique(sade_seurat$final_cell_type))

# Filter out the target population (matching the exact capitalization we used earlier!)
target_pop <- "NK & Innate-like T Cells"
other_types <- current_types[current_types != target_pop & !is.na(current_types)]

# Recombine them, placing your target population at the very end
new_order <- c(other_types, target_pop)

print("2. Applying the new factor levels...")
# Lock in the new order
sade_seurat$final_cell_type <- factor(as.character(sade_seurat$final_cell_type), levels = new_order)

# Crucial step: Tell Seurat to make this new ordered factor the active identity
Idents(sade_seurat) <- "final_cell_type"

print("3. Generating the aligned Violin Plot...")
# Run the Violin Plot
# Note the '&' instead of '+' so the rotated text applies to EVERY gene facet!
VlnPlot(
  sade_seurat, 
  features = c("CD4", "CD8B", "CD3E", "CD14", "CLEC4C", "CD19"), 
  group.by = "final_cell_type", 
  pt.size = 0
) & 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.title.x = element_blank() # Cleans up the plot by removing the redundant x-axis title
  )

