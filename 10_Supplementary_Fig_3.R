# ====================================================================
# UMAPs for each patient
# ====================================================================

# Create a new metadata column: "Outcome_PatientID"
seurat_baseline$patient_outcome_label <- paste0(seurat_baseline$clinical_outcome, "_", seurat_baseline$true_patient_id)

# Sort them so all NR patients appear first, then all R patients
# This makes the "Response vs Non-Response" comparison much easier to see
patient_order <- sort(unique(seurat_baseline$patient_outcome_label))
seurat_baseline$patient_outcome_label <- factor(seurat_baseline$patient_outcome_label, levels = patient_order)

print("Generating UMAPs split by patient and outcome...")

# We use ncol=4 or 5 depending on how many patients you have (approx 19)
DimPlot(seurat_baseline, 
        reduction = "umap", 
        group.by = "cell_type", 
        split.by = "patient_outcome_label", 
        ncol = 5, 
        pt.size = 0.4, 
        label = FALSE) +
  theme(strip.text = element_text(size = 8, face = "bold")) + # Makes headers readable
  ggtitle("Patient-Specific Landscapes: Non-Responders vs Responders")

