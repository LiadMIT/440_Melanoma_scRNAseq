# ====================================================================
# Cell Population Table
# ====================================================================

# Create the table: rows = cell types, columns = patients
patient_table <- table(Idents(sade_seurat), sade_seurat$true_patient_id)

# View the organized table
print(patient_table)

# Create the table
patient_table <- as.data.frame.matrix(table(Idents(sade_seurat), sade_seurat$true_patient_id))

# Save to file
write.csv(patient_table, "Sade_Feldman_Patient_Cell_Counts.csv")

# ====================================================================
# Patient summary table
# ====================================================================

# 1. Extract the unique patient IDs and their outcomes
patient_info <- unique(seurat_baseline@meta.data[, c("true_patient_id", "clinical_outcome")])

# 2. Calculate the total cells for each of those patients
total_counts <- table(seurat_baseline$true_patient_id)

# 3. Combine them into a single clean data frame
baseline_summary <- data.frame(
  Patient_ID = patient_info$true_patient_id,
  Outcome = patient_info$clinical_outcome,
  Total_Cells = as.numeric(total_counts[patient_info$true_patient_id])
)

# 4. Sort by Outcome so NR and R are grouped
baseline_summary <- baseline_summary[order(baseline_summary$Outcome), ]

# 5. Print the result
print(baseline_summary)

# 6. Calculate the Grand Totals for your report
print(paste("Total NR Cells:", sum(baseline_summary$Total_Cells[baseline_summary$Outcome == "NR"])))
print(paste("Total R Cells:", sum(baseline_summary$Total_Cells[baseline_summary$Outcome == "R"])))
