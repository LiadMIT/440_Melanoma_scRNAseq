# =========================================================================
# Fig 3
# =========================================================================

print("Ranking pathways between Responders and Non-Responders...")
# Ranks pathways by how much more active they are in one group vs the other
p_rank <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
print(p_rank)
# Red = Enriched in Responders, Blue = Enriched in Non-Responders

# =========================================================================
# CircosPlots for specific families
# =========================================================================

# 3. Generate the specific CCL plot without the noise
netVisual_chord_gene(chat_R, 
                     signaling = "CCL", 
                     sources.use = all_cells_R, 
                     targets.use = all_cells_R,
                     lab.cex = 1.5,        # Large gene names
                     small.gap = 5,        # High readability
                     #title.name = "Responder CCL Chemokine Signaling"
                     show.legend = FALSE
)

netVisual_chord_gene(chat_NR, 
                     signaling = "CCL", 
                     sources.use = all_cells_NR, 
                     targets.use = all_cells_NR,
                     lab.cex = 1.5,        # Large gene names
                     small.gap = 5,        # High readability
                     show.legend = FALSE
                     
                     #title.name = "Non-responder CCL Chemokine Signaling"
)


# 3. Generate the specific CCL plot without the noise
netVisual_chord_gene(chat_NR, 
                     signaling = "CXCL", 
                     sources.use = all_cells_NR, 
                     targets.use = all_cells_NR,
                     lab.cex = 1.5,        # Large gene names
                     small.gap = 5,        # High readability
                     show.legend = FALSE
                     # title.name = "Non-Responder CXCL Chemokine Signaling"
)


# ====================================================================
# Myeloid patient level chemokine interaction
# ====================================================================

# 1. Defining the complete 8-gene targeted signature
target_genes <- c("CCR5","CCL3","CCR1","CXCL16")

print("1. Fetching data for Myeloid populations...")
expr_data <- FetchData(
  sade_seurat, 
  vars = c(target_genes, "clinical_outcome", "final_cell_type", "true_patient_id")
) %>%
  # Filter to include only relevant Myeloid populations
  filter(grepl("Myeloid APCs", final_cell_type))

# Pivot to long format so we can iterate and plot easily
long_data <- expr_data %>%
  pivot_longer(cols = all_of(target_genes), names_to = "Gene", values_to = "Expression")

print("2. Calculating percentage expressing PER PATIENT (Biological Replicates)...")
patient_gene_props <- long_data %>%
  group_by(final_cell_type, clinical_outcome, true_patient_id, Gene) %>%
  summarise(
    Percent_Expressing = (sum(Expression > 0) / n()) * 100,
    .groups = "drop"
  )

print("3. Calculating Nominal Wilcoxon Stats (Targeted Validation)...")
sig_list <- list()
cell_types <- unique(patient_gene_props$final_cell_type)
genes <- unique(patient_gene_props$Gene)

for (ct in cell_types) {
  for (g in genes) {
    sub_data <- patient_gene_props %>% filter(final_cell_type == ct, Gene == g)
    
    # Run test only if both groups are represented
    if (length(unique(sub_data$clinical_outcome)) == 2) {
      # Suppress warnings for ties, run standard patient-level Wilcoxon
      p_val <- suppressWarnings(wilcox.test(Percent_Expressing ~ clinical_outcome, data = sub_data)$p.value)
    } else {
      p_val <- 1
    }
    
    sig_list[[paste(ct, g)]] <- data.frame(
      final_cell_type = ct, 
      Gene = g, 
      p_val = p_val
    )
  }
}

sig_df <- bind_rows(sig_list)

# Using NOMINAL p-values (no FDR) for this targeted hypothesis test
sig_df$stars <- cut(sig_df$p_val, 
                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                    labels = c("***", "**", "*", ""))

print("4. Formatting for plot and centering stars...")
# Find maximum Y for each facet so the stars hover nicely above the tallest box
y_max <- patient_gene_props %>%
  group_by(final_cell_type, Gene) %>%
  summarise(max_y = max(Percent_Expressing, na.rm = TRUE), .groups = "drop")

sig_df <- merge(sig_df, y_max, by = c("final_cell_type", "Gene"))

# Lock in Gene order so the plot renders exactly how you listed them
patient_gene_props$Gene <- factor(patient_gene_props$Gene, levels = target_genes)
sig_df$Gene <- factor(sig_df$Gene, levels = target_genes)

print("5. Generating the Directed Validation Boxplot...")
full_signature_plot <- ggplot(patient_gene_props, aes(x = Gene, y = Percent_Expressing, fill = clinical_outcome)) +
  # The Boxplots
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black", linewidth = 0.6, 
               position = position_dodge(width = 0.8)) +
  
  # The Patient Replicates (The "Proof")
  geom_point(aes(color = clinical_outcome), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 2, alpha = 0.9, shape = 16, show.legend = FALSE) +
  
  scale_fill_manual(values = c("NR" = "dodgerblue3", "R" = "firebrick3")) +
  scale_color_manual(values = c("NR" = "dodgerblue2", "R" = "firebrick2")) +
  
  # Significance Stars (Nominal)
  geom_text(data = sig_df, aes(x = Gene, y = max_y + (max(patient_gene_props$Percent_Expressing)*0.1), label = stars), 
            inherit.aes = FALSE, 
            size = 6, 
            color = "black", 
            fontface = "bold") +
  
  # Facet by Myeloid Cell Type
  facet_wrap(~ final_cell_type, ncol = 1, scales = "free_y") +
  
  labs(
    title = "Myeloid APC Signaling",
    subtitle = "Patient-Level Wilcoxon Rank Sum Test (Nominal *p < 0.05, **p < 0.01, ***p < 0.001)",
    x = "Gene Name",
    y = "Cells Expressing Gene per Patient (%)",
    fill = "Clinical Outcome"
  ) +
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    legend.position = "top"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3)))

# Print final plot
print(full_signature_plot)

# This will print a clean table of the exact p-values to your R console
knitr::kable(
  sig_df[, c("Gene", "p_val", "stars")], 
  col.names = c("Gene", "Exact P-Value", "Significance Level"),
  format = "markdown"
)
