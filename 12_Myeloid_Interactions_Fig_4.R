# ====================================================================
# Fig 4
# ====================================================================

# ====================================================================
# Incoming CD8+ T cell ligand receptor interaction
# ====================================================================

# 2. Clear the plot device for fresh margins
par(mar = c(1, 1, 1, 1), xpd = TRUE)


# 1. Generate the base plot
p_bubble <- netVisual_bubble(
  cellchat_merged, 
  targets.use = "CD8+ T Cells", 
  comparison = c(1, 2), 
  angle.x = 45, 
  title.name = "Incoming Signaling: NR vs R CD8+ T cells"
)

# 2. Apply the 'Force & Align' fix AND drop the NA columns
p_bubble_fixed <- p_bubble + 
  # na.translate = FALSE tells ggplot to mathematically delete the empty NA columns!
  scale_x_discrete(guide = guide_axis(n.dodge = 1), na.translate = FALSE) + 
  scale_y_discrete(na.translate = FALSE) + # Added here just in case an NA row pops up too
  theme(
    # Small font (6-7) is required to fit all 10+ columns in one view
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    # Bold the NR/R headers
    strip.text.x = element_text(size = 9, face = "bold"),
    # Add vertical grid lines to guide the eye from the bubble to the label
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
    # Ensure the plot doesn't cut off the bottom labels
    plot.margin = margin(b = 20, r = 5, l = 5, t = 5)
  )

print(p_bubble_fixed)


# ====================================================================
# Myeloid patient level immunosuppression table
# ====================================================================

# 1. Defining the complete 8-gene targeted signature
target_genes <- c("CXCL16", "LGALS9", "SIGLEC1", "MIF","SPP1")

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
    title = "Myeloid APC Immunosuppression",
    subtitle = "Patient-Level Wilcoxon Rank Sum Test (Nominal *p < 0.05, **p < 0.01, ***p < 0.001)",
    x = "Gene Name",
    y = "% Myeloid APCs Expressing Gene per Patient",
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

# Clean up the display and sort by significance
p_val_summary <- sig_df %>%
  select(final_cell_type, Gene, p_val, stars) %>%
  arrange(p_val)

print("--- Myeloid APC Statistical Summary ---")
print(p_val_summary)

# ====================================================================
# Myeloid patient level chemokine interaction
# ====================================================================

print("1. Fetching data and filtering for Myeloid APCs...")
target_genes <- c("LGALS9", "SIGLEC1", "CXCL16")

coexp_data <- FetchData(
  sade_seurat, 
  vars = c(target_genes, "clinical_outcome", "final_cell_type", "true_patient_id")
) %>%
  filter(grepl("Myeloid APCs", final_cell_type))

print("2. Determining Triple Co-expression...")
coexp_data$Coexpressed <- ifelse(
  coexp_data$LGALS9 > 0 & coexp_data$SIGLEC1 > 0 & coexp_data$CXCL16 > 0,
  "Yes", "No"
)

print("3. Calculating percentages PER PATIENT...")
patient_coexp <- coexp_data %>%
  group_by(clinical_outcome, true_patient_id) %>%
  summarise(
    Total_Cells = n(),
    Co_expressing_Cells = sum(Coexpressed == "Yes"),
    Patient_Percent = (Co_expressing_Cells / Total_Cells) * 100,
    .groups = 'drop'
  ) %>%
  # *** THE FLIP: NR level first (Left), R level second (Right) ***
  mutate(clinical_outcome = factor(clinical_outcome, levels = c("NR", "R")))

print("4. Calculating Statistics (NR vs R)...")
p_val <- wilcox.test(Patient_Percent ~ clinical_outcome, data = patient_coexp)$p.value

stars <- cut(p_val, 
             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
             labels = c("***", "**", "*", "ns"))

print("5. Generating NR/R Comparison Plot...")
max_y <- max(patient_coexp$Patient_Percent, na.rm = TRUE)

coexp_r_nr_plot <- ggplot(patient_coexp, aes(x = clinical_outcome, y = Patient_Percent, fill = clinical_outcome)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 3, shape = 21, color = "black", stroke = 0.6) +
  
  # *** COLOR ASSIGNMENT: Blue for NR (Left), Red for R (Right) ***
  scale_fill_manual(values = c("NR" = "dodgerblue3", "R" = "firebrick3")) +
  
  # Significance Annotation sits centered at x = 1.5
  annotate("text", x = 1.5, y = max_y * 1.1, 
           label = paste0("p = ", signif(p_val, 3), "\n", stars),
           fontface = "bold", size = 5) +
  
  labs(
    title = "Triple Positive Myeloid APCs (LGALS9+SIGLEC1+CXCL16+)",
    subtitle = "Patient-Level Comparison: Non-Responders vs Responders",
    x = "Clinical Outcome",
    y = "% of Myeloid APCs per Patient"
  ) +
  
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))

print(coexp_r_nr_plot)

