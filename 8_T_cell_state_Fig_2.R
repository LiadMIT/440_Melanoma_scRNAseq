# ====================================================================
# CD8 T cell analysis
# ====================================================================

# =========================================================================
# PLOT 1: Memory / Exhausted Ratio Boxplot
# =========================================================================
print("1. Defining the exact gene signatures from the paper...")
cd8_g_genes <- list(c("IL7R", "TCF7", "REL", "FOXP1", "FOSL2", "STAT4"))
cd8_b_genes <- list(c("CD38", "HAVCR2", "ENTPD1", "PDCD1", "BATF", "LAG3", "CTLA4", "PTPN6"))

print("2. Scoring every CD8+ T-cell for both states...")
cd8_baseline <- AddModuleScore(cd8_baseline, features = cd8_g_genes, name = "CD8_G_Memory")
cd8_baseline <- AddModuleScore(cd8_baseline, features = cd8_b_genes, name = "CD8_B_Exhausted")

print("3. Classifying each cell based on its dominant state...")
# If the cell's Memory score is higher than its Exhausted score, it becomes CD8_G. 
# Otherwise, it becomes CD8_B.
cd8_baseline$Cell_State <- ifelse(
  cd8_baseline$CD8_G_Memory1 > cd8_baseline$CD8_B_Exhausted1, 
  "CD8_G (Memory)", 
  "CD8_B (Exhausted)"
)



print("Generating the CD8_G / CD8_B Ratio Table...")

# Create a clean table of the counts
state_counts <- table(cd8_baseline$clinical_outcome, cd8_baseline$Cell_State)
state_df <- as.data.frame.matrix(state_counts)

# Calculate the exact ratio (Memory / Exhausted)
state_df$Ratio <- state_df$`CD8_G (Memory)` / state_df$`CD8_B (Exhausted)`
print(state_df)

print("1. Calculating Patient-Level Statistics...")
state_counts <- table(cd8_baseline$clinical_outcome, cd8_baseline$Cell_State)
state_df <- as.data.frame.matrix(state_counts)

# Calculate the exact ratio (Memory / Exhausted)
state_df$Ratio <- state_df$`CD8_G (Memory)` / state_df$`CD8_B (Exhausted)`

patient_ratio_data <- cd8_baseline@meta.data %>%
  group_by(true_patient_id, clinical_outcome, Cell_State) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Cell_State, values_from = cell_count, values_fill = 0) %>%
  mutate(Mem_Ex_Ratio = (`CD8_G (Memory)` + 0.1) / (`CD8_B (Exhausted)` + 0.1)) %>%
  mutate(clinical_outcome = factor(trimws(clinical_outcome), levels = c("NR", "R")))

ratio_test <- wilcox.test(Mem_Ex_Ratio ~ clinical_outcome, data = patient_ratio_data)
p_val <- ratio_test$p.value

ratio_stars <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01   ~ "**",
  p_val < 0.05   ~ "*",
  TRUE ~ "ns"
)

print("Generating the Mechanistic Boxplot...")
max_ratio <- max(patient_ratio_data$Mem_Ex_Ratio)

p_ratio_boxplot <- ggplot(patient_ratio_data, aes(x = clinical_outcome, y = Mem_Ex_Ratio, fill = clinical_outcome)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black", width = 0.5) +
  geom_jitter(width = 0.15, size = 3, shape = 21, color = "black", stroke = 0.5) +
  scale_y_log10() +
  
  # *** COLOR ALIGNMENT: NR is blue, R is red ***
  scale_fill_manual(values = c("NR" = "dodgerblue3", "R" = "firebrick3")) +
  
  annotate("text", x = 1.5, y = max_ratio * 1.5, 
           label = paste0("p = ", signif(p_val, 3), "\n", ratio_stars), 
           size = 5, fontface = "bold") +
  labs(x = "Clinical Outcome", y = "Log10 Ratio CD8_G/CD8_B") +
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black", size = 20),
    axis.title = element_text(face = "bold", size = 20),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    legend.position = "none"
  )
print(p_ratio_boxplot)



# =========================================================================
# PLOT 3: TCF7+ to TCF7- Ratio Boxplot
# =========================================================================
print("1. Calculating TCF7+ to TCF7- Ratio per Patient...")
tcf7_data <- FetchData(cd8_baseline, vars = c("TCF7", "clinical_outcome", "true_patient_id"))

patient_tcf7_ratio <- tcf7_data %>%
  group_by(true_patient_id, clinical_outcome) %>%
  summarise(
    pos_cells = sum(TCF7 > 0),
    neg_cells = sum(TCF7 == 0),
    tcf7_ratio = (pos_cells + 0.1) / (neg_cells + 0.1),
    .groups = "drop"
  ) %>%
  # *** THE FIX: NR level first (Left), R level second (Right) ***
  mutate(clinical_outcome = factor(clinical_outcome, levels = c("NR", "R")))

print("2. Running Patient-Level Statistics...")
ratio_test <- wilcox.test(tcf7_ratio ~ clinical_outcome, data = patient_tcf7_ratio)
p_val <- ratio_test$p.value

ratio_stars <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01   ~ "**",
  p_val < 0.05   ~ "*",
  TRUE ~ "ns"
)

print("3. Generating Ratio Boxplot...")
max_y <- max(patient_tcf7_ratio$tcf7_ratio)

tcf7_ratio_plot <- ggplot(patient_tcf7_ratio, aes(x = clinical_outcome, y = tcf7_ratio, fill = clinical_outcome)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 3, shape = 21, color = "black", stroke = 0.6, alpha = 0.9) +
  
  # *** COLOR ALIGNMENT: NR is blue, R is red ***
  scale_fill_manual(values = c("NR" = "dodgerblue3", "R" = "firebrick3")) +
  
  annotate("text", x = 1.5, y = max_y * 1.05, 
           label = paste0("p = ", signif(p_val, 3), "\n", ratio_stars), 
           fontface = "bold", size = 4.5, vjust = 0) +
  labs(x = "Clinical Outcome", y = "Ratio (TCF7+ / TCF7-)") +
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black", size = 20),
    axis.title = element_text(face = "bold", size = 20),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
print(tcf7_ratio_plot)


# =========================================================================
# PLOT 3: ENTPD1+ to ENTPD1- Ratio Boxplot
# =========================================================================

print("1. Calculating ENTPD1+ to ENTPD1- Ratio per Patient...")

# 1. Fetch raw ENTPD1 counts and patient metadata
entpd1_data <- FetchData(
  cd8_baseline, 
  vars = c("ENTPD1", "clinical_outcome", "true_patient_id")
)

# 2. Calculate the Ratio per patient
patient_entpd1_ratio <- entpd1_data %>%
  group_by(true_patient_id, clinical_outcome) %>%
  summarise(
    pos_cells = sum(ENTPD1 > 0),
    neg_cells = sum(ENTPD1 == 0),
    # Ratio of Pos to Neg (with pseudocount to avoid division by zero)
    entpd1_ratio = (pos_cells + 0.1) / (neg_cells + 0.1),
    .groups = "drop"
  ) %>%
  # *** THE FIX: NR level first (Left), R level second (Right) ***
  mutate(clinical_outcome = factor(clinical_outcome, levels = c("NR", "R")))

print("2. Running Patient-Level Statistics...")

# Wilcoxon test comparing Responders vs Non-Responders
ratio_test <- wilcox.test(entpd1_ratio ~ clinical_outcome, data = patient_entpd1_ratio)
p_val <- ratio_test$p.value

ratio_stars <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE ~ "ns"
)

print("3. Generating Ratio Boxplot (NR on Left) - Space Fixed...")

# Calculate max for star placement
max_y <- max(patient_entpd1_ratio$entpd1_ratio)

entpd1_ratio_plot <- ggplot(patient_entpd1_ratio, aes(x = clinical_outcome, y = entpd1_ratio, fill = clinical_outcome)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 3, shape = 21, color = "black", stroke = 0.6, alpha = 0.9) +
  
  # *** COLOR ALIGNMENT: NR is blue, R is red ***
  scale_fill_manual(values = c("NR" = "dodgerblue3", "R" = "firebrick3")) +
  
  # Place significance stars just above the highest data point
  annotate("text", x = 1.5, y = max_y * 1.05, 
           label = paste0("p = ", signif(p_val, 3), "\n", ratio_stars), 
           fontface = "bold", size = 4.5, vjust = 0) +
  
  labs(
    #title = "CD8+ T-Cell Exhaustion ENTPD1 Ratio",
    #subtitle = "Ratio of (ENTPD1+ / ENTPD1-) Cells per Patient",
    x = "Clinical Outcome",
    y = "Ratio (ENTPD1+ / ENTPD1-)"
  ) +
  
  theme_classic() +
  theme(
    axis.text = element_text(face = "bold", color = "black", size = 20),
    axis.title = element_text(face = "bold", size = 20),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    legend.position = "none"
  ) +
  
  # Provide just enough breathing room (10%) at the top for the stars
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

print(entpd1_ratio_plot)
