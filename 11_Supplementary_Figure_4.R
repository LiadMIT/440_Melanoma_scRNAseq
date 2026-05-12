
# ====================================================================
# Unsupervised ligand receptor interactions
# ====================================================================
print("Generating Bubble Plot with reduced font and perfect alignment...")

p_bubble_small <- netVisual_bubble(
  cellchat_merged, 
  comparison = c(1, 2), 
  angle.x = 45, 
  title.name = "Signaling Map: NR vs R",
  font.size = 6              # Lowers the internal CellChat font
) + 
  theme(
    axis.text.x = element_text(
      size = 6,              # Small X-axis labels
      angle = 45, 
      hjust = 1,             # Aligns the end of the word to the dot
      vjust = 1,             # Pulls text up toward the axis
      color = "black"
    ),
    axis.text.y = element_text(size = 5, color = "black"), # Small Y-axis labels
    strip.text.x = element_text(size = 8, face = "bold"),   # Facet headers (NR/R)
    panel.grid.major.x = element_line(color = "grey95", linewidth = 0.3)
  )

print(p_bubble_small)
