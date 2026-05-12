
# =========================================================================
# Unsupervised CircosPlots
# =========================================================================

print("Calculating differential interaction matrix...")
# We subtract NR from R to see what is 'gained' or 'lost'
# Since NR has way more interactions, this plot will be dominated by NR signals
num_link <- cellchat_merged@net$count[,,2] - cellchat_merged@net$count[,,1]

# 1. Set up the plotting area
par(mfrow = c(1,2), xpd = TRUE)

# 1. Get all unique cell types in your NR object
all_cells_NR <- levels(chat_NR@idents)

all_cells_R <- levels(chat_R@idents)

# 4. Clear device and set margins
dev.off() 
par(mar = c(1, 1, 1, 1), xpd = TRUE)

# 5. Generate the Cleaned Gene-Level Chord Plot
netVisual_chord_gene(chat_NR, 
                     sources.use = all_cells_NR,
                     targets.use = all_cells_NR,
                     slot.name = "netP", 
                     lab.cex = 0.6,           # Slightly larger now that space is freed up
                     small.gap = 1.5, 
                     big.gap = 20,            # Larger gaps between the real cell types
                     title.name = "Non-Responder Baseline Interactions",
                     show.legend = TRUE)

# 5. Generate the Cleaned Gene-Level Chord Plot
netVisual_chord_gene(chat_R, 
                     sources.use = all_cells_R,
                     targets.use = all_cells_R,
                     slot.name = "netP", 
                     lab.cex = 0.6,           # Slightly larger now that space is freed up
                     small.gap = 1.5, 
                     big.gap = 20,            # Larger gaps between the real cell types
                     title.name = "Responder Baseline Interactions",
                     show.legend = TRUE)
