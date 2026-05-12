# ====================================================================
# CellChat Ligand Receptor Mapping
# ====================================================================

sade_lean <- subset(seurat_baseline)
#rm(sade_seurat, seurat_baseline); gc() # Wipe massive objects before running CellChat

seurat_R <- subset(sade_lean, clinical_outcome == "R")
seurat_NR <- subset(sade_lean, clinical_outcome == "NR")

print("Extracting data manually to bypass Seurat v5 layer restrictions...")
data_R <- GetAssayData(seurat_R, layer = "data"); meta_R <- seurat_R@meta.data
data_NR <- GetAssayData(seurat_NR, layer = "data"); meta_NR <- seurat_NR@meta.data

print("Initializing independent CellChat objects...")
chat_R <- createCellChat(object = data_R, meta = meta_R, group.by = "cell_type")
chat_NR <- createCellChat(object = data_NR, meta = meta_NR, group.by = "cell_type")

CellChatDB <- CellChatDB.human
chat_R@DB <- CellChatDB
chat_NR@DB <- CellChatDB

print("Running CellChat Core Pipeline for Responders (R)...")
chat_R <- subsetData(chat_R) |> identifyOverExpressedGenes() |> identifyOverExpressedInteractions() |> 
  computeCommunProb(type = "triMean") |> filterCommunication(min.cells = 5) |> 
  computeCommunProbPathway() |> aggregateNet()

print("Running CellChat Core Pipeline for Non-Responders (NR)...")
chat_NR <- subsetData(chat_NR) |> identifyOverExpressedGenes() |> identifyOverExpressedInteractions() |> 
  computeCommunProb(type = "triMean") |> filterCommunication(min.cells = 5) |> 
  computeCommunProbPathway() |> aggregateNet()

print("Merging Baseline Networks for Direct Comparison...")
object.list <- list(NR = chat_NR, R = chat_R)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

print("1. Standardizing nodes across both networks (The 'Lift' Fix)...")
# Get the master list of all cell types from your parent object
unified_levels <- levels(sade_lean$cell_type)

# Force both networks to adopt this exact structure
chat_R <- liftCellChat(chat_R, unified_levels)
chat_NR <- liftCellChat(chat_NR, unified_levels)

print("2. Re-merging the perfectly aligned networks...")
object.list <- list(NR = chat_NR, R = chat_R)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))


# ====================================================================
# 7. CellChat overview
# ====================================================================

print("1. The Barplot (Total Communication Strength)")
# This simple bar chart shows if R or NR has more overall communication
p_bar <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "count")
print(p_bar)

print("2. The Differential Heatmap")
# This is much cleaner than the circle plot. 
# Red squares = pathway turned ON in Responders. 
# Blue squares = pathway turned OFF in Responders.
p_diff_heat <- netVisual_heatmap(cellchat_merged, measure = "weight")
print(p_diff_heat)

print("Ranking pathways between Responders and Non-Responders...")
# Ranks pathways by how much more active they are in one group vs the other
p_rank <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
print(p_rank)
# Red = Enriched in Responders, Blue = Enriched in Non-Responders
