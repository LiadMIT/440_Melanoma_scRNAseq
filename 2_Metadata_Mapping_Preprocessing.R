# ====================================================================
# 2. METADATA MAPPING
# ====================================================================
print("Downloading and Loading clinical metadata from NCBI GEO...")
meta_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz"
download.file(meta_url, destfile = "GSE120575_patient_ID_single_cells.txt.gz")

clinical_meta <- read.table(gzfile("GSE120575_patient_ID_single_cells.txt.gz"), 
                            header = TRUE, sep = "\t", skip = 19, nrows = 16291, stringsAsFactors = FALSE)
rownames(clinical_meta) <- clinical_meta[[1]]

print("Injecting true clinical data into Seurat...")
sade_seurat2 <- AddMetaData(sade_seurat2, metadata = clinical_meta)
sade_seurat2$true_patient_id <- sade_seurat2$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.
sade_seurat2$true_clinical_outcome <- sade_seurat2$characteristics..response

# Parse R vs NR and Pre vs Post
sade_seurat2$clinical_outcome <- ifelse(grepl("non", tolower(sade_seurat2$true_clinical_outcome)), "NR", "R")
sade_seurat2$time_point <- ifelse(grepl("Pre", sade_seurat2$true_patient_id), "Pre", "Post")
sade_seurat2$plate_id <- str_extract(colnames(sade_seurat2), "P[0-9]+")

# ====================================================================
# 3. PREPROCESSING, QC, SUBSETTING & REGRESSION
# ====================================================================
print("Identifying Dead/Dying 'Zombie' Cells...")
sade_seurat2[["percent.mt"]] <- PercentageFeatureSet(sade_seurat2, pattern = "^MT-")

print("Scoring the Cell Cycle...")
sade_seurat2 <- CellCycleScoring(
  sade_seurat2, 
  s.features = cc.genes$s.genes, 
  g2m.features = cc.genes$g2m.genes,
  nbin = 5
)

print("Executing Strict Cell QC...")
sade_seurat2 <- subset(sade_seurat2, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5)

print("Isolating Pre-treatment Biopsies...")
sade_seurat2 <- subset(sade_seurat2, time_point == "Pre")

print("Finding Variable Features & Expanding Noise Filter...")
sade_seurat2 <- FindVariableFeatures(sade_seurat2, selection.method = "vst", nfeatures = 3000)
var_genes <- VariableFeatures(sade_seurat2)
noisy_genes <- grep(
  pattern = paste0(
    "^TR[ABGD][VJ]|^IG[HJKL]|",              # TCR and BCR receptors
    "^RPL|^RPS|^MT-|^HSP|",                  # Ribosomal, Mito, Heat Shock
    "MALAT1|NEAT1|^RN7SK|^RN7SL|^SNORD|Y-RNA|", # Hyper-abundant non-coding & structural
    "^FOS|^JUN|CIRBP|RBM3|EGR1|DUSP1|IER3|ATF3|ZFP36|", # FACS / Dissociation / Cold Stress
    "^KRT|",                                 # Keratin / Skin contamination
    "^HBA|^HBB|^HBD|^HBG|",                  # Red Blood Cell Hemoglobin
    "XIST|RPS4Y1|EIF1AY|DDX3Y|KDM5D"         # Sex-linked batch effect genes
  ), 
  x = var_genes, 
  value = TRUE
)

VariableFeatures(sade_seurat2) <- setdiff(var_genes, noisy_genes)

print("The Ultimate Regression...")
sade_seurat2 <- ScaleData(sade_seurat2, 
                          features = VariableFeatures(sade_seurat2), 
                          vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"), 
                          do.center = TRUE, do.scale = TRUE)

sade_seurat <- sade_seurat2