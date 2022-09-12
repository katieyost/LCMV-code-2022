# Load packages

library(Seurat)
library(SingleR)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 2400 * 1024^2)
options(future.rng.onMisuse="ignore")
library(stringi)

#-----------------
# Functions
#-----------------

CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), tolower(substring(c, 2)),
      sep="", collapse=" ")
}

#-----------------
# Cluster cells
#-----------------

# Read in seurat object

merge_obj <- readRDS("data/seurat_obj.rds")

# Filter doublets, low quality cells and non-T cells

merge_obj <- subset(merge_obj, subset = nFeature_RNA > 200 & percent.mt < 5)
merge_obj <- subset(merge_obj, subset = DF.classifications == "Singlet")
merge_obj <- subset(merge_obj, subset = SingleR_labels == "T cells" | SingleR_labels == "NK cells")

# Add cell cycle score

cc.genes <- readLines(con = "data/regev_lab_cell_cycle_genes.txt")
cc.genes <- sapply(cc.genes, CapStr, USE.NAMES = FALSE)
s.genes <- cc.genes[1:43]
s.genes <- s.genes[!(s.genes %in% c("Mlf1ip", "Fam64a", "Hn1"))]
g2m.genes <- cc.genes[44:97]
g2m.genes <- g2m.genes[!(g2m.genes %in% c("Mlf1ip", "Fam64a", "Hn1"))]

object.list <- SplitObject(merge_obj, split.by = "experiment")
object.list <- lapply(X = object.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
})

merge_obj <- object.list[[1]]
for (i in 2:length(object.list)) {
    merge_obj <- merge(merge_obj, object.list[i])
}

# Split cells by experimental batch and cycling/non-cycling for integration

merge_obj@meta.data$integrate.ident <- merge_obj@meta.data$experiment
merge_obj@meta.data$integrate.ident[merge_obj@meta.data$Phase != "G1"] <- paste0(
    merge_obj@meta.data$integrate.ident[merge_obj@meta.data$Phase != "G1"], "_cycling")
table(merge_obj@meta.data$integrate.ident)
object.list <- SplitObject(merge_obj, split.by = "integrate.ident")
object.list <- lapply(X = object.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

# Select integration features and filter

features <- SelectIntegrationFeatures(object.list = object.list)
features <- features[!grepl("^Tr.v", features)]
features <- features[!grepl("^Ig.v", features)]
features <- features[!(features %in% cc.genes)]
features <- features[!grepl("^mt-", features)]

# Run PCA

object.list <- lapply(X = object.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find integration anchors using non-cycling cells as reference

anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features, reference = c(1,3)
                                  , reduction = "rpca", dims = 1:50)

 # Integrate data

merge_obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

# Find clusters and run UMAP

merge_obj.integrated <- merge_obj.integrated %>% 
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>% FindNeighbors(reduction = "pca", dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = 0.45, verbose = FALSE) %>% 
    RunUMAP(dims = 1:10, min.dist = 0.1, verbose = FALSE)

# Save R object

saveRDS(merge_obj, "data/seurat_obj.rds")