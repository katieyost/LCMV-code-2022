# Load packages

library(Seurat)
library(SingleR)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 2400 * 1024^2)
library(DoubletFinder)
options(future.rng.onMisuse="ignore")

#-----------------
# Functions
#-----------------

add_clonotype <- function(seurat_object, clonotypes_path, contigs_path, lib_type
) {
    contigs <- read.csv(contigs_path)
    contigs <- contigs[!duplicated(contigs$barcode), ]

    contigs <- contigs[,c("barcode", "raw_clonotype_id")]
    names(contigs)[names(contigs) == "raw_clonotype_id"] <- "clonotype_id"

    clonotypes <- read.csv(clonotypes_path)

    contigs <- merge(contigs, clonotypes[, c("clonotype_id", "cdr3s_aa", "cdr3s_nt")], by = "clonotype_id")

    rownames(contigs) <- contigs$barcode
    contigs <- contigs[,c(1,3,4)]
    colnames(contigs) <- paste0(lib_type, "_", colnames(contigs))

    seurat_object <- AddMetaData(object=seurat_object, metadata=contigs)
    return(seurat_object)
}

#-----------------
# Preprocess from cellranger outputs
#-----------------

# Read in gene expression matrices and create list of Seurat objects

objects <- c()
dir <- "data/"
for (folder in list.files(dir)) {
    if (grepl("_gex$", folder)) {
        data <- Read10X(data.dir = paste0(dir, folder, "/outs/filtered_feature_bc_matrix/"))
        if (typeof(data) == 'list') {
            obj <- CreateSeuratObject(counts = data$`Gene Expression`, project = folder, min.cells = 0, min.features = 0)
        } else {
            obj <- CreateSeuratObject(counts = data, project = folder, min.cells = 0, min.features = 0)
        }
        
        objects <- c(obj, objects)
        names(objects)[1] <- folder
    }
}

# Read in TCR data and add to corresponding Seurat object as metadata

for (folder in list.files(dir)) {
    if (grepl("*_tcr*", folder)) {
        name <- folder
        
        if (gsub("tcr", "gex", name) %in% names(objects) & 
            file.exists(paste0(dir, folder,"/outs/clonotypes.csv")) &
           file.exists(paste0(dir, folder,"/outs/filtered_contig_annotations.csv"))){
            obj <- objects[[gsub("tcr", "gex", name)]]
            obj <- add_clonotype(obj,paste0(dir, folder,"/outs/clonotypes.csv"),
                             paste0(dir, folder,"/outs/filtered_contig_annotations.csv"),
                             "tcr")
            objects[[gsub("tcr", "gex", name)]] <- obj
            }
    }
}

# Predict doublets for each object using DoubletFinder

for (i in 1:length(objects)) {
    message(i)
    objects[[i]] <- NormalizeData(objects[[i]])
    objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", nfeatures = 2000)
    objects[[i]] <- ScaleData(objects[[i]], verbose = FALSE)
    objects[[i]] <- RunPCA(objects[[i]], features = VariableFeatures(object = objects[[i]]), verbose = FALSE)
    objects[[i]] <- FindNeighbors(objects[[i]] , dims = 1:10)
    objects[[i]] <- FindNeighbors(objects[[i]], dims = 1:10)
    objects[[i]] <- FindClusters(objects[[i]])
    sweep.res.list <- paramSweep_v3(objects[[i]], PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    homotypic.prop <- modelHomotypic(Idents(objects[[i]])) 
    # Use predicted doublet rate of 0.8% doublets per 1,000 cells
    nExp_poi <- round((0.008/1000)*length(Cells(objects[[i]]))*length(Cells(objects[[i]])))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    objects[[i]] <- doubletFinder_v3(objects[[i]], PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi
                                     , reuse.pANN = FALSE, sct = FALSE)
    objects[[i]] <- doubletFinder_v3(objects[[i]], PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj
                                     , reuse.pANN = colnames(objects[[i]]@meta.data)[
                                         grepl("^pANN_", colnames(objects[[i]]@meta.data))], sct = FALSE)
}

for (i in 1:length(objects)) {
    colnames(objects[[i]]@meta.data)[length(colnames(objects[[i]]@meta.data)) - 2] <- "pANN"
    colnames(objects[[i]]@meta.data)[length(colnames(objects[[i]]@meta.data)) - 1] <- "DF.classifications_1"
    colnames(objects[[i]]@meta.data)[length(colnames(objects[[i]]@meta.data))] <- "DF.classifications_2"
}

# Merge objects and make unique cell names

merge_obj <- objects[[1]]
merge_obj <- RenameCells(merge_obj, add.cell.id = names(objects)[1])
for (i in 2:length(objects)) {
    merge_obj <- merge(merge_obj, objects[i], add.cell.ids = c("", names(objects)[i]))
}
merge_obj <- RenameCells(merge_obj, new.names = gsub("^_*", "", Cells(merge_obj)))

# Add percent mitochondrial reads to metadata

merge_obj[["percent.mt"]] <- PercentageFeatureSet(merge_obj, pattern = "^mt-")

# Add doublet classification metadata

merge_obj@meta.data$DF.classifications <- "Singlet"
merge_obj@meta.data$DF.classifications[merge_obj@meta.data$DF.classifications_1 == "Doublet" |
                                      merge_obj@meta.data$DF.classifications_2 == "Doublet"
                                      ] <- "LowConfidenceDoublet"
merge_obj@meta.data$DF.classifications[merge_obj@meta.data$DF.classifications_1 == "Doublet" &
                                      merge_obj@meta.data$DF.classifications_2 == "Doublet"
                                      ] <- "HighConfidenceDoublet"
table(merge_obj@meta.data$DF.classifications)

# Add experimental batch metadata

merge_obj@meta.data$experiment <- "exp1"
merge_obj@meta.data$experiment[merge_obj@meta.data$orig.ident %in% 
                               c("Spleen_Chronic_D21_PD1Pos"
                                 , "Spleen_Chronic_D21_CX3CR1Pos"
                                 , "Spleen_Chronic_D21_CX3CR1Neg_Slamf6Pos"
                                 , "Spleen_Chronic_D21_CX3CR1Neg"
                                 , "Spleen_Chronic_D21_TetPos_rep1"
                                 , "Spleen_Chronic_D21_TetNeg_rep1")] <- "exp2"

# Add cell type predictions with SingleR

ref <- MouseRNAseqData()
pred <- SingleR(test = merge_obj@assays$RNA@data, ref = ref, labels = ref$label.main)
pred.labels <- data.frame(pred[,"labels", drop = FALSE])
colnames(pred.labels) <- "SingleR_labels"
merge_obj <- AddMetaData(merge_obj, pred.labels)

# Save R object

saveRDS(merge_obj, "data/seurat_obj.rds")