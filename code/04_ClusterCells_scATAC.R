# Load packages

library(ArchR)
library(Seurat)
library(Cairo)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1)
addArchRThreads(threads = 20) 
addArchRGenome("mm10")

#-----------------
# Cluster cells
#-----------------

# Load project

proj <- loadArchRProject("ArchR_output")

# Filter doublets

proj <- filterDoublets(ArchRProj = proj)

# Add iterative LSI

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)

# Add clusters

proj <- addClusters(input = proj, reducedDims = "IterativeLSI", force = TRUE, dimsToUse = 1:10, resolution = 0.4)

# Add UMAP

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE, dimsToUse = 1:10, minDist = 0.1)

# Add gene integration matrix

seRNA <- readRDS("data/seurat_obj.rds")
seRNA[["ident"]] <- Idents(object = seRNA)
# Only include relevant samples/phenotypes
seRNA <- subset(seRNA, cells = Cells(seRNA)[seRNA$lcmv %in% c('Chronic_D8', 'Chronic_D21', 'Acute_D8') & seRNA$tissue == "Spleen"])
seRNA <- subset(seRNA, cells = Cells(seRNA)[!(seRNA$ident %in% c('Memory'))])
seRNA_exp <- SummarizedExperiment(seRNA@assays$RNA@counts, colData = seRNA@meta.data)
proj <- addGeneIntegrationMatrix(force = TRUE,
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA_exp,
    addToArrow = TRUE,
    groupRNA = "ident",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

# Add impute weights

proj <- addImputeWeights(proj)

# Call peaks and add peak matrix

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters", force = TRUE)
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)
proj <- addPeakMatrix(proj)

# Add TF deviation matrix

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

# Add co-accessibility and peak to gene links
proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)
proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)

# Save ArchR object

proj <- saveArchRProject(ArchRProj = proj)