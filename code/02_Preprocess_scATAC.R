# Load packages

library(ArchR)
library(Seurat)
library(Cairo)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1)
addArchRThreads(threads = 20) 
addArchRGenome("mm10")

#-----------------
# Preprocess from cellranger outputs (fragment files)
#-----------------

# Input fragment files
inputFiles <- list.files("data/", pattern = "*.tsv.gz$", full.names = TRUE)
names(inputFiles) <- gsub("_fragments", "", gsub(".tsv.gz", "", basename(inputFiles)))

setwd(dir = "ArchR_output")

# Create arrow files

ArrowFiles <- createArrowFiles(force = TRUE,
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Add doublet scores

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

# Create ArchR project

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR_output",
  copyArrows = TRUE 
)

# Save ArchR object

proj <- saveArchRProject(ArchRProj = proj)