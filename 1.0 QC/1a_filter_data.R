print(getwd())
source("0_Functions.R")

seu <- readRDS("/Users/bs23jmw/Library/CloudStorage/OneDrive-UniversityofLeeds/PhD/Data Analysis/Chapter 1/Seurat Objects/Same_Samples/2.0_Meta_Update/seurat_object_6k.Rds")

# First need to identify mitochondrial genes by alternative gene names 
# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Get mitochondrial genes along with synonyms
mito_genes <- getBM(
  attributes = c("external_gene_name", "hgnc_symbol", "external_synonym"),
  filters = "chromosome_name",
  values = "MT",
  mart = mart
)

# Match mitochondrial genes using exteNanostringl_synonym
mito_synonyms <- unique(mito_genes$external_synonym)

# Identify mitochondrial genes in Seurat objects
mito_genes <- rownames(seu)[rownames(seu) %in% mito_synonyms]

# Compute mitochondrial gene percentage
seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genes)

seu$Patient_P_R <- paste(seu$Patient_ID, seu$P_R, sep = "_")

sObj <- SplitObject(seu, split.by = "Patient_P_R")

# Create a data frame to store sample information
sampleInfo <- data.frame(
  fileName = names(sObj),  # Names of the split objects (patients)
  beforeFiltering  = sapply(sObj, function(obj) ncol(obj))  # Cell count for each patient
)

sObj <- lapply(sObj, function(obj) {
  if (!is.null(obj@images)) {  # Check if @images exists and is not NULL
    obj@images <- list()       # Clear the @images slot by setting it to an empty list
  }
  return(obj)
})

dir.create("Seurat Objects/Transfer Report Pipeline Mk3/1.0 QC", 
           showWarnings = FALSE, 
           recursive = TRUE)
saveRDS(sObj, "Seurat Objects/Transfer Report Pipeline Mk3/1.0 QC/sObj_pre_filter.rds")

# Filter metrics 
minCounts <- 100
minFeatures <- 75
maxMit <- 20


for (i in 1:nrow(sampleInfo)){
  filter <- sObj[[i]]$nCount_Nanostring > minCounts & 
    sObj[[i]]$nFeature_Nanostring > minFeatures & 
    sObj[[i]]$percent.mt < maxMit 

  sObj[[i]] <- sObj[[i]][ , filter]
  
}

sObj <- lapply(sObj, function(obj) {
  if (!is.null(obj@images)) {  # Check if @images exists and is not NULL
    obj@images <- list()       # Clear the @images slot by setting it to an empty list
  }
  return(obj)
})

# Remove FOVs with fewer than 50 cells AFTER filtering
cell_counts_after <- sapply(sObj, ncol)
sObj <- sObj[cell_counts_after >= 100]

# Update sampleInfo to reflect filtering and FOV exclusion
sampleInfo$afterFiltering <- cell_counts_after
sampleInfo <- sampleInfo[sampleInfo$afterFiltering >= 100, ]

# Save filtered objects
saveRDS(sObj, "Seurat Objects/Transfer Report Pipeline Mk3/1.0 QC/sObj.rds")

dir.create("Files/Transfer Report Pipeline Mk3/1.0 SampleInfo", 
           showWarnings = FALSE, 
           recursive = TRUE)

saveRDS(sampleInfo, "Files/Transfer Report Pipeline Mk3/1.0 SampleInfo/sampleInfo.rds")

