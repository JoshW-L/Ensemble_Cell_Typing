print(getwd())
source("0_Functions.R")

# Load Seurat objects
seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/1.0 QC/sObj.Rds")

# -----------------  Filter Seurats for genes present in Nomura list -----------------------
# Load cell type marker data
ref <- read.csv("Files/Ref/Shoaib_ref.csv")

# Get the list of reference genes
ref_genes <- unique(ref$Gene)

# Function to filter one Seurat object
filter_seurat_by_ref <- function(seu, ref_genes) {
  common_genes <- intersect(rownames(seu), ref_genes)
  seu[common_genes, ]
}

seu_split <- lapply(seu_split, filter_seurat_by_ref, ref_genes = ref_genes)

# ----------------------------  Normalize, Scale and Cluster ---------------------------------

## pre-process data function
clusterCells <- function(sObj, sample_name){
  sObj <- NormalizeData(sObj)
  sObj <- ScaleData(object = sObj, features = rownames(sObj))
  sObj <- FindVariableFeatures(sObj)
  sObj <- RunPCA(object = sObj, npcs = 40)
  
  sObj <- FindNeighbors(object = sObj, dims = 1:40,  reduction = "pca")
  resolutions <- 2.0
  for (res in resolutions) {
    sObj <- FindClusters(object = sObj, resolution = res, algorithm = 4)
    sObj <- RunUMAP(object = sObj, dims = 1:40)
    
    # UMAP Plot for each resolution
    umap_plot <- DimPlot(sObj, reduction = "umap") + plot_annotation(title = paste("UMAP Plot -", sample_name))
    ggsave(filename = paste0("Figures/Transfer Report Pipeline Mk3/6.0 Cluster Methods/6b Initial_Clustering_Filtered/UMAP_Plot_", sample_name,".png"), plot = umap_plot)
  }
  
  return(sObj)
}

## Apply the function to each Seurat object in seu_split
for (sample_name in names(seu_split)) {
  seu_split[[sample_name]] <- clusterCells(seu_split[[sample_name]], sample_name)
}

## save data
dir.create("Seurat Objects/Transfer Report Pipeline Mk3/6.0 Cluster Methods", 
           showWarnings = FALSE, 
           recursive = TRUE)

saveRDS(seu_split, "Seurat Objects/Transfer Report Pipeline Mk3/6.0 Cluster Methods/6b Initial_Clustering_Filter.Rds")



