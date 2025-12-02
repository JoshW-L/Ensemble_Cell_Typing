print(getwd())
source("0_Functions.R")

# Load Seurat objects
seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/3.0 Method_2/3a Immune_Endo_cells_labelled.Rds")

# ----------------------------  Subcluster Endo, Immune and Other separately ---------------------------------

# Directory for UMAPs
out_dir <- "Figures/Transfer Report Pipeline Mk3/3.0 Method_2/3b Immune_Endo_Others_subclustering/UMAPs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

process_category <- function(seu, category, label_col, sample_name) {
  n_cells <- sum(seu[[label_col]] == category)
  message("Processing: ", sample_name, " - ", category, " (", n_cells, " cells)")
  
  # ⬅️ Skip if no cells
  if (n_cells == 0) {
    message("⚠️ Skipping ", category, " for ", sample_name, " (0 cells)")
    return(NULL)
  }
  
  clust_col_name <- paste0(category, "_Clusters")
  
  # Subset cells
  seu_sub <- subset(seu, subset = !!as.name(label_col) == category)
  
  # Default: all cluster = "1"
  clust_col <- rep(NA, ncol(seu))
  clust_col[colnames(seu) %in% colnames(seu_sub)] <- "1"
  seu[[clust_col_name]] <- clust_col
  
  # Adjust clustering params
  if (n_cells < 30) {
    dims_to_use <- 1:2; npcs <- 2; res <- 0.1
  } else if (n_cells < 50) {
    dims_to_use <- 1:2; npcs <- 2; res <- 1
  } else if (n_cells < 60) {
    dims_to_use <- 1:2; npcs <- 2; res <- 2
  } else if (n_cells < 100) {
    dims_to_use <- 1:10; npcs <- 10; res <- 2
  } else if (n_cells < 150) {
    dims_to_use <- 1:20; npcs <- 20; res <- 2
  } else if (n_cells < 500) {
    dims_to_use <- 1:30; npcs <- 30; res <- 2
  } else {
    dims_to_use <- 1:50; npcs <- 50; res <- 2
  }
  
  # Try clustering pipeline
  res <- tryCatch({
    seu_sub <- NormalizeData(seu_sub, verbose = FALSE)
    seu_sub <- FindVariableFeatures(seu_sub, verbose = FALSE)
    seu_sub <- ScaleData(seu_sub, verbose = FALSE)
    seu_sub <- RunPCA(seu_sub, npcs = npcs, verbose = FALSE)
    seu_sub <- FindNeighbors(seu_sub, dims = dims_to_use, verbose = FALSE)
    seu_sub <- FindClusters(seu_sub, resolution = res, algorithm = 4, verbose = FALSE)
    
    # Handle n.neighbors properly
    n_neighbors <- min(30, n_cells - 1)
    seu_sub <- RunUMAP(seu_sub, dims = dims_to_use, n.neighbors = n_neighbors, verbose = FALSE)
    
    # Store clusters
    clust_col[colnames(seu) %in% colnames(seu_sub)] <- as.character(Idents(seu_sub))
    seu[[clust_col_name]] <- clust_col
    
    # Make plot
    p <- DimPlot(seu_sub, label = TRUE) + 
      ggtitle(paste0(category, " (", n_cells, " cells)"))
    list(seu = seu, plot = p)
  }, error = function(e) {
    warning("⚠️ Skipping PCA/UMAP/clustering for ", category, " in ", sample_name, 
            ": ", conditionMessage(e))
    # Fallback: fake UMAP with all cluster 1
    seu_sub$DummyCluster <- "1"
    p <- tryCatch(
      DimPlot(seu_sub, group.by = "DummyCluster") + 
        ggtitle(paste0(category, " (", n_cells, " cells; fallback=1)")),
      error = function(e2) ggplot() + ggtitle(paste(category, "(fallback=1)"))
    )
    list(seu = seu, plot = p)
  })
  
  return(res)
}

# Loop through Seurat list
label_col <- "Immune_Endo_Status"

for (sample_name in names(seu_split)) {
  seu <- seu_split[[sample_name]]
  plots <- list()
  
  for (cat in c("Endothelial", "Immune", "Other")) {
    res <- process_category(seu, cat, label_col, sample_name)
    if (!is.null(res)) {
      seu <- res$seu
      plots[[cat]] <- res$plot
    }
  }
  
  # Combine only available plots
  if (length(plots) > 0) {
    combined_plot <- wrap_plots(plots, ncol = length(plots))
    ggsave(
      filename = file.path(out_dir, paste0(sample_name, "_UMAPs.png")),
      plot = combined_plot,
      width = 7 * length(plots), height = 5, dpi = 300
    )
  } else {
    message("⚠️ No categories available for ", sample_name, " — skipping plot")
  }
  
  seu_split[[sample_name]] <- seu
}

# ----------------------------  Save Seurats ---------------------------------
saveRDS(seu_split, "Seurat Objects/Transfer Report Pipeline Mk3/3.0 Method_2/3b Endo_Immune_Other_subclustered.Rds")
