print(getwd())
source("0_Functions.R")

# Load Seurat objects
seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/1.0 QC/sObj.Rds")

# ----------------------------  Normalise and Scale ---------------------------------

Normalise <- function(sObj, sample_name){
  sObj <- NormalizeData(sObj)
  sObj <- ScaleData(object = sObj, features = rownames(sObj))
  return(sObj)
}

# Apply the function to each Seurat object in seu_split
for (sample_name in names(seu_split)) {
  cat("Normalising:", sample_name, "\n")
  seu_split[[sample_name]] <- Normalise(seu_split[[sample_name]], sample_name)
}

# ----------------------------  Define marker genes by category ---------------------------------
marker_list <- list(
  Immune = c("PTPRC"),
  Endothelial = c("CDH5", "VWF")
)

marker_df <- do.call(rbind, lapply(names(marker_list), function(ct) {
  data.frame(Gene = marker_list[[ct]], CellType = ct, stringsAsFactors = FALSE)
}))

# ----------------------------  Nonmalignant marker expression ---------------------------------

plot_umap_and_barcharts_interleaved <- function(seu, seu_name, marker_df, save_dir) {
  genes_in_obj <- rownames(seu)
  match_genes <- marker_df[sapply(marker_df$Gene, function(g) {
    any(tolower(g) == tolower(genes_in_obj))
  }), ]
  
  if (nrow(match_genes) == 0) {
    message(paste("No marker genes found in", seu_name))
    return(NULL)
  }
  
  # Build interleaved plots (bar then UMAP for each marker)
  interleaved_plots <- lapply(1:nrow(match_genes), function(i) {
    gene <- match_genes$Gene[i]
    ct <- match_genes$CellType[i]
    matched_gene <- genes_in_obj[match(tolower(gene), tolower(genes_in_obj))]
    
    # Barchart
    expr <- FetchData(seu, vars = matched_gene, slot = "data")
    expr_df <- data.frame(Cell = rownames(expr), Expression = expr[,1])
    expr_df <- expr_df %>% arrange(Expression) %>% mutate(Cell = factor(Cell, levels = Cell))
    
    bar_plot <- ggplot(expr_df, aes(x = Cell, y = Expression)) +
      geom_bar(stat = "identity") +
      labs(title = paste0(gene, " (", ct, ") - Barchart"),
           x = "Cells (ordered by expression)", y = "Expression") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    # UMAP
    umap_plot <- FeaturePlot(seu, features = matched_gene, reduction = "umap") +
      ggtitle(paste0(gene, " (", ct, ") - UMAP"))
    
    list(bar_plot, umap_plot)
  }) |> unlist(recursive = FALSE)  # Flatten list
  
  # Arrange 4 across
  final_plot <- wrap_plots(interleaved_plots, ncol = 4)
  
  ggsave(
    file.path(save_dir, paste0(seu_name, "_Interleaved_Barchart_UMAP.png")),
    plot = final_plot,
    width = 16,
    height = 4 * ceiling(length(interleaved_plots) / 4),
    dpi = 300
  )
}

# Output directory
combo_dir <- "Figures/Transfer Report Pipeline Mk3/5.0 Method_4/5a Immune_Endo_Labelling/Barcharts_UMAPs"
dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)

# Run for all samples
for (sample_name in names(seu_split)) {
  cat("Plotting interleaved UMAPs + Barcharts for:", sample_name, "\n")
  plot_umap_and_barcharts_interleaved(seu_split[[sample_name]], sample_name, marker_df, combo_dir)
}

# ----------------------------  Label cells ---------------------------------

expr_threshold <- 0.1  # "expressed" threshold
coexp_min <- list(
  Immune = 1,
  Endothelial = 1
)

classify_cells_immune_endo <- function(seu_obj) {
  
  expr_mat <- GetAssayData(seu_obj, slot = "data")
  cell_names <- colnames(expr_mat)
  
  # Count markers above threshold for each category
  type_counts <- sapply(names(marker_list), function(type) {
    markers <- marker_list[[type]]
    markers <- markers[markers %in% rownames(expr_mat)]
    apply(expr_mat, 2, function(x) sum(x[markers] > expr_threshold, na.rm = TRUE))
  })
  
  if (is.vector(type_counts)) {
    type_counts <- matrix(type_counts, ncol = 1, dimnames = list(cell_names, names(marker_list)))
  }
  
  # Logical matches
  matches <- sapply(names(marker_list), function(type) {
    type_counts[, type] >= coexp_min[[type]]
  })
  
  # Priority: Immune > Endothelial > Other
  assignment <- ifelse(matches[, "Immune"], "Immune",
                       ifelse(matches[, "Endothelial"], "Endothelial", "Other"))
  
  seu_obj$Immune_Endo_Status <- assignment
  
  return(seu_obj)
}

# Apply to all Seurat objects
for (sample_name in names(seu_split)) {
  cat("Classifying Immune/Endo for:", sample_name, "\n")
  seu_split[[sample_name]] <- classify_cells_immune_endo(seu_split[[sample_name]])
}

# ----------------------------  Labelled UMAPs ---------------------------------


# Create save directory
umap_dir <- "Figures/Transfer Report Pipeline Mk3/5.0 Method_4/5a Immune_Endo_Labelling/UMAPs"
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)

# Loop over each Seurat object and save UMAP plot
for (sample_name in names(seu_split)) {
  
  seu_obj <- seu_split[[sample_name]]
  
  # Ensure UMAP exists
  if (!"umap" %in% names(seu_obj@reductions)) {
    seu_obj <- RunUMAP(seu_obj, dims = 1:30)  # Adjust dims as needed
  }
  
  # Make UMAP plot
  p <- DimPlot(
    seu_obj,
    reduction = "umap",
    group.by = "Immune_Endo_Status",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste0("UMAP - ", sample_name, " (immune - endothelial)")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  # Save as PNG
  ggsave(
    filename = file.path(umap_dir, paste0(sample_name, "_UMAP_immune_endothelial.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}
# ----------------------------  Save Seurats ---------------------------------
dir.create("Seurat Objects/Transfer Report Pipeline Mk3/5.0 Method_4", 
           showWarnings = FALSE, 
           recursive = TRUE)

saveRDS(seu_split, "Seurat Objects/Transfer Report Pipeline Mk3/5.0 Method_4/5a Immune_Endo_cells_labelled.Rds")
