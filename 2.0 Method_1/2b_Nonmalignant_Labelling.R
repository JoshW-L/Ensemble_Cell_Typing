print(getwd())
source("0_Functions.R")

# Load List of Seurat objects that have been clustered
seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/2.0 Method_1/2a Initial_Clustering.Rds")

# ----------------------------  Define marker genes by category ---------------------------------
# List of Markers for Immune and Endo
marker_list <- list(
  Immune = c("PTPRC"),
  Endothelial = c("CDH5", "VWF")
)
# ----------------------------  Find Markers ---------------------------------

# Directory to save plots
output_dir <- "Files/Transfer Report Pipeline Mk3/2.0 Method_1/2b Nonmalignant_clusters_labelled/Nonmalignant_markers"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Run marker identification
for (obj_name in names(seu_split)) {
  cat("Processing:", obj_name, "\n")
  obj <- seu_split[[obj_name]]
  
  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    logfc.threshold = 0.1,
    min.pct = 0.1,
    test.use = "wilcox",
    assay = DefaultAssay(obj),
    p.adjust.method = "bonferroni"
  ) %>%
    filter(p_val_adj < 0.1)
  
  write.csv(markers, file.path(output_dir, paste0(obj_name, "_cluster_markers.csv")), row.names = FALSE)
}

# ----------------------------  Label each cluster by Nonmalignant marker ---------------------------------

# Reverse map marker genes to cell types
gene_to_type <- setNames(rep(names(marker_list), times = sapply(marker_list, length)),
                         unlist(marker_list))

# Data frame to collect multi-marker conflict info
conflict_log <- data.frame(
  sample = character(),
  cluster = character(),
  gene = character(),
  cell_type = character(),
  avg_log2FC = numeric(),
  p_val_adj = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(seu_split)) {
  cat("Labelling:", name, "\n")
  obj <- seu_split[[name]]
  marker_path <- file.path("Files/Transfer Report Pipeline Mk3/2.0 Method_1/2b Nonmalignant_clusters_labelled//Nonmalignant_markers", paste0(name, "_cluster_markers.csv"))
  markers <- read.csv(marker_path)
  
  pos_markers <- markers %>%
    filter(avg_log2FC > 0, p_val_adj < 0.1)
  
  cluster_labels <- rep("Other", length(unique(pos_markers$cluster)))
  names(cluster_labels) <- unique(pos_markers$cluster)
  
  for (clust in unique(pos_markers$cluster)) {
    clust_genes <- pos_markers %>% filter(cluster == clust)
    matched_genes <- clust_genes$gene[clust_genes$gene %in% names(gene_to_type)]
    
    if (length(matched_genes) > 0) {
      matched_types <- gene_to_type[matched_genes]
      unique_types <- unique(matched_types)
      
      if (length(unique_types) == 1) {
        cluster_labels[as.character(clust)] <- unique_types
      } else {
        top_gene <- clust_genes %>%
          filter(gene %in% names(gene_to_type)) %>%
          arrange(desc(avg_log2FC)) %>%
          head(1)
        
        selected_type <- gene_to_type[top_gene$gene]
        cluster_labels[as.character(clust)] <- selected_type
        
        message(sprintf("⚠️ Cluster %s in %s expresses multiple markers (%s), assigned to '%s' based on highest avg_log2FC.",
                        clust, name, paste(unique_types, collapse = ", "), selected_type))
        
        # Log conflicts for export
        clust_conflict_genes <- clust_genes %>%
          filter(gene %in% matched_genes) %>%
          mutate(
            sample = name,
            cluster = clust,
            cell_type = gene_to_type[gene]
          ) %>%
          dplyr::select(sample, cluster, gene, cell_type, avg_log2FC, p_val_adj)
        
        conflict_log <- rbind(conflict_log, clust_conflict_genes)
      }
    }
  }
  
  cluster_col <- "Nanostring_snn_res.2"
  cluster_ids <- as.character(obj[[cluster_col]][, 1])
  obj$Immune_Endo_Status <- ifelse(cluster_ids %in% names(cluster_labels),
                                cluster_labels[cluster_ids],
                                "Other")
  
  seu_split[[name]] <- obj
}

# Write conflict table to disk
conflict_path <- "Files/Transfer Report Pipeline Mk3/2.0 Method_1/2b Nonmalignant_clusters_labelled/multi_marker_cluster_conflicts.csv"
write.csv(conflict_log, conflict_path, row.names = FALSE)

# ----------------------------  UMAPs ---------------------------------

output_dir <- "Figures/Transfer Report Pipeline Mk3/2.0 Method_1/2b Nonmalignant_clusters_labelled/UMAP_Nonmalignant_markers"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (name in names(seu_split)) {
  cat("Creating UMAPs for:", name, "\n")
  obj <- seu_split[[name]]
  
  # Cluster UMAP
  p1 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Original Clusters")
  
  # FeaturePlots for each marker gene
  markers_to_plot <- unique(unlist(marker_list))
  available_genes <- markers_to_plot[markers_to_plot %in% rownames(obj)]
  
  fp_list <- lapply(available_genes, function(gene) {
    FeaturePlot(obj, features = gene, reduction = "umap") + ggtitle(gene)
  })
  
  # Final labeled UMAP
  Idents(obj) <- obj$Immune_Endo_Status
  p2 <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Cluster Labeling")
  
  # Arrange 3 columns: Clusters | Features | Final Labels
  all_plots <- wrap_plots(c(list(p1), fp_list, list(p2)), ncol = 3)
  
  ggsave(file.path(output_dir, paste0(name, "_nonmalignant_umap_grid.png")),
         plot = all_plots, width = 15, height = 4 * ceiling((length(fp_list)+2)/3))
}



# ----------------------------  Save Seurats ---------------------------------
saveRDS(seu_split, "Seurat Objects/Transfer Report Pipeline Mk3/2.0 Method_1/2b Nonmalignant_labelled.Rds")


# ----------------- Look at cell counts for Immune_Endo_Status ---------------------

# Create a table of counts per Immune_Endo_Status per Seurat object
malignancy_counts <- lapply(names(seu_split), function(sample_name) {
  obj <- seu_split[[sample_name]]
  counts <- table(obj$Immune_Endo_Status)
  df <- as.data.frame(counts)
  colnames(df) <- c("Immune_Endo_Status", "cell_count")
  df$sample <- sample_name
  return(df)
}) %>%
  bind_rows()

# Pivot wider so each Immune_Endo_Status is a column
malignancy_counts_wide <- malignancy_counts %>%
  pivot_wider(
    names_from = Immune_Endo_Status,
    values_from = cell_count,
    values_fill = 0
  )

# View results
malignancy_counts
malignancy_counts_wide
