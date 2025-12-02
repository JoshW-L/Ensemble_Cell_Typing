print(getwd())
source("0_Functions.R")

# Load Seurat objects
seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/2.0 Method_1/2c Endo_Immune_Other_subclustered.Rds")

# ----------------------------  make new column to allow for FindMarkers ---------------------------------
# Function to add unified subcluster column
add_combined_subcluster <- function(seu) {
  subcol_names <- c("Immune_Clusters", "Endothelial_Clusters", "Other_Clusters")
  
  # Start with NA
  seu$subcluster_comb <- NA
  
  # Assign values for each cluster type
  if ("Immune_Clusters" %in% colnames(seu@meta.data)) {
    seu$subcluster_comb[!is.na(seu$Immune_Clusters)] <- paste0("Immune_", seu$Immune_Clusters[!is.na(seu$Immune_Clusters)])
  }
  if ("Endothelial_Clusters" %in% colnames(seu@meta.data)) {
    seu$subcluster_comb[!is.na(seu$Endothelial_Clusters)] <- paste0("Endothelial_", seu$Endothelial_Clusters[!is.na(seu$Endothelial_Clusters)])
  }
  if ("Other_Clusters" %in% colnames(seu@meta.data)) {
    seu$subcluster_comb[!is.na(seu$Other_Clusters)] <- paste0("Other_", seu$Other_Clusters[!is.na(seu$Other_Clusters)])
  }
  
  return(seu)
}

# Apply to all Seurat objects in the list
seu_split <- lapply(seu_split, add_combined_subcluster)

# ----------------------------  Check new column works ---------------------------------

table(seu_split$Imp2_P@meta.data$subcluster_comb)
table(seu_split$Imp2_P@meta.data$Other_Clusters)
table(seu_split$Imp2_P@meta.data$Immune_Clusters)
table(seu_split$Imp2_P@meta.data$Endothelial_Clusters)

# ----------------------------  Find Markers ---------------------------------

# Directory to save plots
output_dir <- "Files/Transfer Report Pipeline Mk3/2.0 Method_1/2d_Subcluster_Markers"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Run marker identification
for (obj_name in names(seu_split)) {
  cat("Processing:", obj_name, "\n")
  obj <- seu_split[[obj_name]]
  
  Idents(obj) <- obj$subcluster_comb
  
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

# ----------------------------  Reference Markers Nomura Filtered  ---------------------------------

# Load cell type marker data
ref <- read.csv("Files/Ref/Shoaib_ref.csv")

# ----------------------------  Fisher scoring  ---------------------------------

# Directories for output
out_fisher  <- "Files/Transfer Report Pipeline Mk3/2.0 Method_1/2d_Subcluster_Alignment_Fisher"
out_png     <- "Figures/Transfer Report Pipeline Mk3/2.0 Method_1/2d_Subcluster_Alignment"

dir.create(out_fisher, recursive = TRUE, showWarnings = FALSE)
dir.create(out_png, recursive = TRUE, showWarnings = FALSE)

# Helper function to compute Fisher enrichment 
fisher_score <- function(set1, set2, background) {
  a <- length(intersect(set1, set2)) # overlap 
  b <- length(set1) - a # in set1 only 
  c <- length(set2) - a # in set2 only 
  d <- length(background) - (a+b+c) # neither 
  mat <- matrix(c(a, b, c, d), nrow = 2) 
  fisher.test(mat, alternative = "greater")$p.value } 

# Map cluster prefixes to categories 
category_map <- list( 
  Other = c("Cancer", "Normal"), 
  Immune = c("Immune"), 
  Endothelial = c("Vasculature") )

# Loop through Seurat objects
for (obj_name in names(seu_split)) {
  cat("Processing:", obj_name, "\n")
  obj <- seu_split[[obj_name]]
  
  # Read its markers
  markers_file <- file.path("Files/Transfer Report Pipeline Mk3/2.0 Method_1/2d_Subcluster_Markers",
                            paste0(obj_name, "_cluster_markers.csv"))
  markers <- read.csv(markers_file)
  
  # Unique gene list per cluster
  cluster_genes <- markers %>%
    group_by(cluster) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  
  results_fisher  <- list()
  
  # Background = all genes in the Seurat object
  background <- rownames(obj)
  
  for (i in seq_len(nrow(cluster_genes))) {
    clust <- cluster_genes$cluster[i]
    gene_set <- cluster_genes$genes[[i]]
    
    prefix <- strsplit(clust, "_")[[1]][1]
    valid_cats <- category_map[[prefix]]
    ref_sub <- ref %>% filter(Category %in% valid_cats)
    
    for (ctype in unique(ref_sub$Broad_Cell_Type)) {
      ref_genes <- ref_sub %>% filter(Broad_Cell_Type == ctype) %>% pull(Gene) %>% unique()
      pval <- fisher_score(gene_set, ref_genes, background)
      
      results_fisher[[length(results_fisher)+1]] <- data.frame(
        Cluster = clust, Broad_Cell_Type = ctype, Pval = pval,
        NegLog10P = -log10(pval + 1e-300)
      )
    }
  }
  
  # Save CSVs
  df_fisher  <- bind_rows(results_fisher)
  
  write.csv(df_fisher,  file.path(out_fisher,  paste0(obj_name, "_fisher.csv")), row.names = FALSE)
  
  p_fish <- ggplot(df_fisher, aes(x = Cluster, y = Broad_Cell_Type, fill = NegLog10P)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Fisher")
  
  # Combine heatmaps side-by-side
  combined <- p_fish
  
  outfile <- file.path(out_png, paste0(obj_name, "_heatmaps.png"))
  ggsave(outfile, combined, width = 7, height = 6, dpi = 300)
  message("Saved combined heatmap: ", outfile)
}

# ----------------------------  Label Seurats based on best scores  ---------------------------------

# Directories
out_fisher  <- "Files/Transfer Report Pipeline Mk3/2.0 Method_1/2d_Subcluster_Alignment_Fisher"

# Function to read CSVs and extract best assignment per cluster
get_best_assignments <- function(seu_names, out_dir, score_type = c("Fisher")) {
  score_type <- match.arg(score_type)
  best_assignments <- list()
  
  for (obj_name in seu_names) {
    file_path <- file.path(out_dir, paste0(obj_name, "_fisher.csv"))    
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    df <- read.csv(file_path)
    
    if (score_type == "Fisher") {
      best <- df %>%
        dplyr::group_by(Cluster) %>%
        dplyr::slice_min(Pval, n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::select(Cluster, Broad_Cell_Type)
    }
    
    best_assignments[[obj_name]] <- best
  }
  
  return(best_assignments)
}

# Get Seurat object names
seu_names <- names(seu_split)

# Extract best assignments per object
best_fisher  <- get_best_assignments(seu_names, out_fisher, score_type = "Fisher")

# Function to label Seurat metadata safely
label_seurat_clusters <- function(seu_list, best_assignments, score_type = c("Fisher")) {
  score_type <- match.arg(score_type)
  
  for (obj_name in names(seu_list)) {
    obj <- seu_list[[obj_name]]
    assigned <- best_assignments[[obj_name]]
    
    if (is.null(assigned)) next
    
    cluster_to_cell <- setNames(as.character(assigned$Broad_Cell_Type), as.character(assigned$Cluster))
    
    # Assign new column in metadata using match() to avoid duplicate subscripts
    # Safe assignment
    col_name <- paste0(score_type, "_assignment")
    obj@meta.data[[col_name]] <- sapply(obj$subcluster_comb, function(cl) {
      if(cl %in% names(cluster_to_cell)) cluster_to_cell[cl] else NA
    })
    
    seu_list[[obj_name]] <- obj
  }
  
  return(seu_list)
}

# Apply labeling
seu_split <- label_seurat_clusters(seu_split, best_fisher, "Fisher")

# ----------------------------  UMAP of labelling  ---------------------------------

outdir <- "Figures/Transfer Report Pipeline Mk3/2.0 Method_1/2d_UMAP_cell_types"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Unified palettes 
# (1) Subcluster palette: show only 3 categories in the legend
subcluster_palette <- c(
  "Other" = "#1f78b4",       # blue
  "Immune" = "#31a354",      # green
  "Endothelial" = "#d7301f", # red
  "NA" = "grey80"
)

# (2) Full cell-type palette for Jaccard/Fisher (malignant+neural = purple→blue/green,
#     immune = green→yellow, endothelial/pericyte = red→orange). Also align "Other" to blue.
celltype_colors <- c(
  # Immune cells (green → yellow)
  "B cell" = "#2ca25f",
  "Plasma B" = "#a1d99b",
  "DC" = "#238b45",
  "Mast cell" = "#c7e9c0",
  "Monocyte" = "#41ab5d",
  "NK cell" = "#addd8e",
  "T cell" = "#006d2c",
  "Macrophage" = "#78c679",
  
  # Brain-resident / neural (blue → purple)
  "Microglia" = "#4eb3d3",
  "Astrocyte" = "#0571b0",
  "Oligodendrocyte" = "#74add1",
  "OPC normal" = "#9e9ac8",
  "Neuron" = "#3182bd",
  "Radial glial" = "#6baed6",
  
  # Malignant / neural stem-like (purple shades)
  "MES1" = "#6a51a3",
  "MES2" = "#807dba",
  "AC" = "#8da0cb",
  "OPC" = "#88419d",
  "NPC1" = "#225ea8",
  "NPC2" = "#1d91c0",
  
  # Endothelial / perivascular
  "Mural cell" = "#f16913",
  "Endothelial" = "#d7301f"
)

# ---- Map subcluster labels ("Other_*", "Immune_*", "Endothelial_*") to the 3 categories
map_subcluster_to_three <- function(x) {
  out <- ifelse(grepl("^Immune", x, ignore.case = TRUE), "Immune",
                ifelse(grepl("^Endothelial", x, ignore.case = TRUE), "Endothelial",
                       ifelse(grepl("^Other", x, ignore.case = TRUE), "Other", NA_character_)))
  out[is.na(out)] <- "Other"  # if anything unexpected appears, treat it as "Other" by default
  out
}

# ---- Safe DimPlot with named palettes (uses manual scales, not DimPlot's 'cols' argument)
safe_dimplot <- function(seu, group_col, title_txt, mode = c("subcluster", "celltype")) {
  mode <- match.arg(mode)
  
  if (!group_col %in% colnames(seu@meta.data)) {
    warning(sprintf("Column '%s' not found in meta.data for %s",
                    group_col, seu@project.name %||% "SeuratObj"))
    return(ggplot() + ggtitle(sprintf("%s (missing)", title_txt)) + theme_void())
  }
  
  vals <- seu[[group_col]][,1]
  vals[is.na(vals)] <- "NA"
  
  if (mode == "subcluster") {
    mapped <- map_subcluster_to_three(vals)
    present_levels <- unique(mapped)
    palette <- subcluster_palette[names(subcluster_palette) %in% present_levels]
  } else {
    mapped <- vals
    present_levels <- unique(mapped)
    # keep only colors for present levels; fallback to grey for anything unexpected
    palette <- celltype_colors[names(celltype_colors) %in% present_levels]
    missing <- setdiff(present_levels, names(palette))
    if (length(missing) > 0) {
      pal_extra <- rep("grey80", length(missing))
      names(pal_extra) <- missing
      palette <- c(palette, pal_extra)
    }
  }
  
  mapped[is.na(mapped)] <- "NA"
  seu[[paste0(group_col, "_mapped")]] <- factor(mapped, levels = present_levels)
  
  umap_name <- if ("umap" %in% Reductions(seu)) "umap" else if ("UMAP" %in% Reductions(seu)) "UMAP" else NULL
  if (is.null(umap_name)) {
    warning(sprintf("No UMAP reduction found for %s", seu@project.name %||% "SeuratObj"))
    return(ggplot() + ggtitle("UMAP missing") + theme_void())
  }
  
  p <- DimPlot(
    seu,
    reduction = umap_name,
    group.by = paste0(group_col, "_mapped"),
    label = TRUE,
    repel = TRUE
  ) + ggtitle(title_txt) +
    theme(plot.title = element_text(hjust = 0.5))
  
  p + scale_color_manual(values = palette, drop = TRUE) +
    scale_fill_manual(values = palette, drop = TRUE)
}

# Triptych creator
make_umap_triptych <- function(seu_list, outdir) {
  for (nm in names(seu_list)) {
    seu <- seu_list[[nm]]
    
    p1 <- safe_dimplot(seu, "subcluster_comb", "Subclusters", mode = "subcluster")
    p3 <- safe_dimplot(seu, "Fisher_assignment", "Fisher", mode = "celltype")
    
    combined <- (p1 | p3) + plot_layout(guides = "collect") &
      theme(legend.position = "right")
    
    outfile <- file.path(outdir, paste0(nm, "_UMAP_subcluster_Fisher.png"))
    ggsave(outfile, combined, width = 18, height = 6, dpi = 300)
    message("Saved: ", outfile)
  }
}

#  Run
make_umap_triptych(seu_split, outdir)


# ----------------------------  Save Seurats ---------------------------------
saveRDS(seu_split, "Seurat Objects/Transfer Report Pipeline Mk3/2.0 Method_1/2d Cell_labelled.Rds")

