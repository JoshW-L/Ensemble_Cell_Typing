print(getwd())
source("0_Functions.R")

# ---------------------------- Helpers ---------------------------------
get_celltype_order <- function(cell_type_groups) {
  unlist(cell_type_groups)
}

force_matrix_order <- function(mat, ordered_types) {
  new_rows <- intersect(ordered_types, rownames(mat))
  new_cols <- intersect(ordered_types, colnames(mat))
  mat[new_rows, new_cols, drop = FALSE]
}

is_plottable_matrix <- function(mat) {
  !is.null(mat) && is.matrix(mat) && nrow(mat) > 0 && ncol(mat) > 0 && any(is.finite(mat))
}

fisher_score <- function(set1, set2, background) {
  a <- length(intersect(set1, set2))
  b <- length(set1) - a
  c <- length(set2) - a
  d <- max(length(background) - (a + b + c), 0)
  
  mat <- matrix(c(a, b, c, d), nrow = 2)
  if (any(mat < 0) || all(mat == 0)) return(NA_real_)
  fisher.test(mat, alternative = "greater")$p.value
}

# ---------------------------- Compute unified ordering ----------------------------
get_ordered_types <- function(ref, seu, cell_type_groups) {
  base_order <- get_celltype_order(cell_type_groups)
  all_types <- union(unique(ref$Broad_Cell_Type), unique(seu@meta.data$Resolved_Label))
  extra_types <- setdiff(all_types, base_order)
  c(base_order, sort(extra_types))
}

# ---------------------------- Load data ----------------------------
seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7b Clash_cell_labelled.Rds")
ref <- read.csv("Files/Ref Comparison/New Ref/Nomura Ref for 6k.csv")

marker_dirs <- list(
  Fisher_Method_1 = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Method_1",
  Fisher_Method_2 = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Method_2",
  Fisher_Method_3 = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Method_3",
  Fisher_Method_4 = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Method_4",
  Fisher_Draw_Breaker_1 = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Draw_Breaker_1",
  Fisher_Draw_Breaker_2 = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Draw_Breaker_1",
  ensemble_label_w_o_clash = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Ensemble_Vote",
  Resolved_Label  = "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients/Resolved"
)

outdir_base <- "Figures/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8e_Markers_By_Patients/Heatmaps"
dir.create(outdir_base, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- Cell type ordering ----------------------------
cell_type_groups <- list(
  Malignant   = c("AC", "MES", "LQ", "Hypoxia", "RP"),
  Neural      = c("ExN", "InN", "IPC", "CC", "Cilia", "Astrocyte"),
  Immune      = c("B cell", "Macrophage", "MIC", "T cell", "NK cell"),
  Vasculature = c("Endothelial", "Pericyte")
)

# ---------------------------- Compute similarity ----------------------------
compute_similarity_matrix <- function(marker_file, ref, ordered_types, top_n = 20) {
  if (!file.exists(marker_file)) return(NULL)
  
  markers <- read.csv(marker_file) %>% filter(cluster != "Unknown")
  if (nrow(markers) == 0) return(NULL)
  
  clusters <- unique(markers$cluster)
  sim_list <- list()
  
  for (cl in clusters) {
    tbl <- markers %>% filter(cluster == cl, gene %in% ref$Gene) %>% arrange(p_val_adj)
    if (!identical(top_n, "all")) tbl <- slice_head(tbl, n = as.integer(top_n))
    genes <- tbl$gene
    if (length(genes) == 0) next
    
    sim_row <- sapply(unique(ref$Broad_Cell_Type), function(ct) {
      rg <- ref %>% filter(Broad_Cell_Type == ct) %>% pull(Gene)
      if (length(rg) == 0) return(NA_real_)
      length(intersect(genes, rg)) / length(genes) * 100
    }, USE.NAMES = TRUE)
    
    sim_list[[cl]] <- sim_row
  }
  
  if (length(sim_list) == 0) return(NULL)
  
  sim_mat <- do.call(rbind, sim_list)
  rownames(sim_mat) <- names(sim_list)
  sim_mat <- force_matrix_order(sim_mat, ordered_types)
  if (!is_plottable_matrix(sim_mat)) return(NULL)
  
  sim_mat
}

# ---------------------------- Compute Fisher ----------------------------
compute_fisher_matrix <- function(marker_file, ref, ordered_types, top_n = 20) {
  if (!file.exists(marker_file)) return(NULL)
  
  markers <- read.csv(marker_file) %>% filter(cluster != "Unknown")
  if (nrow(markers) == 0) return(NULL)
  
  background <- union(unique(ref$Gene), unique(markers$gene))
  clusters <- unique(markers$cluster)
  res_list <- list()
  
  for (cl in clusters) {
    tbl <- markers %>% filter(cluster == cl, gene %in% ref$Gene) %>% arrange(p_val_adj)
    if (!identical(top_n, "all")) tbl <- slice_head(tbl, n = as.integer(top_n))
    genes <- tbl$gene
    if (length(genes) == 0) next
    
    scores <- sapply(unique(ref$Broad_Cell_Type), function(ct) {
      rg <- ref %>% filter(Broad_Cell_Type == ct) %>% pull(Gene)
      if (length(rg) == 0) return(NA_real_)
      p <- fisher_score(genes, rg, background)
      -log10(p + 1e-300)
    }, USE.NAMES = TRUE)
    
    res_list[[cl]] <- scores
  }
  
  if (length(res_list) == 0) return(NULL)
  
  mat <- do.call(rbind, res_list)
  rownames(mat) <- names(res_list)
  mat <- force_matrix_order(mat, ordered_types)
  if (!is_plottable_matrix(mat)) return(NULL)
  
  mat
}

# ---------------------------- Main loop over seu_split ----------------------------
top_ns <- list(1, 3, 5, 20, 50, "all")

for (top_n in top_ns) {
  
  cat("=== Processing top_n =", top_n, "===\n")
  
  outdir_top <- file.path(outdir_base, paste0("Top", top_n))
  dir.create(outdir_top, recursive = TRUE, showWarnings = FALSE)
  
  for (sample_name in names(seu_split)) {
    seu <- seu_split[[sample_name]]
    cat("Processing sample:", sample_name, "\n")
    
    # Determine ordered cell types for this Seurat
    ordered_types <- get_ordered_types(ref, seu, cell_type_groups)
    
    # ------------------ Similarity ------------------
    sim_dir <- file.path(outdir_top, "Similarity")
    dir.create(sim_dir, showWarnings = FALSE)
    
    sim_plots <- list()
    for (method in names(marker_dirs)) {
      marker_file <- file.path(marker_dirs[[method]], paste0(sample_name, "_cluster_markers.csv"))
      sim_mat <- compute_similarity_matrix(marker_file, ref, ordered_types, top_n)
      
      if (!is.null(sim_mat)) {
        ph <- pheatmap(
          sim_mat,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          display_numbers = TRUE,
          number_format = "%.1f%%",
          main = paste(sample_name, method, "Top", ifelse(top_n == "all", "all", as.character(top_n))),
          silent = TRUE
        )
        sim_plots[[method]] <- ph[[4]]
      }
    }
    
    if (length(sim_plots) > 0) {
      ncols <- ceiling(length(sim_plots)/2)
      out_png <- file.path(sim_dir, paste0(sample_name, "_similarity.png"))
      png(out_png, width = 8000, height = 3000, res = 200)
      gridExtra::grid.arrange(grobs = sim_plots, nrow = 2, ncol = ncols)
      dev.off()
      message("Saved similarity heatmap: ", out_png)
    }
    
    # ------------------ Fisher ------------------
    fisher_dir <- file.path(outdir_top, "Fisher")
    dir.create(fisher_dir, showWarnings = FALSE)
    
    fisher_plots <- list()
    for (method in names(marker_dirs)) {
      marker_file <- file.path(marker_dirs[[method]], paste0(sample_name, "_cluster_markers.csv"))
      mat <- compute_fisher_matrix(marker_file, ref, ordered_types, top_n)
      
      if (!is.null(mat)) {
        ph <- pheatmap(
          mat,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          display_numbers = TRUE,
          number_format = "%.1f",
          main = paste(sample_name, method, "Top", ifelse(top_n == "all", "all", as.character(top_n))),
          silent = TRUE
        )
        fisher_plots[[method]] <- ph[[4]]
      }
    }
    
    if (length(fisher_plots) > 0) {
      ncols <- ceiling(length(fisher_plots)/2)
      out_png <- file.path(fisher_dir, paste0(sample_name, "_fisher.png"))
      png(out_png, width = 8000, height = 3000, res = 200)
      gridExtra::grid.arrange(grobs = fisher_plots, nrow = 2, ncol = ncols)
      dev.off()
      message("Saved Fisher heatmap: ", out_png)
    }
    
  }
}
