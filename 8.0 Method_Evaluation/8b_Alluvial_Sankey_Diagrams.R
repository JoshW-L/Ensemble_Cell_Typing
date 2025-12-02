print(getwd())
source("0_Functions.R")

seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7b Clash_cell_labelled.Rds")

# Directories
dir_alluvial <- "Figures/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8b_Alluvial_UpSet_Plots/Alluvial"
dir_upset    <- "Figures/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8b_Alluvial_UpSet_Plots/UpSet"
dir_overall  <- "Figures/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8b_Alluvial_UpSet_Plots/Merged"

dir.create(dir_alluvial, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_upset, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_overall, recursive = TRUE, showWarnings = FALSE)

# Use the actual Fisher method columns
methods <- c("Fisher_Method_1","Fisher_Method_2","Fisher_Method_3",
             "Fisher_Method_4", "Fisher_Draw_Breaker_1", "Fisher_Draw_Breaker_2")

# ------------------ Alluvial ------------------
plot_alluvial <- function(seu, sample_name) {
  df <- seu@meta.data[, methods] %>%
    as.data.frame() %>%
    mutate(Cell = rownames(seu@meta.data))
  
  df_long <- df %>%
    pivot_longer(cols = all_of(methods), names_to = "Method", values_to = "Label")
  
  p <- ggplot(df_long,
              aes(x = Method, stratum = Label, alluvium = Cell, fill = Label, label = Label)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.5) +
    geom_stratum() +
    theme_minimal(base_size = 14) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA)
    ) +
    ggtitle(paste("Alluvial: ", sample_name))
  
  ggsave(file.path(dir_alluvial, paste0(sample_name, "_Alluvial.png")), p, width = 15, height = 6, dpi = 300, bg = "white")
  return(p)
}

# ------------------ UpSet ------------------
plot_upset <- function(seu, sample_name) {
  df <- seu@meta.data[, methods] %>% as.data.frame()
  
  # Reference: resolved label
  ref <- seu@meta.data$Resolved_Label
  df_bin <- sapply(df, function(x) x == ref) %>% as.data.frame()
  
  # ComplexUpset requires df + at least one numeric/categorical column
  df_bin$Cell <- rownames(seu@meta.data)  # keep ID for safety
  
  # Plot
  p <- ComplexUpset::upset(
    df_bin,
    intersect = methods  # instead of sets=
  ) + ggtitle(paste("UpSet: ", sample_name))
  
  ggsave(file.path(dir_upset, paste0(sample_name, "_UpSet.png")), p, width = 20, height = 6, dpi = 300)
  return(p)
}

# ------------------ Apply per Seurat ------------------
for(nm in names(seu_split)){
  message("Processing ", nm)
  plot_alluvial(seu_split[[nm]], nm)
  plot_upset(seu_split[[nm]], nm)
}


# ------------------ Overall ------------------

# ----------------------------  Merge Seurats ---------------------------------

# Merge all Seurat objects in the list into one
seu_merged <- Reduce(function(x, y) merge(x, y), seu_split)
seu_merged <- JoinLayers(seu_merged)

seu_merged <- NormalizeData(seu_merged, assay = "Nanostring", verbose = TRUE)
seu_merged <- FindVariableFeatures(seu_merged)
seu_merged <- ScaleData(seu_merged)
seu_merged <- RunPCA(seu_merged, npcs = 50)
seu_merged <- RunUMAP(seu_merged, dims = 1:50)




# Use the actual Fisher method columns
methods <- c("Fisher_Method_1","Fisher_Method_2","Fisher_Method_3",
             "Fisher_Method_4", "Fisher_Draw_Breaker_1", "Fisher_Draw_Breaker_2")

# ------------------ Alluvial ------------------
plot_alluvial <- function(seu, sample_name = "Overall") {
  df <- seu@meta.data[, methods] %>%
    as.data.frame() %>%
    mutate(Cell = rownames(seu@meta.data))
  
  df_long <- df %>%
    pivot_longer(cols = all_of(methods), names_to = "Method", values_to = "Label")
  
  p <- ggplot(df_long,
              aes(x = Method, stratum = Label, alluvium = Cell, fill = Label, label = Label)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.5) +
    geom_stratum() +
    theme_minimal(base_size = 14) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA)
    ) +
    ggtitle(paste("Alluvial: ", sample_name))
  
  ggsave(file.path(dir_overall, paste0(sample_name, "_Alluvial.png")), 
         p, width = 15, height = 6, dpi = 300, bg = "white")
  return(p)
}

# ------------------ UpSet ------------------
plot_upset <- function(seu, sample_name = "Overall") {
  df <- seu@meta.data[, methods] %>% as.data.frame()
  
  # Reference: resolved label
  ref <- seu@meta.data$Resolved_Label
  df_bin <- sapply(df, function(x) x == ref) %>% as.data.frame()
  
  # ComplexUpset requires df + at least one numeric/categorical column
  df_bin$Cell <- rownames(seu@meta.data)
  
  # Plot
  p <- ComplexUpset::upset(
    df_bin,
    intersect = methods
  ) + ggtitle(paste("UpSet: ", sample_name))
  
  ggsave(file.path(dir_overall, paste0(sample_name, "_UpSet.png")), 
         p, width = 25, height = 6, dpi = 300)
  return(p)
}

# ------------------ Run once for overall ------------------
message("Processing overall")
plot_alluvial(seu_merged, "Overall")
plot_upset(seu_merged, "Overall")

# ----------------------------  Save Merged Seurat ---------------------------------
dir.create("Seurat Objects/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/", 
           showWarnings = FALSE, 
           recursive = TRUE)
saveRDS(seu_merged, "Seurat Objects/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8b Merged_Seurat.Rds")
