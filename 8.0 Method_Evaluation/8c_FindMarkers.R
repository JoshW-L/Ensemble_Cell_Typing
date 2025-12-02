print(getwd())
source("0_Functions.R")

seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7b Clash_cell_labelled.Rds")

# ----------------------------  Find Markers ---------------------------------

# Define the label columns and output directories
methods <- list(
  Fisher_Method_1 = "Method_1",
  Fisher_Method_2 = "Method_2",
  Fisher_Method_3 = "Method_3",
  Fisher_Method_4 = "Method_4",
  Fisher_Draw_Breaker_1 = "Draw_Breaker_1",
  Fisher_Draw_Breaker_2 = "Draw_Breaker_2",
  ensemble_label_w_o_clash = "Ensemble_Vote",
  Resolved_Label = "Resolved")

base_dir <- "Files/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8c_Markers_By_Patients"

# Loop over each method/label
for (label_col in names(methods)) {
  outdir <- file.path(base_dir, methods[[label_col]])
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n=== Processing", label_col, "===\n")
  
  for (obj_name in names(seu_split)) {
    cat("  -", obj_name, "\n")
    obj <- seu_split[[obj_name]]
    
    # Set cluster identities to the chosen column
    Idents(obj) <- obj[[label_col, drop = TRUE]]
    
    # Run FindAllMarkers
    markers <- FindAllMarkers(
      obj,
      only.pos = TRUE,
      logfc.threshold = 0.1,
      min.pct = 0.1,
      test.use = "wilcox",
      assay = DefaultAssay(obj),
      p.adjust.method = "bonferroni"
    ) %>%
      dplyr::filter(p_val_adj < 0.1)
    
    # Save results
    outpath <- file.path(outdir, paste0(obj_name, "_cluster_markers.csv"))
    write.csv(markers, outpath, row.names = FALSE)
  }
}
