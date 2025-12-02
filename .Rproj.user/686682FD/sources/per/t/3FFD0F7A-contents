print(getwd())
source("0_Functions.R")

# Load Seurat objects
seu_split_Method_1 <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/2.0 Method_1/2d Cell_labelled.Rds")
seu_split_Method_2 <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/3.0 Method_2/3c Cell_labelled.Rds")
seu_split_Method_3 <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/4.0 Method_3/4d Cell_labelled.Rds")
seu_split_Method_4 <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/5.0 Method_4/5c Cell_labelled.Rds")

seu_split_Draw_Breakers_1 <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/6.0 Cluster Methods/6c Cell_labelled_Full.Rds")
seu_split_Draw_Breakers_2 <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/6.0 Cluster Methods/6d Cell_labelled_Filtered.Rds")

seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/1.0 QC/sObj.Rds")


# ----------------------------  Copy Fisher assignments only ---------------------------------

# Function to copy over Fisher assignments
copy_fisher_assignments <- function(base_list, donor_list, new_label) {
  for (nm in names(base_list)) {
    if (nm %in% names(donor_list)) {
      donor_meta <- donor_list[[nm]]@meta.data
      
      if ("Fisher_assignment" %in% colnames(donor_meta)) {
        # Make a single column with custom label
        new_col <- donor_meta[rownames(base_list[[nm]]@meta.data), "Fisher_assignment", drop = FALSE]
        colnames(new_col) <- new_label
        
        # Add to metadata
        base_list[[nm]]@meta.data <- cbind(base_list[[nm]]@meta.data, new_col)
      } else {
        message(sprintf("[%s] No Fisher_assignment column found in donor list.", nm))
      }
    }
  }
  return(base_list)
}

# Apply for each method list
seu_split <- copy_fisher_assignments(seu_split, seu_split_Method_1, "Fisher_Method_1")
seu_split <- copy_fisher_assignments(seu_split, seu_split_Method_2, "Fisher_Method_2")
seu_split <- copy_fisher_assignments(seu_split, seu_split_Method_3, "Fisher_Method_3")
seu_split <- copy_fisher_assignments(seu_split, seu_split_Method_4, "Fisher_Method_4")
seu_split <- copy_fisher_assignments(seu_split, seu_split_Draw_Breakers_1, "Fisher_Draw_Breaker_1")
seu_split <- copy_fisher_assignments(seu_split, seu_split_Draw_Breakers_2, "Fisher_Draw_Breaker_2")

# Check result
head(seu_split$Wal40_P@meta.data[, grep("Fisher_", colnames(seu_split$Wal40_P@meta.data))])

# ----------------------------  Ensemble vote for final cell label ---------------------------------

# Directory for clashes
dir.create("Files/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7a_final_assigment_clashes", 
           recursive = TRUE, showWarnings = FALSE)

ensemble_vote <- function(seu_list, draw_breakers_1, draw_breakers_2) {
  clash_list <- list()
  
  for (nm in names(seu_list)) {
    meta <- seu_list[[nm]]@meta.data
    
    # Main assignment columns (methods 1–4 only)
    assign_cols <- c("Fisher_Method_1", "Fisher_Method_2",
                     "Fisher_Method_3", "Fisher_Method_4")
    assign_cols <- assign_cols[assign_cols %in% colnames(meta)]
    
    # Get breaker assignments if available
    breaker1 <- if (nm %in% names(draw_breakers_1)) {
      draw_breakers_1[[nm]]@meta.data$Fisher_assignment
    } else rep(NA, nrow(meta))
    
    breaker2 <- if (nm %in% names(draw_breakers_2)) {
      draw_breakers_2[[nm]]@meta.data$Fisher_assignment
    } else rep(NA, nrow(meta))
    
    final_labels <- character(nrow(meta))
    clashes <- list()
    
    for (i in seq_len(nrow(meta))) {
      votes <- meta[i, assign_cols, drop = FALSE]
      non_na_votes <- unlist(votes[!is.na(votes)])
      
      if (length(non_na_votes) == 0) {
        final_labels[i] <- "Unknown"
        next
      }
      
      vote_counts <- table(non_na_votes)
      max_count <- max(vote_counts)
      winners <- names(vote_counts[vote_counts == max_count])
      
      if (length(winners) == 1) {
        # Clear winner
        final_labels[i] <- winners
      } else {
        # Tie -> check breakers
        b1 <- breaker1[i]
        b2 <- breaker2[i]
        breaker_votes <- c(b1, b2)
        breaker_votes <- breaker_votes[!is.na(breaker_votes)]
        
        # See if either breaker supports one of the tied labels
        breaker_support <- breaker_votes[breaker_votes %in% winners]
        
        if (length(breaker_support) > 0) {
          # If both breakers agree on the same label → take that
          if (length(unique(breaker_support)) == 1) {
            final_labels[i] <- unique(breaker_support)
          } else {
            # If they disagree → still a tie
            clash_label <- paste(sort(winners), collapse = ":")
            final_labels[i] <- clash_label
          }
        } else {
          # No help from breakers → keep tie
          clash_label <- paste(sort(winners), collapse = ":")
          final_labels[i] <- clash_label
        }
        
        # Record the clash info
        clash_entry <- data.frame(
          Seurat = nm,
          Cell = rownames(meta)[i],
          votes,
          Breaker1 = b1,
          Breaker2 = b2,
          Tied_Labels = paste(sort(winners), collapse = ":"),
          Final_Label = final_labels[i],
          stringsAsFactors = FALSE
        )
        clashes[[length(clashes) + 1]] <- clash_entry
      }
    }
    
    # Add final label to Seurat metadata
    seu_list[[nm]]$ensemble_cell_label <- final_labels
    
    # Save clashes
    if (length(clashes) > 0) {
      clash_df <- do.call(rbind, clashes)
      clash_list[[nm]] <- clash_df
      clash_path <- file.path("Files/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7a_final_assigment_clashes",
                              paste0(nm, "_clashes.csv"))
      write.csv(clash_df, clash_path, row.names = FALSE)
    }
  }
  
  return(seu_list)
}

# Run with your breaker lists
seu_split <- ensemble_vote(seu_split, seu_split_Draw_Breakers_1, seu_split_Draw_Breakers_2)

# Check result
head(seu_split$Wal40_P@meta.data)

table(seu_split$Imp2_P@meta.data$ensemble_cell_label)

# ----------------------------  Add ensemble_label_w/o_clash column ---------------------------------
for (nm in names(seu_split)) {
  meta <- seu_split[[nm]]@meta.data
  
  if ("ensemble_cell_label" %in% colnames(meta)) {
    new_col <- as.character(meta$ensemble_cell_label)
    # Any label with ":" is a clash → set to Unknown
    new_col[grepl(":", new_col)] <- "Unknown"
    seu_split[[nm]]$ensemble_label_w_o_clash <- new_col
  } else {
    message(sprintf("[%s] No ensemble_cell_label column found.", nm))
  }
}

# ----------------------------  Normalise and Scale ---------------------------------

Normalise <- function(sObj, sample_name){
  sObj <- NormalizeData(sObj)
  sObj <- ScaleData(object = sObj, features = rownames(sObj))
  sObj <- FindVariableFeatures(sObj)
  sObj <- RunPCA(object = sObj, npcs = 40)
  
  sObj <- FindNeighbors(object = sObj, dims = 1:40,  reduction = "pca")
  resolutions <- 2.0
  for (res in resolutions) {
    sObj <- RunUMAP(object = sObj, dims = 1:40)
  }
  
  return(sObj)
}

# Apply the function to each Seurat object in seu_split
for (sample_name in names(seu_split)) {
  cat("Normalising:", sample_name, "\n")
  seu_split[[sample_name]] <- Normalise(seu_split[[sample_name]], sample_name)
}

# ----------------------------  UMAP panels for assignment columns ---------------------------------

# Define color scheme
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

outdir_umaps <- "Figures/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7a_Ensemble_Votes/UMAPs"
dir.create(outdir_umaps, recursive = TRUE, showWarnings = FALSE)

make_umap_panels <- function(seu_list,
                             cols = c("Fisher_Method_1","Fisher_Method_2","Fisher_Method_3",
                                      "Fisher_Method_4", "Fisher_Draw_Breaker_1", "Fisher_Draw_Breaker_2",
                                      "ensemble_label_w_o_clash"),
                             outdir = outdir_umaps,
                             max_cols = 4,
                             point_size = 0.5) {
  for (nm in names(seu_list)) {
    obj <- seu_list[[nm]]
    
    # Ensure UMAP exists
    if (!"umap" %in% names(obj@reductions)) {
      if ("pca" %in% names(obj@reductions)) {
        message(sprintf("[%s] UMAP missing — computing from PCA...", nm))
        obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
        seu_list[[nm]] <- obj
      } else {
        message(sprintf("[%s] Skipping: no UMAP or PCA present.", nm))
        next
      }
    }
    
    # Build per-column UMAPs
    plot_list <- list()
    present_cols <- cols[cols %in% colnames(obj@meta.data)]
    if (length(present_cols) == 0) {
      message(sprintf("[%s] No requested columns found. Skipping.", nm))
      next
    }
    
    # Clean NAs -> "Unknown"
    for (col in present_cols) {
      vec <- obj@meta.data[[col]]
      vec <- as.character(vec)
      vec[is.na(vec) | vec == ""] <- "Unknown"
      obj@meta.data[[col]] <- factor(vec)
    }
    
    # Generate plots without legends
    for (col in present_cols) {
      p <- DimPlot(
        obj,
        reduction = "umap",
        group.by = col,
        raster = FALSE,
        pt.size = point_size
      ) +
        scale_color_manual(values = celltype_colors, drop = FALSE) +   # enforce palette
        ggtitle(col) +
        theme_bw(base_size = 11) +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid = element_blank()
        )
      plot_list[[col]] <- p
    }
    
    # One shared legend
    legend_col <- present_cols[1]
    p_leg <- DimPlot(
      obj,
      reduction = "umap",
      group.by = legend_col,
      raster = FALSE,
      pt.size = point_size
    ) +
      scale_color_manual(values = celltype_colors, drop = FALSE) +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
      theme_bw(base_size = 11) +
      theme(legend.position = "right", panel.grid = element_blank())
    shared_legend <- cowplot::get_legend(p_leg)
    
    # Arrange plots into grid
    ncol_plots <- min(max_cols, length(plot_list))
    nrow_plots <- ceiling(length(plot_list) / ncol_plots)
    grid <- wrap_plots(plot_list, ncol = ncol_plots) +
      plot_annotation(title = paste0(nm, " — Ensemble Vote UMAPs"))
    
    final_plot <- cowplot::plot_grid(grid, shared_legend, rel_widths = c(1, 0.24))
    
    # Save
    width <- 4 * ncol_plots + 3
    height <- 4 * nrow_plots
    outfile <- file.path(outdir, paste0(nm, "_Ensemble_Votes_UMAPs.png"))
    ggsave(outfile, final_plot, width = width, height = height, dpi = 300)
    message(sprintf("[%s] Saved: %s", nm, outfile))
  }
  
  invisible(seu_list)
}

# Run it on your updated list
seu_split <- make_umap_panels(seu_split)


# ----------------------------  Save Seurats ---------------------------------
dir.create("Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/", 
           showWarnings = FALSE, 
           recursive = TRUE)

saveRDS(seu_split, "Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7a Ensemble_cell_labelled.Rds")

