print(getwd())
source("0_Functions.R")

seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7a Ensemble_cell_labelled.Rds")
ref <- read.csv("Files/Ref/Shoaib_ref.csv")


# --------------------  Resolve the clashes (otherwise "Unknown") ----------------------------

resolve_clashes <- function(seu, label_col = "ensemble_cell_label", ref, new_col = "Resolved_Label") {
  # Copy existing labels into new column
  seu@meta.data[[new_col]] <- seu@meta.data[[label_col]]
  
  # Find clash cells (those with ":" in label)
  clash_cells <- rownames(seu@meta.data)[grep(":", seu@meta.data[[label_col]])]
  if (length(clash_cells) == 0) {
    message("No clash cells found in ", label_col)
    return(seu)
  }
  
  message("Found ", length(clash_cells), " clash cells to resolve...")
  
  # Iterate through clash cells
  for (cell in clash_cells) {
    clash_label <- seu@meta.data[cell, label_col]
    candidates <- strsplit(clash_label, ":")[[1]]
    
    # Get expression profile for this cell
    expr <- GetAssayData(seu, slot = "data")[, cell]
    
    # Score each candidate using mean of expressed markers only
    scores <- sapply(candidates, function(ct) {
      markers <- ref$Gene[ref$Cell_Type == ct]
      markers <- markers[markers %in% rownames(seu)]
      markers <- markers[expr[markers] > 0]  # only expressed markers
      if (length(markers) == 0) return(NA)
      mean(expr[markers], na.rm = TRUE)
    })
    
    # If no informative markers, assign "Unknown"
    if (all(is.na(scores))) {
      best_ct <- "Unknown"
    } else {
      best_ct <- names(which.max(scores))
    }
    
    # Assign back to metadata
    seu@meta.data[cell, new_col] <- best_ct
  }
  
  return(seu)
}

# Loop through each sample in your split list
seu_split <- lapply(seu_split, function(x) {
  resolve_clashes(x, label_col = "ensemble_cell_label", ref = ref)
})


# --------------------  Resolve the clashes (otherwise retain clash) ----------------------------

# Function (no change needed inside resolve_clashes except sample_name handling)
resolve_clashes <- function(seu, label_col = "ensemble_cell_label", ref, 
                            new_col = "Resolved_Label_with_clash", 
                            out_dir = "Files/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7b_marker_resolution",
                            sample_name = NULL) {
  # Make sure output dir exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if (is.null(sample_name)) {
    stop("sample_name must be provided when running inside lapply")
  }
  
  # Copy existing labels into new column
  seu@meta.data[[new_col]] <- seu@meta.data[[label_col]]
  
  # Find clash cells (those with ":" in label)
  clash_cells <- rownames(seu@meta.data)[grep(":", seu@meta.data[[label_col]])]
  if (length(clash_cells) == 0) {
    message("No clash cells found in ", sample_name)
    return(seu)
  }
  
  message("Found ", length(clash_cells), " clash cells to resolve in ", sample_name, "...")
  
  # Store resolution details
  resolution_log <- list()
  
  for (cell in clash_cells) {
    clash_label <- seu@meta.data[cell, label_col]
    candidates <- strsplit(clash_label, ":")[[1]]
    
    expr <- GetAssayData(seu, slot = "data")[, cell]
    
    candidate_scores <- list()
    
    for (ct in candidates) {
      markers <- ref$Gene[ref$Cell_Type == ct]
      markers <- markers[markers %in% rownames(seu)]
      
      if (length(markers) > 0) {
        exp_vals <- expr[markers]
        expressed <- exp_vals[exp_vals > 0]
        score <- if (length(expressed) > 0) mean(expressed) else NA
        
        candidate_scores[[ct]] <- list(
          markers = names(expressed),
          values = expressed,
          score = score
        )
      } else {
        candidate_scores[[ct]] <- list(
          markers = character(0),
          values = numeric(0),
          score = NA
        )
      }
    }
    
    # Decide winner
    scores <- sapply(candidate_scores, function(x) x$score)
    if (all(is.na(scores))) {
      best_ct <- clash_label
    } else {
      best_ct <- names(which.max(scores))
    }
    
    seu@meta.data[cell, new_col] <- best_ct
    
    # Expand log into rows
    for (ct in names(candidate_scores)) {
      if (length(candidate_scores[[ct]]$markers) > 0) {
        tmp <- data.frame(
          cell = cell,
          clash = clash_label,
          candidate = ct,
          marker = candidate_scores[[ct]]$markers,
          expression = as.numeric(candidate_scores[[ct]]$values),
          score = candidate_scores[[ct]]$score,
          final_label = best_ct,
          stringsAsFactors = FALSE
        )
      } else {
        tmp <- data.frame(
          cell = cell,
          clash = clash_label,
          candidate = ct,
          marker = NA,
          expression = NA,
          score = candidate_scores[[ct]]$score,
          final_label = best_ct,
          stringsAsFactors = FALSE
        )
      }
      resolution_log[[length(resolution_log) + 1]] <- tmp
    }
  }
  
  # Combine and save log
  resolution_df <- do.call(rbind, resolution_log)
  out_file <- file.path(out_dir, paste0(sample_name, "_marker_resolution.csv"))
  write.csv(resolution_df, out_file, row.names = FALSE)
  
  message("Resolution log saved to: ", out_file)
  
  return(seu)
}

# Loop through each sample in your split list, passing the name
seu_split <- mapply(function(x, nm) {
  resolve_clashes(x, label_col = "ensemble_cell_label", ref = ref, sample_name = nm)
}, seu_split, names(seu_split), SIMPLIFY = FALSE)


# ----------------------------  Save Seurats ---------------------------------
saveRDS(seu_split, "Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7b Clash_cell_labelled.Rds")
