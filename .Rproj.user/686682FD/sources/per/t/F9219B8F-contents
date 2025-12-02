print(getwd())
source("0_Functions.R")

seu_split <- readRDS("Seurat Objects/Transfer Report Pipeline Mk3/7.0 Ensemble_Vote/7b Clash_cell_labelled.Rds")

# ----------------------------  Which assignment is used most for clashes ---------------------------------

check_resolutions <- function(seu){
  meta <- seu@meta.data
  
  # Only look at cells with clashes (":" in ensemble_cell_label)
  clash_cells <- meta[grepl(":", meta$ensemble_cell_label), ]
  
  # Remove clashes that were resolved as "Unknown"
  clash_cells <- clash_cells[clash_cells$Resolved_Label != "Unknown", ]
  
  if(nrow(clash_cells) == 0) return(NULL)
  
  method_cols <- c("Fisher_Method_1","Fisher_Method_2","Fisher_Method_3",
                   "Fisher_Method_4", "Fisher_Draw_Breaker_1", "Fisher_Draw_Breaker_2")
  
  # Accuracy per method vs Resolved_Label
  results <- sapply(method_cols, function(col){
    mean(clash_cells[[col]] == clash_cells$Resolved_Label)
  })
  
  results <- c(results, n_clashes = nrow(clash_cells))
  return(results)
}

resolution_stats <- lapply(seu_split, check_resolutions)
resolution_stats_df <- do.call(rbind, resolution_stats) %>% as.data.frame()
resolution_stats_df$Sample <- rownames(resolution_stats_df)

resolution_stats_df <- resolution_stats_df[, c("Sample", 
                                               "Fisher_Method_1","Fisher_Method_2","Fisher_Method_3",
                                               "Fisher_Method_4", "Fisher_Draw_Breaker_1", "Fisher_Draw_Breaker_2",
                                               "n_clashes")]
head(resolution_stats_df)

# ----------------------------  Visualise clash accuracy ---------------------------------

df_long <- resolution_stats_df %>%
  pivot_longer(cols = c(starts_with("Fisher_Method"), starts_with("Fisher_Draw_Breaker")),
               names_to = "method", values_to = "accuracy") %>%
  filter(!is.na(accuracy))

method_type_stats <- df_long %>%
  group_by(method) %>%
  summarise(mean_accuracy = mean(accuracy, na.rm = TRUE),
            median_accuracy = median(accuracy, na.rm = TRUE),
            sd_accuracy = sd(accuracy, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(mean_accuracy))

outdir <- "Figures/Transfer Report Pipeline Mk3/8.0 Method_Evaluation/8a_Clash_Method_Comparison"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

p2 <- ggplot(method_type_stats, aes(x = reorder(method, mean_accuracy), 
                                    y = mean_accuracy, fill = method)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_accuracy - sd_accuracy, 
                    ymax = mean_accuracy + sd_accuracy), 
                width = 0.2, colour = "black") +
  geom_text(aes(label = round(mean_accuracy, 3)), hjust = -0.1, size = 3.5) +
  coord_flip(ylim = c(0,1)) +
  labs(title = "Clash Accuracy by Fisher Method",
       x = "Method", y = "Mean Accuracy ± SD") +
  theme_minimal(base_size = 14) + 
  theme(plot.background = element_rect(fill = 'white'),
        legend.position = "none")

ggsave(file.path(outdir, "Clash_Method_Comparison.png"), p2, width = 6, height = 5, dpi = 300)

# ----------------------------  Which assignment is used most (all cells) ---------------------------------

check_resolutions_all <- function(seu){
  meta <- seu@meta.data
  all_cells <- meta[meta$Resolved_Label != "Unknown", ]
  if(nrow(all_cells) == 0) return(NULL)
  
  method_cols <- c("Fisher_Method_1","Fisher_Method_2","Fisher_Method_3",
                   "Fisher_Method_4", "Fisher_Draw_Breaker_1", "Fisher_Draw_Breaker_2")
  
  results <- sapply(method_cols, function(col){
    mean(all_cells[[col]] == all_cells$Resolved_Label)
  })
  
  results <- c(results, n_cells = nrow(all_cells))
  return(results)
}

resolution_stats_all <- lapply(seu_split, check_resolutions_all)
resolution_stats_df_all <- do.call(rbind, resolution_stats_all) %>% as.data.frame()
resolution_stats_df_all$Sample <- rownames(resolution_stats_df_all)

resolution_stats_df_all <- resolution_stats_df_all[, c("Sample", 
                                                       "Fisher_Method_1","Fisher_Method_2","Fisher_Method_3",
                                                       "Fisher_Method_4", "Fisher_Draw_Breaker_1", "Fisher_Draw_Breaker_2",
                                                       "n_cells")]
head(resolution_stats_df_all)

# ----------------------------  Visualise all cells accuracy ---------------------------------

df_long_all <- resolution_stats_df_all %>%
  pivot_longer(cols = c(starts_with("Fisher_Method"), starts_with("Fisher_Draw_Breaker")),
               names_to = "method", values_to = "accuracy") %>%
  filter(!is.na(accuracy))

method_type_stats_all <- df_long_all %>%
  group_by(method) %>%
  summarise(mean_accuracy = mean(accuracy, na.rm = TRUE),
            median_accuracy = median(accuracy, na.rm = TRUE),
            sd_accuracy = sd(accuracy, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(mean_accuracy))

p4 <- ggplot(method_type_stats_all, aes(x = reorder(method, mean_accuracy), 
                                        y = mean_accuracy, fill = method)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_accuracy - sd_accuracy, 
                    ymax = mean_accuracy + sd_accuracy), 
                width = 0.2, colour = "black") +
  geom_text(aes(label = round(mean_accuracy, 3)), hjust = -0.1, size = 3.5) +
  coord_flip(ylim = c(0,1)) +
  labs(title = "All Cells 'Wins' by Fisher Method",
       x = "Method", y = "Mean 'Accuracy 'Wins' ± SD") +
  theme_minimal(base_size = 14) + 
  theme(plot.background = element_rect(fill = 'white'),
        legend.position = "none")

ggsave(file.path(outdir, "AllCells_Method_Comparison.png"), p4, width = 6, height = 5, dpi = 300)
