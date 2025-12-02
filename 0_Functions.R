setwd("/Users/bs23jmw/Library/CloudStorage/OneDrive-UniversityofLeeds/Shoaib Analysis/Cell typing")

library(Seurat)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(patchwork)
library(cowplot)
library(tidyr)
library(stringr)
library(tidyverse)
library(readxl)
library(openxlsx)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(qpdf)
library(readr)
library(gridExtra)
library(grid)
library(ggalluvial)
library(ComplexUpset)

### Functions
howManyCells <- function(sObj){dim(sObj)[2]}




