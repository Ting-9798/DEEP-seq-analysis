## Visualization
rm(list=ls())


############## GSVA Analysis ########
library(ggplot2)
library(dplyr)
library(msigdbr)
library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)

#BiocManager::install("GSVA")

setwd("/out/")
outdir <- "/out/"


## Output directory
output <- file.path(outdir, "GSVA")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

## Read object
file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)
summary(scRNAsub$Cancer_3class)
view(scRNAsub@meta.data)

#################### 1. Prepare expression matrix and meta ####################
# expr: genes x cells
expr <- as.data.frame(scRNAsub@assays$RNA@data)
expr <- expr[rowSums(expr)>0,]  # select non-zero genes

# meta: required grouping information, assuming meta.data already contains celltype / stage / treatment
meta <- scRNAsub@meta.data[, c("Metastatic_3class", "Digestive_tract_3class","Ovarian_3class","Cancer_3class","treatment")]


#################### 2. Extract KEGG gene sets and perform GSVA ####################
# Get all C2 gene sets first, then filter KEGG-related subsets
m_all_c2 <- msigdbr(
  species    = "Homo sapiens",
  collection = "C2"               # replaces original category = "C2"
)

# Keep KEGG-related subsets: CP:KEGG_LEGACY and CP:KEGG_MEDICUS
m_df <- dplyr::filter(
  m_all_c2,
  gs_subcollection %in% c("CP:KEGG_LEGACY")
)

msigdbr_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)

expr_mat <- as.matrix(expr)
# Note: For UMI-data, strictly speaking kcdf="Poisson" is more common, keeping Gaussian as per example
# 1) Build parameter object
gsva_par <- gsvaParam(
  exprData = expr_mat,
  geneSets = msigdbr_list,
  kcdf     = "Gaussian"
  # or "Poisson" depending on your data type
  # Other parameters like minSize / maxSize / tau can also be set here
)

# 2) Parallel parameters (recommended MulticoreParam for Linux/server)
bp <- MulticoreParam(workers = 20)

# 3) Run GSVA
kegg <- gsva(gsva_par, verbose = TRUE)

rownames(kegg) <- gsub("^KEGG_", "", rownames(kegg))

write.csv(
  kegg,
  file = file.path(output, "GSVA_scores.csv"),
  quote = FALSE
)




#################### 3. Average GSVA scores by treatment + celltype combination groups ####################
# Create combined grouping
meta$combined_group <- paste0(meta$treatment, "_", meta$celltype)
cluster_vec <- meta$combined_group
names(cluster_vec) <- rownames(meta)

# Ensure column names correspond to meta row names
kegg11 <- kegg[, names(cluster_vec), drop = FALSE]

# Calculate mean score for each combined group
kegg_combined <- sapply(
  X = unique(cluster_vec),
  FUN = function(cl){
    cells <- names(cluster_vec)[cluster_vec == cl]
    rowMeans(kegg11[, cells, drop = FALSE])
  }
)

# Convert to matrix, rows are pathways, columns are combined groups
kegg_combined <- as.matrix(kegg_combined)

#################### 4. Reorder columns by celltype, then by treatment within each celltype ####################
# Extract all unique celltypes and treatments
all_celltypes <- sort(unique(meta$celltype))
all_treatments <- sort(unique(meta$treatment))

# Create new column order: first by celltype, then by treatment within each celltype
new_col_order <- c()
for(celltype in all_celltypes) {
  for(treat in all_treatments) {
    group_name <- paste0(treat, "_", celltype)
    if(group_name %in% colnames(kegg_combined)) {
      new_col_order <- c(new_col_order, group_name)
    }
  }
}

# Reorder matrix columns
kegg_combined <- kegg_combined[, new_col_order, drop = FALSE]

#################### 5. Select top 50 most variable pathways ####################
pathway_var <- apply(kegg_combined, 1, var, na.rm = TRUE)
pathway_var <- sort(pathway_var, decreasing = TRUE)

top50_pathways <- names(pathway_var)[1:6]
kegg_top50 <- kegg_combined[top50_pathways, , drop = FALSE]



# Define list of pathways of interest, remove "KEGG_" prefix
kegg_interest_ids <- c(
  "JAK_STAT_SIGNALING_PATHWAY",
  "PI3K_AKT_SIGNALING_PATHWAY", "ABC_TRANSPORTERS", "GLUTATHIONE_METABOLISM", 
  "METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450", "DRUG_METABOLISM_OTHER_ENZYMES", 
  "PATHWAYS_IN_CANCER", "MAPK_SIGNALING_PATHWAY", 
  "APOPTOSIS", "MTOR_SIGNALING_PATHWAY", "FOCAL_ADHESION", 
  "TGF_BETA_SIGNALING_PATHWAY", "HEDGEHOG_SIGNALING_PATHWAY", "WNT_SIGNALING_PATHWAY", 
  "CELL_ADHESION_MOLECULES_CAMS", "GLYCOLYSIS_GLUCONEOGENESIS","TIGHT_JUNCTION","OXIDATIVE_PHOSPHORYLATION"
  
)

# Find pathways that actually exist in the matrix
existing_pathways <- kegg_interest_ids[kegg_interest_ids %in% rownames(kegg_combined)]

# Output pathways not found
kegg_interest <- kegg_combined[existing_pathways, , drop = FALSE]






#################### 6. Prepare heatmap annotation information ####################
# Extract grouping information for column annotation
col_annot <- data.frame(
  treatment = factor(gsub("_.*", "", colnames(kegg_top50)), levels = all_treatments),
  celltype = factor(gsub(".*_", "", colnames(kegg_top50)), levels = all_celltypes)
)
rownames(col_annot) <- colnames(kegg_top50)

# Set color scheme
library(RColorBrewer)

# Generate colors for treatment - use Set2 palette
treatment_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(all_treatments))
names(treatment_colors) <- all_treatments

# Generate colors for celltype - use Set3 palette (provides more colors)
celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_celltypes))
names(celltype_colors) <- all_celltypes

# Create annotation color list
annotation_colors <- list(
  treatment = treatment_colors,
  celltype = celltype_colors
)

#################### 7. Calculate differences between treatments within each celltype ####################
# Optional: Calculate difference matrix between different treatments within each celltype
# This helps identify which pathways are most sensitive to treatment changes in specific celltypes

celltype_treatment_differences <- list()

for(celltype in all_celltypes) {
  # Extract all columns for current celltype
  celltype_cols <- grep(paste0("_", celltype, "$"), colnames(kegg_interest), value = TRUE)
  
  if(length(celltype_cols) > 1) {
    # Calculate average difference between treatments within current celltype
    celltype_data <- kegg_top50[, celltype_cols, drop = FALSE]
    
    # If there are multiple treatments, calculate pairwise average absolute differences
    if(ncol(celltype_data) > 1) {
      # Calculate average absolute difference across all treatment combinations for all pathways
      diff_matrix <- matrix(0, nrow = nrow(celltype_data), ncol = ncol(celltype_data),
                            dimnames = list(rownames(celltype_data), colnames(celltype_data)))
      
      for(i in 1:ncol(celltype_data)) {
        for(j in 1:ncol(celltype_data)) {
          if(i != j) {
            diff_matrix[, i] <- diff_matrix[, i] + abs(celltype_data[, i] - celltype_data[, j])
          }
        }
        diff_matrix[, i] <- diff_matrix[, i] / (ncol(celltype_data) - 1)
      }
      
      # Calculate average difference for each pathway
      pathway_avg_diff <- rowMeans(diff_matrix)
      pathway_avg_diff <- sort(pathway_avg_diff, decreasing = TRUE)
      
      celltype_treatment_differences[[celltype]] <- list(
        data = celltype_data,
        diff_matrix = diff_matrix,
        pathway_avg_diff = pathway_avg_diff
      )
    }
  }
}


#################### 7. Draw heatmap (display top 50 pathways) ####################
library(ComplexHeatmap)
library(circlize)

# Normalize data (row-wise scaling)
kegg_top50_scaled <- t(scale(t(kegg_top50)))

# Determine grouping separator positions - add separator lines where celltype changes
celltype_labels <- col_annot$celltype
gaps_positions <- c()
for(i in 1:(length(celltype_labels)-1)) {
  if(celltype_labels[i] != celltype_labels[i+1]) {
    gaps_positions <- c(gaps_positions, i)
  }
}

# Create column annotation
col_ha <- HeatmapAnnotation(
  df = col_annot,
  col = annotation_colors,
  annotation_name_side = "left",
  gap = unit(2, "mm")
)

# Create row annotation (optional, display pathway classification)
# row_ha <- rowAnnotation(
#   pathway = rownames(kegg_top50_scaled),
#   show_annotation_name = FALSE
# )

# Draw heatmap
ht <- Heatmap(
  kegg_top50_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  
  # Row settings
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 12),
  row_names_max_width = unit(10, "cm"),
  row_title = "KEGG Pathways",
  
  # Column settings
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = col_annot$celltype,
  column_gap = unit(3, "mm"),
  column_title = "Cell Types",
  column_title_side = "top",
  
  # Annotations
  top_annotation = col_ha,
  
  # Heatmap parameters
  border = TRUE,
  border_gp = gpar(col = "gray", lty = 1),
  heatmap_legend_param = list(
    title = "Z-score",
    title_position = "leftcenter-rot"
  ),
  
  # Row clustering parameters
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  # Display parameters
  show_heatmap_legend = TRUE,
  use_raster = ifelse(ncol(kegg_top50_scaled) * nrow(kegg_top50_scaled) > 20000, TRUE, FALSE)
)

# Draw and save
pdf(file.path(output, "GSVA_KEGG_top50_by_celltype_treatment.pdf"), 
    width = max(10, ncol(kegg_top50) * 0.4 + 4),
    height = max(8, nrow(kegg_top50) * 0.15 + 3))

draw(ht, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)

dev.off()


#################### 7. Draw heatmap (display top 50 pathways) - optimized version ####################
library(ComplexHeatmap)
library(circlize)

# Normalize data (row-wise scaling)
kegg_top50_scaled <- t(scale(t(kegg_top50)))

# Create color mapping function
col_fun <- colorRamp2(
  seq(min(kegg_top50_scaled, na.rm = TRUE), 
      max(kegg_top50_scaled, na.rm = TRUE), 
      length = 100),
  colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)

# Create column annotation - optimized legend parameters
col_ha <- HeatmapAnnotation(
  df = col_annot,
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_legend_param = list(
    # Set legend to display in two rows
    nrow = 2,
    ncol = ceiling(length(annotation_colors) / 2),
    direction = "horizontal",
    
    # Increase legend title size
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    # Increase legend label size
    labels_gp = gpar(fontsize = 12),
    
    # Increase legend graphic size
    legend_height = unit(1.2, "cm"),
    legend_width = unit(1.2, "cm"),
    
    # Set spacing between legends
    gap = unit(0.3, "cm"),
    
    title_position = "topcenter"
  ),
  gap = unit(2, "mm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12)
)

# Draw heatmap
ht <- Heatmap(
  kegg_top50_scaled,
  name = "Z-score",
  col = col_fun,
  
  # Row settings
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_names_max_width = unit(14, "cm"),
  row_title = "KEGG Pathways",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # Column settings
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = col_annot$celltype,
  column_gap = unit(3, "mm"),
  column_title = "Cell Types",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # Annotations
  top_annotation = col_ha,
  
  # Heatmap parameters
  border = TRUE,
  border_gp = gpar(col = "gray", lty = 1),
  
  # Heatmap legend settings - optimized version
  heatmap_legend_param = list(
    title = "GSVA Z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(10, "cm"),
    legend_height = unit(1.5, "cm"),
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    at = round(seq(min(kegg_top50_scaled), max(kegg_top50_scaled), length.out = 5), 1),
    grid_border = "black",
    color_bar = "continuous",
    labels_rot = 0,
    border = "black",
    padding = unit(c(8, 8, 8, 8), "mm")
  ),
  
  # Row clustering parameters
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  # Display parameters
  show_heatmap_legend = TRUE,
  use_raster = ifelse(ncol(kegg_top50_scaled) * nrow(kegg_top50_scaled) > 20000, TRUE, FALSE),
  
  # Overall size settings - modified height parameter, increased heatmap height
  width = unit(max(15, ncol(kegg_top50_scaled) * 0.4 + 4), "cm"),
  height = unit(max(4, nrow(kegg_top50_scaled) * 0.1 + 4), "cm")  # Increased from 0.15 to 0.25, base value from 3 to 4
)

# Draw and save - modified PDF height settings
pdf(file.path(output, "GSVA_KEGG_top50_by_celltype_treatment.pdf"), 
    width = max(10, ncol(kegg_top50) * 0.3 + 6),
    height = max(6, nrow(kegg_top50) * 0.1 + 5))  # Increased from 0.12 to 0.2, base value from 4 to 5

# Draw heatmap, legend at bottom, displayed in two rows
draw(ht, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legends = TRUE,
     ht_gap = unit(2, "cm"),
     # Increase bottom margin to accommodate two rows of legends
     padding = unit(c(20, 20, 60, 20), "mm")  # Top, right, bottom, left margins
)

dev.off()


