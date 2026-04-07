## Visualization
rm(list=ls())

#BiocManager::install("RcisTarget")
#devtools::install_github("aertslab/SCopeLoomR")
devtools::install_github("aertslab/SCENIC")

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)

setwd("/out/cancer cell/pyscenic/")


####### Prepare files for pySCENIC #########
data = readRDS("./celltype.rds")

# Note: the matrix must be transposed, otherwise an error will occur
write.csv(t(as.matrix(data@assays$RNA@counts)), file = "for.scenic.data.csv")



#### 1. Extract information from out_SCENIC.loom
loom <- open_loom('out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)



# Read sequencing data
seurat.data = readRDS("./celltype.rds")
seurat.data

seurat.data@meta.data$Barcode = colnames(seurat.data)
seurat.data = subset(seurat.data, Barcode %in% colnames(regulonAUC)) # Extract sampled cells
seurat.data
table(seurat.data@meta.data$Seu_Clusters)

DimPlot(seurat.data, reduction = "umap", label = T) 

## AUC visualization
sub_regulonAUC <- regulonAUC[, match(colnames(seurat.data), colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data

# Check whether they are identical
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$Meta_3class))

## Append group name after each cell
seurat.data@meta.data$Celltype_Group <- paste0(seurat.data@meta.data$Meta_3class, "_", seurat.data@meta.data$treatment)
table(seurat.data@meta.data$Celltype_Group)

cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$Meta_3class)

Celltype_Group <- data.frame(row.names = colnames(seurat.data), 
                             celltype = seurat.data$Celltype_Group)

head(cellTypes)
head(Celltype_Group)
sub_regulonAUC[1:4,1:4] 

# Save intermediate objects
save(sub_regulonAUC, cellTypes, Celltype_Group, cellClusters, seurat.data,
     file = 'for_rss_and_visual.Rdata')



### 4.1. Mean TF activity
# Check the average transcription factor activity across different single-cell subclusters
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[, selectedResolution])

# Remove extended regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[, cells]))

# Scale expression.
# The scale() function normalizes by column, so regulonActivity_byGroup should be transposed
# to make cells rows and genes columns
# Reference: https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale = T)) 
# Scale the same regulon across different clusters

regulonActivity_byGroup_Scaled = na.omit(regulonActivity_byGroup_Scaled)


Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2, to=2, length=11), rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows       = "ward.D2",
  clustering_method_columns    = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)


### 4.2. Use RSS to inspect specific TFs
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), 
               cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss = na.omit(rss) 
rssPlot <- plotRSS(rss)

rss_treatment <- calcRSS(AUC = getAUC(sub_regulonAUC), 
                         cellAnnotation = Celltype_Group[colnames(sub_regulonAUC), "celltype"])

#rss_treatment=na.omit(rss_treatment) 
#rssPlot <- plotRSS(rss_treatment)
#plotly::ggplotly(rssPlot$plot)




### Dot plot
regulon_AUC <- regulonAUC@NAMES
seurat.data@meta.data = cbind(seurat.data@meta.data, t(assay(sub_regulonAUC[regulon_AUC, ])))
# View(seurat.data@meta.data)

top10_TFs <- top10_TFs_df$TF

# Method 1: directly calculate average activity from metadata (recommended)
# Because TF activity data have already been stored in metadata
calculate_tf_activity_from_metadata <- function(seurat_obj, tf_list, group_var = "Meta_3class") {
  # Extract TF activity data from metadata
  metadata <- seurat_obj@meta.data
  tf_columns <- tf_list[tf_list %in% colnames(metadata)]
  
  if(length(tf_columns) == 0) {
    stop("Specified TFs were not found in metadata")
  }
  
  # Calculate average activity by group
  groups <- unique(metadata[[group_var]])
  result <- matrix(NA, nrow = length(tf_columns), ncol = length(groups))
  rownames(result) <- tf_columns
  colnames(result) <- groups
  
  for(tf in tf_columns) {
    for(group in groups) {
      group_cells <- which(metadata[[group_var]] == group)
      result[tf, group] <- mean(metadata[group_cells, tf], na.rm = TRUE)
    }
  }
  
  return(result)
}

# Calculate TF activity
tf_activity <- calculate_tf_activity_from_metadata(seurat.data, top10_TFs, "Meta_3class")

# Check the result
print("TF activity matrix:")
print(tf_activity)

# Calculate Tumor/Normal differences and sort
if("Tumor" %in% colnames(tf_activity) & "Normal" %in% colnames(tf_activity)) {
  tf_differences <- tf_activity[, "Tumor"] - tf_activity[, "Normal"]
  sorted_tfs <- names(sort(tf_differences, decreasing = TRUE))
  
  # Set treatment order (Tumor first)
  treatment_order <- c("Tumor", "Normal")
  seurat.data$treatment <- factor(seurat.data$treatment, levels = treatment_order)
  
  cat("TF ranking information (sorted by Tumor-Normal difference from high to low):\n")
  for(i in 1:length(sorted_tfs)) {
    tf <- sorted_tfs[i]
    tumor_act <- tf_activity[tf, "Tumor"]
    normal_act <- tf_activity[tf, "Normal"]
    diff <- tumor_act - normal_act
    cat(sprintf("%2d. %s: Tumor=%.3f, Normal=%.3f, Difference=%.3f\n", 
                i, tf, tumor_act, normal_act, diff))
  }
} else {
  # If the column names are not Tumor/Normal, use automatic detection
  cat("Detected treatment groups:\n")
  print(colnames(tf_activity))
  sorted_tfs <- top10_TFs  # Use original order
}

# Draw dot plot
p_sorted <- DotPlot(seurat.data, 
                    features = sorted_tfs, 
                    group.by = 'Meta_3class') +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, size = 12),
    axis.text.y = element_text(face = "italic", size = 12),  # TF names in italics
    #legend.direction = "horizontal", 
    #legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, 
    y = NULL,
    title = NULL,
    subtitle = NULL,
  ) +
  coord_flip() +  # Swap X and Y axes
  scale_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#C71000FF",
    midpoint = 0
  )

# Display plot
print(p_sorted)

# Save sorted plot
ggsave("TF_activity_by_treatment_sorted_dotplot.pdf", 
       plot = p_sorted, 
       width = 5, 
       height = 5)





### Directly select the top five positively regulated transcription factors in each treatment group and draw a dot plot

# Calculate RSS values of TFs within each treatment group
treatment_rss <- calcRSS(AUC = getAUC(sub_regulonAUC), 
                         cellAnnotation = Celltype_Group[colnames(sub_regulonAUC), "celltype"])

treatment_rss <- na.omit(treatment_rss)

# Extract positively regulated transcription factors
positive_tfs <- rownames(treatment_rss)[grepl("\\(\\+\\)", rownames(treatment_rss))]
treatment_rss_pos <- treatment_rss[positive_tfs, ]

# Get all treatment groups
treatment_groups <- unique(seurat.data@meta.data$Meta_3class)

# Select top 5 positive TFs for each treatment group
top_tfs_per_treatment <- list()

for(treatment_group in treatment_groups) {
  # Extract RSS values for the current treatment group
  treatment_data <- treatment_rss_pos[, grep(treatment_group, colnames(treatment_rss_pos)), drop = FALSE]
  
  if(ncol(treatment_data) > 0) {
    # Calculate average RSS value (if multiple celltypes belong to the same treatment)
    if(ncol(treatment_data) > 1) {
      treatment_means <- rowMeans(treatment_data, na.rm = TRUE)
    } else {
      treatment_means <- treatment_data[, 1]
    }
    
    # Sort by RSS value and take top 5
    sorted_tfs <- names(sort(treatment_means, decreasing = TRUE))
    top_tfs <- head(sorted_tfs, 5)
    
    top_tfs_per_treatment[[treatment_group]] <- top_tfs
    
    cat("Treatment:", treatment_group, "\n")
    cat("Top 5 positive TFs:", paste(top_tfs, collapse = ", "), "\n")
    cat("RSS values:", paste(round(treatment_means[top_tfs], 3), collapse = ", "), "\n\n")
  }
}

# Combine top TFs from all treatment groups (remove duplicates)
all_top_tfs <- unique(unlist(top_tfs_per_treatment))

cat("Total unique top TFs across all treatments:", length(all_top_tfs), "\n")
print(all_top_tfs)

# If the number of selected TFs is insufficient, supplement with some highly active TFs
if(length(all_top_tfs) < 5) {
  # Calculate average activity of all positively regulated TFs
  all_positive_tfs <- rownames(treatment_rss_pos)
  overall_means <- rowMeans(treatment_rss_pos, na.rm = TRUE)
  additional_tfs <- names(sort(overall_means, decreasing = TRUE))[1:min(10, length(overall_means))]
  additional_tfs <- setdiff(additional_tfs, all_top_tfs)
  all_top_tfs <- c(all_top_tfs, head(additional_tfs, 5 - length(all_top_tfs)))
}

# Add TF activity data to the metadata of the Seurat object
regulon_AUC <- regulonAUC@NAMES
seurat.data@meta.data = cbind(seurat.data@meta.data, 
                              t(assay(sub_regulonAUC[regulon_AUC, ])))

# Set treatment order (alphabetical or custom order)
treatment_order <- sort(unique(seurat.data@meta.data$treatment))
seurat.data$treatment <- factor(seurat.data$treatment, levels = treatment_order)

# Draw dot plot
p_treatment <- DotPlot(seurat.data, 
                       features = all_top_tfs, 
                       group.by = 'Meta_3class') +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(hjust = 1, angle = 45, vjust = 1, size = 12),
    axis.text.y = element_text(size = 14),  # TF names in italics: face = "italic",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    #legend.direction = "horizontal", 
    #legend.position = "bottom",
  ) +
  labs(
    x = NULL, 
    y = NULL,
    title = NULL,
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  scale_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#C71000FF",
    midpoint = 0,
    name = "Average\nExpression"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = "Percent\nExpressed"
  )

# Display plot
print(p_treatment)

# Save plot
ggsave("TF_activity_by_treatment_top5_dotplot.pdf", 
       plot = p_treatment, 
       width = 8, 
       height = 2.8)



## -------------------------------
## Draw UMAPs (AUC) of top5 TFs for each Meta_3class
## -------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)  # Used for plot_annotation; if unavailable, install with: install.packages("patchwork")

# Ensure UMAP coordinates are available (usually reduction = "umap")
# If your dimensional reduction name is not "umap", for example "umap_scRNA",
# modify the reduction name below accordingly

umap_plots <- list()

for (treatment_group in names(top_tfs_per_treatment)) {
  
  tf_use <- top_tfs_per_treatment[[treatment_group]]
  tf_use <- tf_use[tf_use %in% colnames(seurat.data@meta.data)]
  
  if (length(tf_use) == 0) {
    message("No matching TF columns in metadata for group: ", treatment_group)
    next
  }
  
  p_group <- FeaturePlot(
    seurat.data,
    features  = tf_use,
    reduction = "umap",
    cols      = c("lightgrey", "red"),
    order     = TRUE,
    combine   = TRUE,
    ncol      = length(tf_use)
  ) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 18)
    )
  
  umap_plots[[treatment_group]] <- p_group
  print(p_group)
  
  ggsave(
    filename = paste0("UMAP_AUC_", treatment_group, "_top5_TFs.pdf"),
    plot     = p_group,
    width    = 4 * min(length(tf_use), 5),
    height   = 4
  )
}




