# Clear environment
rm(list = ls())

# Load required packages
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
#install.packages("tidydr")


setwd("/data/")
# Define folder path
data_dir <- "/data/"

folders <- c("0808-WF-25/","0808-WF-40/","0808-WF-50/")

# Initialize Seurat object list
seurat.list <- list()

# Loop through each subfolder to read expression matrix files and perform Seurat analysis
for (folder in folders) {
  folder_path <- file.path(data_dir, folder)
  
  # Read 10X data format
  sample_data <- Read10X(data.dir = folder_path)
  sample_name <- basename(folder)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 3, min.features = 200)
  
  # Add mitochondrial gene percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^GRCh38_MT-|^GRCm39_mt-")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hbb-|^Hba-")
  
  # Quality control filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 &percent.mt < 25)   #nFeature_RNA > 200 & nFeature_RNA < 8000 & 
  
  
  # Add processed Seurat object to list
  seurat.list[[sample_name]] <- seurat_obj
}

# Integrate data
#anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)
#combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:20)

combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)
combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
View(combined_seurat@meta.data)

# Save merged Seurat object
saveRDS(combined_seurat, file = "/out/combined_seurat.rds")


#######################Seurat Analysis#####################
# Set output directory
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("/out/")
outdir <- "/out/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# Construct full file path
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("WF-25", "WF-40","WF-50"))
ScRNA$`orig.ident` <- factor(ScRNA$`orig.ident`, levels = c("WF-25", "WF-40","WF-50"))

## Calculate red blood cell proportion
ScRNA[["percent.hb"]] <- PercentageFeatureSet(ScRNA, pattern = "^Hbb-|^Hba-")

# Generate violin plot showing QC metrics
#pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 5)
#VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4, group.by = "treatment", pt.size = 0,cols = col)
#dev.off()

## Remove red blood cell contamination
ScRNA <- subset(ScRNA, subset = percent.hb < 1)

# Generate violin plot showing QC metrics
svg(paste(OUTPUT, "QC-BoxPlot.svg"), width = 8, height = 6)
p1 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nFeature_RNA, color = treatment)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nFeature_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p2 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nCount_RNA, color = treatment)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nCount_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

CombinePlots(plots = list(p1, p2))
dev.off()


pdf(paste(OUTPUT, "QC-ViolinPlot-sample.pdf"), width = 16, height = 6)
# Violin plot 1: nFeature_RNA
p1 <- ggplot(data = ScRNA@meta.data, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values =c("#99CCCC","#A4CDE1",'#3399CC','#57C3F3')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Violin plot 2: nCount_RNA
p2 <- ggplot(data = ScRNA@meta.data, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values =c("#99CCCC","#A4CDE1",'#3399CC','#57C3F3')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p3 <- ggplot(data = ScRNA@meta.data, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "percent.mt", x = "", y = "") +
  scale_fill_manual(values =c("#99CCCC","#A4CDE1",'#3399CC','#57C3F3')) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Combine plots
library(patchwork)
(p1 | p2 |p3) + plot_layout(ncol = 3)
dev.off()






########## Add mean labels ##########
library(ggplot2)
library(dplyr)
library(patchwork)
#install.packages("shadowtext")
#library(shadowtext)

# —— Color mapping (ensuring same treatment group has consistent color) ——
treatment_levels <- unique(ScRNA@meta.data$treatment)
treatment_colors <- c()



## ---- Calculate mean values for each treatment (gene count and UMI) ----
means_df <- ScRNA@meta.data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(
    mean_features = mean(nFeature_RNA, na.rm = TRUE),
    mean_counts   = mean(nCount_RNA,   na.rm = TRUE),
    .groups = "drop"
  )


# Unify treatment factor levels in main data and mean table to the same levels
#ScRNA@meta.data <- ScRNA@meta.data %>% dplyr::mutate(treatment = factor(treatment, levels = c("10X_v3", "10X_v4","SNJ","WF")))

#means_df <- means_df %>% dplyr::mutate(treatment = factor(treatment, levels = c("10X_v3", "10X_v4","SNJ","WF")))



# Helper function: round numbers to no decimals or 1 decimal (adjustable)
fmt_num <- function(x) ifelse(x >= 100, sprintf("%.0f", x), sprintf("%.1f", x))

# —— Output PDF —— 
pdf(paste(OUTPUT, "QC-ViolinPlot.pdf"), width = 9, height = 5)

# Violin plot 1: nFeature_RNA
p1 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nFeature_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  # Add mean text (consistent with violin color)
  geom_label(
    data = means_df,
    aes(x = treatment, y = mean_features, label = fmt_num(mean_features)),
    color = "black", fill = "grey90",
    vjust = -0.6, size = 5, inherit.aes = FALSE,
    fontface = "bold", 
    label.size = 0,
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines")
  ) +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col, guide = "none") +
  #coord_cartesian(ylim = c(0, 11000)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title  = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Violin plot 2: nCount_RNA
p2 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nCount_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  # Add mean labels (grey background)
  geom_label(
    data = means_df,
    aes(x = treatment, y = mean_counts, label = fmt_num(mean_counts)),
    color = "black", fill = "grey90",
    vjust = -0.6, size = 5, inherit.aes = FALSE,
    fontface = "bold", 
    label.size = 0,
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines")
  ) +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col, guide = "none") +
  #coord_cartesian(ylim = c(0, 80000)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title  = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Combine plots
(p1 | p2) + plot_layout(ncol = 2)
dev.off()






## ========= Downsample based on equivalent 20,000 reads and generate QC violin plots =========
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## Target reads (per cell)
target_reads <- 20000L
set.seed(1234)

## Get raw count matrix (sparse)
counts <- ScRNA@assays$RNA@counts
stopifnot(inherits(counts, "dgCMatrix"))

## Calculate total UMI per cell
libsizes <- Matrix::colSums(counts)

## --- Prioritize using DropletUtils::downsampleMatrix (if available) ---
downsample_with_dropletutils <- function(x, target, libs) {
  if (requireNamespace("DropletUtils", quietly = TRUE)) {
    ## Proportional downsampling only for columns >= target; columns < target remain unchanged
    prop <- pmin(1, target / libs)
    ## DropletUtils::downsampleMatrix performs random dilution by multiplying by prop per column
    y <- DropletUtils::downsampleMatrix(x, prop = prop)
    return(y)
  } else {
    return(NULL)
  }
}

## --- Pure R sparse implementation (fallback): perform multinomial sampling to min(total, target) per column ---
downsample_sparse_multinom <- function(x, target) {
  x <- as(x, "dgCMatrix")
  p <- x@p; i <- x@i; xv <- x@x
  ncol_x <- ncol(x)
  out_list_i <- vector("list", ncol_x)
  out_list_v <- vector("list", ncol_x)
  
  for (cc in seq_len(ncol_x)) {
    idx <- seq.int(p[cc] + 1L, p[cc + 1L])
    if (length(idx) == 0L) {
      out_list_i[[cc]] <- integer(0)
      out_list_v[[cc]] <- numeric(0)
      next
    }
    genes <- i[idx]              # row indices (starting from 0)
    vals  <- xv[idx]             # non-zero counts for this column
    tot   <- sum(vals)
    draw  <- min(tot, target)
    
    if (tot == 0 || draw == tot) {
      ## tot==0 or originally <= target (no downsampling, keep unchanged)
      keep <- vals > 0
      out_list_i[[cc]] <- genes[keep]
      out_list_v[[cc]] <- vals[keep]
    } else {
      ## Multinomial sampling to "draw" reads
      prob <- vals / tot
      sampled <- as.vector(rmultinom(1L, size = draw, prob = prob))
      keep <- sampled > 0L
      out_list_i[[cc]] <- genes[keep]
      out_list_v[[cc]] <- sampled[keep]
    }
  }
  
  ## Assemble back to dgCMatrix
  new_p <- integer(ncol_x + 1L)
  new_p[1] <- 0L
  for (cc in seq_len(ncol_x)) new_p[cc + 1L] <- new_p[cc] + length(out_list_i[[cc]])
  
  new_i <- do.call(c, out_list_i)
  new_x <- as.numeric(do.call(c, out_list_v))
  
  new_mat <- new("dgCMatrix",
                 i = as.integer(new_i),
                 p = as.integer(new_p),
                 x = new_x,
                 Dim = dim(x),
                 Dimnames = dimnames(x))
  return(new_mat)
}

## Perform downsampling
down_counts <- downsample_with_dropletutils(counts, target_reads, libsizes)
if (is.null(down_counts)) {
  message("DropletUtils not detected, using built-in multinomial implementation (may be slower)...")
  down_counts <- downsample_sparse_multinom(counts, target_reads)
}

## Recalculate QC metrics based on downsampled matrix
down_nCount   <- Matrix::colSums(down_counts)
down_nFeature <- Matrix::colSums(down_counts > 0)

## Assemble metadata for plotting (retain original treatment, replace QC with downsampled values)
meta_down <- ScRNA@meta.data %>%
  dplyr::mutate(
    nCount_RNA   = as.numeric(down_nCount[colnames(ScRNA)]),
    nFeature_RNA = as.numeric(down_nFeature[colnames(ScRNA)])
  )

## Calculate group means (after downsampling)
means_df_down <- meta_down %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(
    mean_features = mean(nFeature_RNA, na.rm = TRUE),
    mean_counts   = mean(nCount_RNA,   na.rm = TRUE),
    .groups = "drop"
  )

fmt_num <- function(x) ifelse(x >= 100, sprintf("%.0f", x), sprintf("%.1f", x))

## Colors follow the col vector defined above
## Output directory and file names

# downsampled to 20000 reads/cell

pdf(file = paste0(OUTPUT, "QC-ViolinPlot-downsampled(20000 readings per cell).pdf"), width = 9, height = 5)

p1 <- ggplot(data = meta_down, aes(x = treatment, y = nFeature_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  geom_label(
    data = means_df_down,
    aes(x = treatment, y = mean_features, label = fmt_num(mean_features)),
    color = "black", fill = "grey90",
    vjust = -0.6, size = 5, inherit.aes = FALSE,
    fontface = "bold", 
    label.size = 0,
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines")
  ) +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col, guide = "none") +
  coord_cartesian(ylim = c(0, 9000)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title  = element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p2 <- ggplot(data = meta_down, aes(x = treatment, y = nCount_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  geom_label(
    data = means_df_down,
    aes(x = treatment, y = mean_counts, label = fmt_num(mean_counts)),
    color = "black", fill = "grey90",
    vjust = -0.6, size = 5, inherit.aes = FALSE,
    fontface = "bold", 
    label.size = 0,
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines")
  ) +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col, guide = "none") +
  coord_cartesian(ylim = c(-20, 25000)) +  ## After downsampling, upper limit ≈20,000
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title  = element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

(p1 | p2)
dev.off()





#### 5. Expression normalization ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Calculate highly variable genes
ScRNA <- FindVariableFeatures(ScRNA, selection.method = "vst",
                              nfeatures = 2000)

# Display highly variable genes
pdf(paste(OUTPUT,"variable gene.pdf"),width = 9,height = 6)
top10 <- head(VariableFeatures(ScRNA), 10) 
plot1 <- VariableFeaturePlot(ScRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=3)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()

#### 6. Normalization and PCA dimensionality reduction ####
# Normalization
ScRNA<-ScaleData(ScRNA)

# Run PCA
ScRNA<-RunPCA(ScRNA,npcs = 30)

pdf(paste(OUTPUT,"Dimplot.pdf"),width = 9,height = 6)
DimPlot(object = ScRNA, reduction = "pca", pt.size = .1, group.by = "treatment",cols = col)
dev.off()

pdf(paste(OUTPUT,"vlnplot.pdf"),width = 9,height = 6)
VlnPlot(object = ScRNA, features = "PC_1", group.by = "treatment", pt.size = 0,cols = col)
dev.off()

# PCA visualization
pdf(paste(OUTPUT, "DimHeatmap.pdf"),width = 9,height = 6)
DimHeatmap(ScRNA, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# Evaluate PC dimensions
pdf(paste0(OUTPUT,"PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(ScRNA)
dev.off()


# Batch correction
#install.packages("harmony")
library(harmony)
ScRNA<-RunHarmony(ScRNA,group.by.vars = c("treatment"), 
                  plot_convergence = TRUE, verbose = TRUE)

pdf(paste(OUTPUT,  "Dimplot-corret.pdf"),width = 12,height = 6)
DimPlot(object = ScRNA, reduction = "harmony",
        pt.size = 0.1, group.by = "treatment")
dev.off()

pdf(paste(OUTPUT, "vlnplot-corret.pdf"),width = 12,height = 6)
VlnPlot(object = ScRNA, features = "harmony_1", 
        group.by = "treatment", pt.size =0)
dev.off()


save(ScRNA, file = "ScRNA（after batch correction before clustering）.RData")



#### 7. Cell clustering and annotation ####

col <- c('#6A4C93',"#A4CDE1",'#FF9999',"#CC99CC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#6A4C93',
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA（after batch correction before clustering）.RData")
# Cell clustering
#ScRNA <- ScRNA %>% 
#  RunUMAP(dims = 1:20) %>% 
#  RunTSNE(dims = 1:20) %>%
#  FindNeighbors(dims = 1:20)

ScRNA <- ScRNA %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,spread = 2) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)

ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))
#pdf(paste(OUTPUT, "clustree.pdf"),width=10,height=9)
#library(clustree)
#clustree(ScRNA)
#dev.off()

#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident##Select your resolution based on the clustering tree
table(Idents(ScRNA))

#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))

# Ensure "treatment" factor levels are ordered as Non-infected and Infected
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Display clustering, ordered by Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 6*length(unique(ScRNA$treatment)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=0.1, label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Display clustering, ordered by Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 6*length(unique(ScRNA$orig.ident)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=0.1, label = TRUE, repel = TRUE, split.by = "orig.ident", cols = col)
dev.off()

# Generate separate umap plot
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap",pt.size=0.1, label = TRUE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

DimPlot(ScRNA, reduction = "umap", pt.size=0.1, label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())
dev.off()

# Display tumor and normal together
pdf(paste(OUTPUT, "cluster-diff_umap.pdf"),width=6,height=6)
DimPlot(ScRNA, repel = TRUE,pt.size=0.1, 
        reduction = "umap",
        group.by ="treatment")+
  scale_color_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")
dev.off()


saveRDS(ScRNA, "ScRNA（after clustering）.rds")



# Select top 6 cell clusters to display marker gene expression trends
output <- paste(outdir,'cell_markers', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
ScRNA <- readRDS(file_path)

# Find marker genes
ScRNA.markers <- FindAllMarkers(ScRNA, only.pos = TRUE,   ### only.pos = TRUE: only find upregulated genes
                                min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ScRNA.markers,paste0(output,"./ScRNA.all.markers.csv"))

dim.use <-1:30
top5 <- ScRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top5,file=paste0(output,"/top20_marker_genes_tsne_",max(dim.use),"PC.csv"))

pdf(paste0(output,"/Heatmap_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 25,height = 20)
DoHeatmap(ScRNA, features = top5$gene,size = 2)+
  scale_fill_gradientn(colors = c("#437eb8", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18),  # Increase Y-axis title text size
        legend.text = element_text(size = 14),  # Modify legend text size
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 30))  # Modify grouping label color and size
dev.off()

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 80,height = 7)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18))  # Increase Y-axis title text size
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=80*dpi,h=7*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 20),  # Increase X-axis title text size
        axis.title.y = element_text(size = 20))  # Increase Y-axis title text size
dev.off()

