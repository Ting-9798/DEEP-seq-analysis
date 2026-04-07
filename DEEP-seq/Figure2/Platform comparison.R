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
library(ggpubr)

# Set working directory
setwd("/data/")

# Define folder path
data_dir <- "/data/"

# Define folder groups
folders_10x <- c("A","10X_v3/","10X_v4/","InDrops_v2/")

# Initialize Seurat object list
seurat.list <- list()

# Loop through folders containing 10X format files
for (folder in folders_10x) {
  # Define current folder path
  folder_path <- file.path(data_dir, folder)
  
  # Manually read files
  barcodes <- read.table(file.path(folder_path, "barcodes.tsv"), header = FALSE)
  features <- read.table(file.path(folder_path, "features.tsv"), header = FALSE)
  matrix <- Matrix::readMM(file.path(folder_path, "matrix.mtx"))
  
  # Create counts matrix required for Seurat object
  # Process feature names to avoid Seurat errors
  gene_names <- make.unique(gsub("_", "-", features[, 2]))
  rownames(matrix) <- gene_names
  colnames(matrix) <- barcodes[, 1]
  print(rownames(matrix))
  
  # Use folder name as sample name prefix
  sample_name <- basename(folder)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = matrix, project = sample_name, min.cells = 10, min.features = 200)
  
  # Add sample prefix to cell names to ensure uniqueness
  seurat_obj <- RenameCells(seurat_obj, new.names = paste(sample_name, colnames(seurat_obj), sep = "_"))
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^hg-MT-|^mm-mt-|^GRCh38-MT-|^GRCm39-mt-")
  
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & nCount_RNA <100000 & percent.mt < 15)
  
  # Add annotations based on sample information
  seurat_obj$orig.ident <- sample_name
  
  # Assign treatment groups
  seurat_obj$treatment <- sample_name
  
  # Add Seurat object to list
  seurat.list[[sample_name]] <- seurat_obj
}


# Merge all Seurat objects
combined_seurat1 <- Reduce(function(x, y) merge(x, y), seurat.list)

# View merged metadata
View(combined_seurat1@meta.data)


# Save merged Seurat object
saveRDS(combined_seurat1, file = "/out/combined_seurat1.rds")


folders_txt <- c("DisCo")
# Initialize Seurat object list
seurat.list <- list()
# Loop through folders containing .txt format files
for (folder in folders_txt) {
  # Define current folder path
  folder_path <- file.path(data_dir, folder)
  # Use file name as sample name prefix
  sample_name <- basename(folder)
  # Get list of .txt files
  txt_files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)
  
  for (txt_file in txt_files) {
    # Read .txt file as expression matrix
    sample_data <- read.table(txt_file, header = TRUE, row.names = 1)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = as.matrix(sample_data), project = sample_name, min.cells = 10, min.features = 200)
    
    # Add sample prefix to cell names to ensure uniqueness
    seurat_obj <- RenameCells(seurat_obj, new.names = paste(sample_name, colnames(seurat_obj), sep = "_"))
    
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^GRCh38_MT-|^GRCm39_mt-")
    
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & nCount_RNA <100000 & percent.mt < 15)
    
    # Add annotations based on sample information
    seurat_obj$orig.ident <- sample_name
    
    # Assign treatment groups
    seurat_obj$treatment <- sample_name
    
    # Add Seurat object to list
    seurat.list[[sample_name]] <- seurat_obj
  }
}

# Merge all Seurat objects
combined_seurat2 <- Reduce(function(x, y) merge(x, y), seurat.list)

# View merged metadata
View(combined_seurat2@meta.data)

# Save merged Seurat object
saveRDS(combined_seurat2, file = "/out/combined_seurat2.rds")



#############################Integration####################################
# Read Seurat objects for 8C and 4C
seurat1 <- readRDS("/out/combined_seurat1.rds")
seurat2<- readRDS("/out/combined_seurat2.rds")

# Create object list
seurat.list <- list(seurat1, seurat2)

combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)

# Normalization and dimensionality reduction
combined_seurat <- ScaleData(combined_seurat, verbose = FALSE)

View(combined_seurat@meta.data)


# Save merged Seurat object
saveRDS(combined_seurat, file = "/out/combined_seurat.rds")



#######################Seurat Analysis#####################
# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)

# Set colors
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# Set working directory
setwd("/out/")
outdir <- "/out/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# Read Seurat object
ScRNA <- readRDS("/out/combined_seurat.rds")

ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("DisCo","InDrops_v2", "SpinDrop","PIP","10X_v3","10X_v4", "A"))
ScRNA$`orig.ident` <- factor(ScRNA$`orig.ident`, levels =  c("DisCo","InDrops_v2", "SpinDrop","PIP","10X_v3","10X_v4", "A"))

# Generate violin plot showing QC metrics
#pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 9, height = 6)
#VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol =2 , pt.size = 0)
#dev.off()



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
  coord_cartesian(ylim = c(0, 25000)) +  ## After downsampling, upper limit ≈20,000
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


