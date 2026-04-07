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
#install.packages("tidydr")

# Set working directory
setwd("/data/")

# Define folder path
data_dir <- "/data/"

# Initialize Seurat object list
seurat.list <- list()

# Define subfolder list
folders <- c("A")

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
  
  # Data normalization and identification of highly variable genes
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # Add treatment column
  seurat_obj$treatment <- sample_name
  
  # Add processed Seurat object to the list
  seurat.list[[sample_name]] <- seurat_obj
}

# Save merged Seurat object
saveRDS(seurat_obj, file = "/out/seurat_obj.rds")

####################### Seurat Analysis #####################
col <- c("#FFFFCC",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# Set output directory
setwd("/out/")

outdir <- "/out/"
MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

ScRNA <- readRDS("/out/seurat_obj.rds")

# Calculate the UMI ratio for human and mouse
# Get row indices matching cell names
human_genes <- grepl("GRCh38", rownames(ScRNA@assays$RNA@counts))
mouse_genes <- grepl("GRCm39", rownames(ScRNA@assays$RNA@counts))

# Calculate UMI counts for each cell
human_umi <- Matrix::colSums(ScRNA@assays$RNA@counts[human_genes, ])
mouse_umi <- Matrix::colSums(ScRNA@assays$RNA@counts[mouse_genes, ])
total_umi <- human_umi + mouse_umi
human_umi_ratio <- human_umi / total_umi

# Convert the calculated UMI information to a data frame and ensure row names are cell names
umi_data <- data.frame(human_umi = human_umi, mouse_umi = mouse_umi, total_umi = total_umi, human_umi_ratio = human_umi_ratio)
umi_data <- umi_data[rownames(ScRNA@meta.data), ]  # Ensure row names match the cell names in the Seurat object

# Add UMI information to the meta.data of the Seurat object
ScRNA <- AddMetaData(ScRNA, metadata = umi_data)

# Determine species type based on ratio
ScRNA$species <- ifelse(ScRNA$human_umi_ratio > 0.8, "Human", 
                        ifelse(ScRNA$human_umi_ratio < 0.2, "Mouse", "Multiplet"))

# Calculate mixture ratio
mixed_ratio <- mean(ScRNA$species == "Multiplet")

# Draw scatter plot
pdf(paste(OUTPUT, "Mixture Ratio.pdf"), width = 7, height = 7)
ggplot(data = as.data.frame(ScRNA@meta.data), aes(x = human_umi, y = mouse_umi, color = species)) +
  geom_point(alpha = 0.5, size = 0.2) +
  scale_color_manual(values = c("Human" = "#966392", "Mouse" ="#3C8487", "Multiplet" = "gray")) +
  labs(x = "Human UMI", y = "Mouse UMI", color = "Species") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("Mixture Ratio: ", round(mixed_ratio * 100, 2), "%"), 
           hjust = 1.1, vjust = 2, size = 7, fontface = "bold") +
  theme_minimal() +
  theme(axis.title = element_text(size = 25),  
        axis.text = element_text(size = 25),   
        legend.position = "top",
        legend.title = element_text(size = 24), 
        legend.text = element_text(size = 22), 
        panel.background = element_blank(),        
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  guides(color = guide_legend(override.aes = list(size = 2))) 

dev.off()

svg(paste(OUTPUT, "Mixture Ratio.svg"), width = 7, height = 7)
ggplot(data = as.data.frame(ScRNA@meta.data), aes(x = human_umi, y = mouse_umi, color = species)) +
  geom_point(alpha = 0.5, size = 0.2) +
  scale_color_manual(values = c("Human" = "#966392", "Mouse" ="#3C8487", "Multiplet" = "gray")) +
  labs(x = "Human UMI", y = "Mouse UMI", color = "Species") +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("Mixture Ratio: ", round(mixed_ratio * 100, 2), "%"), 
           hjust = 1.1, vjust = 2, size = 7, fontface = "bold") +
  theme_minimal() +
  theme(axis.title = element_text(size = 25),  
        axis.text = element_text(size = 25),   
        legend.position = "top",
        legend.title = element_text(size = 24), 
        legend.text = element_text(size = 22), 
        panel.background = element_blank(),        
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  guides(color = guide_legend(override.aes = list(size = 2))) 

dev.off()


