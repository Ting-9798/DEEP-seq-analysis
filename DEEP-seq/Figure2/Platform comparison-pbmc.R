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

folders <- c("A/","PBMC-10X/","PBMC-MZ/")
# Initialize Seurat object list
seurat.list <- list()

# Loop through each subfolder to read expression matrix files and perform Seurat analysis
for (folder in folders) {
  folder_path <- file.path(data_dir, folder)
  
  # Read 10X data format
  sample_data <- Read10X(data.dir = folder_path, gene.column = 2)
  sample_name <- basename(folder)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 3, min.features = 200)
  
  # Add mitochondrial gene percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Quality control filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
  
  # Add processed Seurat object to list
  seurat.list[[sample_name]] <- seurat_obj
}

combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)
combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
View(combined_seurat@meta.data)

# Save merged Seurat object
saveRDS(combined_seurat, file = "/out/combined_seurat.rds")


#######################Seurat Analysis#####################
# Set output directory
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
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


# Set output directory
plot_dir <- paste(outdir,"scatter_plots", sep='/')
dir.create(plot_dir, showWarnings = FALSE)

# Extract orig.ident for all samples
metadata <- ScRNA@meta.data
samples <- unique(metadata$treatment)


# Get counts (sparse matrix)
expr_data <- ScRNA@assays$RNA@counts
# Calculate total counts per column (cell)
total_counts <- Matrix::colSums(expr_data)
# Construct a sparse matrix form of TP10K
TP10K <- t(t(expr_data) / total_counts) * 10000
# Convert to log2(TP10K + 1), maintaining sparsity
log2_TP10K <- log1p(TP10K) / log(2)  # Using log1p avoids precision issues with decimals


library(cowplot)
library(ggplot2)

# Group and generate pairwise correlation scatter plots
generate_combined_plots <- function(samples, plot_dir) {
  plot_list <- list()
  plot_count <- 1
  
  # Iterate through pairwise combinations of samples
  for (i in 1:(length(samples) - 1)) {
    for (j in (i + 1):length(samples)) {
      sample1 <- samples[i]
      sample2 <- samples[j]
      
      # Get gene expression values for both samples
      sample1_data <- log2_TP10K[, metadata$treatment == sample1]
      sample2_data <- log2_TP10K[, metadata$treatment == sample2]
      
      # Calculate mean expression values
      sample1_mean <- rowMeans(sample1_data)
      sample2_mean <- rowMeans(sample2_data)
      
      # Create scatter plot data frame
      plot_data <- data.frame(Sample1 = sample1_mean, Sample2 = sample2_mean)
      
      # Calculate Pearson correlation coefficient and its p-value
      correlation_test <- cor.test(sample1_mean, sample2_mean, method = "pearson")
      correlation <- correlation_test$estimate
      p_value <- correlation_test$p.value
      
      # Format p-value display
      format_p_value <- function(p) {
        if (p < 0.001) {
          return("P < 0.001")
        } else if (p < 0.01) {
          return(paste0("P = ", format(p, digits = 1, scientific = FALSE)))
        } else {
          return(paste0("P = ", round(p, 3)))
        }
      }
      
      p_value_text <- format_p_value(p_value)
      
      
      # Fit linear regression model
      fit <- lm(Sample2 ~ Sample1, data = plot_data)
      slope <- coef(fit)[2]
      intercept <- coef(fit)[1]
      r2 <- summary(fit)$r.squared
      
      # Get F-test p-value from linear regression
      fit_summary <- summary(fit)
      f_statistic <- fit_summary$fstatistic
      lm_p_value <- pf(f_statistic[1], f_statistic[2], f_statistic[3], lower.tail = FALSE)
      lm_p_value_text <- format_p_value(lm_p_value)
      
      # Create text label containing statistical information
      stats_label <- paste0(
        "r = ", round(correlation, 3), 
        " (", p_value_text, ")"
      )
      
      # Scatter plot
      plot <- ggplot(data = plot_data, aes(x = Sample1, y = Sample2)) +
        geom_hex(bins = 100) +
        scale_fill_gradientn(
          colors = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF"),
          name = "Count"
        ) +
        # Regression line
        geom_smooth(method = "lm", se = TRUE, color = "#440154FF", size = 1.2, alpha = 0.2) +
        theme_minimal(base_size = 14) +
        labs(
          title = paste(sample1, "vs", sample2),
          subtitle = paste(stats_label),
          x = paste(sample1, "log2(TP10K)"),
          y = paste(sample2, "log2(TP10K)")
        ) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
          plot.subtitle = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 2)
        )
      
      # Add current plot to list
      plot_list[[plot_count]] <- plot
      plot_count <- plot_count + 1
    }
  }
  
  # Calculate appropriate figure dimensions
  n_plots <- length(plot_list)
  n_cols <- 3
  n_rows <- ceiling(n_plots / n_cols)
  
  # Combine all plots
  combined_plot <- plot_grid(plotlist = plot_list, ncol = n_cols, align = "v")
  
  # Save output
  output_file_base <- paste0(plot_dir, "/All_Samples_combined_scatter_plots")
  
  ggsave(
    filename = paste0(output_file_base, ".pdf"),
    plot = combined_plot,
    width = 16,  # Increased width to accommodate text labels
    height = 6 * n_rows,  # Adjust height based on number of rows
    limitsize = FALSE
  )
  
  ggsave(
    filename = paste0(output_file_base, ".svg"),
    plot = combined_plot,
    width = 16,
    height = 6 * n_rows,
    limitsize = FALSE
  )
  
  
  return(combined_plot)
}

# Generate pairwise scatter plots for all samples
generate_combined_plots(samples, plot_dir)



########## Add mean labels ##########
library(ggplot2)
library(dplyr)
library(patchwork)
#install.packages("shadowtext")
#library(shadowtext)

# —— Color mapping (ensuring same treatment group has consistent color) ——
treatment_levels <- unique(ScRNA@meta.data$treatment)
treatment_colors <- c()


treatment_colors <- sapply(treatment_levels, function(lvl) {
  if (grepl("10X$", lvl)) {
    "#1E90FF"
  } else if (grepl("MZ$", lvl)) {
    "#33CCCC"
  } else if (grepl("WF$", lvl)) {
    "#FF3366"
  } else {
    "#9E9E9E"
  }
})

names(treatment_colors) <- treatment_levels

## ---- Calculate mean values for each treatment (gene count and UMI) ----
means_df <- ScRNA@meta.data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(
    mean_features = mean(nFeature_RNA, na.rm = TRUE),
    mean_counts   = mean(nCount_RNA,   na.rm = TRUE),
    .groups = "drop"
  )


# Unify treatment factor levels in main data and mean table to the same levels
ScRNA@meta.data <- ScRNA@meta.data %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("10X", "MZ","WF")))

means_df <- means_df %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("10X", "MZ","WF")))



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
    aes(x = treatment, y = mean_features,
        label = paste0(fmt_num(mean_features)),
        fill = NULL, color = treatment),
    vjust = -0.6, size = 5, inherit.aes = FALSE,
    fontface = "bold", 
    label.size = 0,             # Remove border
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines"),
    fill = "grey90") +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 7000)) +
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
    aes(x = treatment, y = mean_counts,
        label = paste0(fmt_num(mean_counts)),
        fill = NULL, color = treatment),
    vjust = -0.6, size = 5, inherit.aes = FALSE,
    fontface = "bold", 
    label.size = 0,             # Remove border
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines"),
    fill = "grey90") +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 40000)) +
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



#QC: Correlation of gene count with mitochondrial genes and RNA count
pdf(paste(OUTPUT,"cor-plot.pdf"),width = 15,height = 6)
plot1 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = col)
plot2 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols = col)
CombinePlots(plots = list(plot1, plot2),legend = "right")
dev.off()

# Calculate cell cycle scores
pdf(paste(OUTPUT, "cellcycle.pdf"), width = 9, height = 6)
s.genes <- cc.genes$s.genes    ##S phase
g2m.genes <- cc.genes$g2m.genes    ##G2/M phase
ScRNA <- CellCycleScoring(ScRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
VlnPlot(ScRNA, features = c("S.Score", "G2M.Score"), group.by = "treatment", pt.size = 1,cols = col)
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

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA（after batch correction before clustering）.RData")

ScRNA <- ScRNA %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)

ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1, 
                                          by = 0.1))

#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.0.7"
ScRNA$seurat_clusters <- ScRNA@active.ident##Select your resolution based on the clustering tree
table(Idents(ScRNA))


# Display clustering, ordered by Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Generate separate umap plot
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))

DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
dev.off()

# Display tumor and normal together
pdf(paste(OUTPUT, "cluster-diff_umap.pdf"),width=6,height=6)
DimPlot(ScRNA, repel = TRUE,
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

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 100,height = 10)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18))  # Increase Y-axis title text size
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=100*dpi,h=10*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 20),  # Increase X-axis title text size
        axis.title.y = element_text(size = 20))  # Increase Y-axis title text size
dev.off()



###########Manual cell annotation########

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("/out/")
outdir <- "/out/"


output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
scedata <- readRDS(file_path)

# Define marker genes for different cell types
cellmarker <- c(
  "IL7R", "CCR7", # Naive CD4+ T
  "CD14", "LYZ",  # CD14+ Mono
  "IL7R", "CD44","S100A4","NEAT1",    # Memory CD4+ T
  "MS4A1",  # B cell
  "IGLC2","IGHA1","MZB1","JCHAIN",    # Plasma Cells
  "CD8A","CD8B", # CD8+ T
  "FCGR3A", "MS4A7", # FCGR3A+ Mono
  "GNLY", "NKG7", # NK
  "JARID2", "PLCB1",          ##NK T cells 
  "FCER1A", "CST3","LGALS2",  # DCs
  "CLEC9A","C1orf54",   # cDC1
  "ITGAM",     # cDC2
  "MYB","STMN1","CLEC4C", "CUX2","IL3RA", "LILRA4", "TCF4",  # pDCs
  "RNF220","SOX4",   # Multilymphoid progenitor
  "PPBP", # Platelets
  "KLRB1",     # Th17
  "FOXO1","FOXP3","RTKN2","IL2RA",   # Treg
  "CD25","CTLA4","TNFRSF18" , "IL10", "TGFB","GITR", "IKZF2",
  "CD3D","CD3E", "CD3G",  # T Cells
  "FOSB","LCP1","ZFP36","LRP",   # Monocytes
  "PCNA","TOP2A", "CCNA2", "CDK1"   # Proliferating Cells
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Visualize expression of immune cell marker genes using DotPlot
library(ggplot2)
plot <- DotPlot(scedata, features = unique(cellmarker))+
  theme_bw()+theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 18),  # Increase X-axis text size
    axis.text.y = element_text(size = 18),  # Increase Y-axis text size
    legend.title  = element_text(size = 18),
    legend.text = element_text(size = 16)   # Increase legend text size
  ) +
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 14, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width =14, height = 8)




library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)


col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#094867",'#E59CC4',"#6666CC",
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue → green → light green gradient
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq light colors
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
scedata <- readRDS(file_path)

# Annotate cell types for clusters
scedata <- RenameIdents(scedata, c(
  "0"="CD14+ Mono",
  "1"="CD8+ T",
  "2"="Memory CD4+ T", 
  "3"="CD14+ Mono", 
  "4"="Naive CD4+ T", 
  "5"="NK", 
  "6"="CD14+ Mono",
  "7"="CD8+ T",
  "8"="CD14+ Mono",
  "9"="B",
  "10"="Th7",
  "11"="CD14+ Mono",
  "12"="FCGR3A+ Mono",
  "13"="B",
  "14"= "Platelets", 
  "15"="CD14+ Mono",
  "16"="DCs", 
  "17"="Treg",
  "18"="NK",
  "19"="NK",
  "20"="pDCs",
  "21"="NK T",
  "22"="Plasma Cells",
  "23"="Proliferating Cells")
)

# Add cell type to metadata
scedata$celltype <- scedata@active.ident
head(scedata@meta.data)

# Extract UMAP coordinates and cell types
umap_coords <- as.data.frame(Embeddings(scedata, reduction = "umap"))
umap_coords$celltype <- scedata$celltype
#umap_coords$CB <- scedata$CB

# Save as CSV file
write.csv(umap_coords, "celltype_umap.csv", row.names = TRUE)

saveRDS(scedata,  "celltype.rds")




library(ggsci)
# Plot umap with cell types
pdf(paste(output, "ann_umap.pdf",sep = '/'), width = 7, height = 6)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=0.1,label=TRUE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.position = "none",
        #legend.title = element_text(size = 18), 
        #legend.text = element_text(size = 18),
        plot.title = element_blank())

dev.off() 

#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom"))


# Plot umap with cell types
svg(paste(output, "ann_umap.svg",sep = '/'), width = 10, height = 6)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=0.1,label=TRUE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.position = "none",
        #legend.title = element_text(size = 18), 
        #legend.text = element_text(size = 18),
        plot.title = element_blank())
#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=12,height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)
dev.off()

svg(paste(output, "ann-diff-umap.svg",sep = '/'),width=12,height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)
dev.off()



########## Add overall total cell count ##########
library(ggsci)

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

# Plot umap with cell types and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 7, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + # Center the title
  ggtitle(title_label)
dev.off()

# Plot umap with cell types and save as SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 7, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + # Center the title
  ggtitle(title_label)
dev.off()



###### Add cell labels ##########
library(ggsci)
library(Seurat)
library(dplyr)
library(ggrepel)

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

# Set cell type colors, assuming col is a predefined color vector corresponding to cell types
celltype_colors <- col  # Assuming col is a predefined color vector

# Extract UMAP coordinates
umap_coords <- Embeddings(scedata, "umap")  # Extract UMAP coordinates
umap_data <- as.data.frame(umap_coords)
umap_data$celltype <- scedata$celltype  # Add cell type information

# Create a data frame with center coordinates for each cell type
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# Calculate center coordinates for each cell type
celltype_centers <- umap_df %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    UMAP1 = mean(UMAP1, na.rm = TRUE),
    UMAP2 = mean(UMAP2, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate boundary range for each cell type to determine label placement direction
celltype_ranges <- umap_df %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    min_x = min(UMAP1, na.rm = TRUE),
    max_x = max(UMAP1, na.rm = TRUE),
    min_y = min(UMAP2, na.rm = TRUE),
    max_y = max(UMAP2, na.rm = TRUE),
    width = max_x - min_x,
    height = max_y - min_y,
    .groups = "drop"
  )

# Merge center coordinates and boundary information
celltype_labels <- left_join(celltype_centers, celltype_ranges, by = "celltype")

# Determine label placement direction: based on cell cluster shape
celltype_labels <- celltype_labels %>%
  mutate(
    # Decide label placement direction based on cluster aspect ratio
    direction_x = ifelse(width > height, 1, 0),
    direction_y = ifelse(height > width, 1, 0),
    # Set nudge parameters for labels: push labels outward from cluster
    nudge_x = ifelse(UMAP1 > mean(UMAP1), 2, -2),  # Push right clusters left, left clusters right
    nudge_y = ifelse(UMAP2 > mean(UMAP2), 2, -2)   # Push top clusters down, bottom clusters up
  )

# Plot umap with cell types and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1) +
  ggrepel::geom_text_repel(
    data = celltype_labels, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 7, 
    fontface = "bold",
    box.padding = 1.5,        # Increase border padding to keep labels away from points
    point.padding = 0.8,      # Increase point padding
    nudge_x = celltype_labels$nudge_x,  # Horizontal offset
    nudge_y = celltype_labels$nudge_y,  # Vertical offset
    min.segment.length = 1,   # Increase minimum segment length
    segment.size = 0.5,       # Segment thickness
    segment.color = "grey40", # Segment color
    segment.alpha = 0.7,      # Segment transparency
    force = 2,                # Increase repulsion force
    max.iter = 10000,         # Increase maximum iterations
    direction = "both",       # Allow bidirectional adjustment
    seed = 123,               # Set random seed for reproducibility
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  ggtitle(title_label)
dev.off()

# Plot umap with cell types and save as SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1) +
  ggrepel::geom_text_repel(
    data = celltype_labels, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 7, 
    fontface = "bold",
    box.padding = 1.5,
    point.padding = 0.8,
    nudge_x = celltype_labels$nudge_x,
    nudge_y = celltype_labels$nudge_y,
    min.segment.length = 1,
    segment.size = 0.5,
    segment.color = "grey40",
    segment.alpha = 0.7,
    force = 2,
    max.iter = 10000,
    direction = "both",
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  ggtitle(title_label)
dev.off()






##### Add total cell count for each group #####
# Count number of cells per treatment
cell_counts <- scedata@meta.data %>%
  as_tibble() %>%
  dplyr::count(treatment, name = "n") %>%
  mutate(label = sprintf("%s (n = %s cells)", treatment, format(n, big.mark=",")))


# Build named vector for replacing facet labels
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# Plot PDF
pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),
    width = 6*length(unique(scedata$treatment)), height = 5)

DimPlot(scedata, reduction = "umap", split.by = "treatment",
        pt.size = 0.1, label = FALSE, cols = col) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # Subplot title bold black
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
dev.off()

# Similarly, SVG output
svg(paste(output, "ann-diff-umap.svg",sep = '/'),
    width = 6*length(unique(scedata$treatment)), height = 5)

DimPlot(scedata, reduction = "umap", split.by = "treatment",
        pt.size = 0.1, label = FALSE, cols = col) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
dev.off()



####### Calculate cell proportions ###########

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1",'#FF9999',"#66CCCC",'#4F6272',"#FF3366","#CC0066","#00CC66","#CC99CC","#FFCCCC","#9999FF","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#094867",'#E59CC4',"#6666CC",
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue → green → light green gradient
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq light colors
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

table(scedata$seurat_clusters)

######## Calculate cell counts for all samples by cell type
cell_counts <- as.data.frame(table(Idents(scedata)))
colnames(cell_counts) <- c("CellType", "Counts")

# Sort by cell count in descending order
cell_counts <- cell_counts[order(-cell_counts$Counts), ]
# Save cell counts for all cell types to specified directory
write.csv(cell_counts, paste(output, "cell_counts.csv", sep='/'), row.names = FALSE)

# Select top 11 cell types and save to specified directory
cell_counts_top9 <- head(cell_counts, 11)
write.csv(cell_counts_top9, paste(output, "cell_counts_top9.csv", sep='/'), row.names = FALSE)

# Load required packages
library(ggplot2)
# Create bar plot
p <- ggplot(cell_counts, aes(x = reorder(CellType, -Counts), y = Counts, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Cell Type", y = "Counts", title = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 16), 
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.position = "none",  
        panel.grid = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(values = col)

# Save plot as png format
ggsave(paste(output, "cell_type_distribution.pdf", sep='/'), plot = p, width = 7, height = 6, dpi = 800)
ggsave(paste(output, "cell_type_distribution.svg", sep='/'), plot = p, width = 7, height = 6, dpi = 800)


# Calculate cell counts by cell type for each group
# Calculate number of each cell type by sample group
cell_counts_group <- as.data.frame(table(scedata$orig.ident, Idents(scedata)))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# Add grouping information (assuming grouping variable is `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # Ensure uniqueness of grouping information
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# Calculate proportion of each cell type within each sample
cell_counts_group <- cell_counts_group %>%
  group_by(Sample) %>%
  mutate(Ratio = Counts / sum(Counts))

p <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Counts') +
  scale_fill_manual(values = col) +
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # Modify X-axis labels
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size
# Add cell count text labels
p <- p + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 8, dpi = 800)


p <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Ratio') +
  scale_fill_manual(values = col) +
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # Modify X-axis labels
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size
# Add cell proportion text labels
p <- p + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 8, dpi = 800)



############Group by treatment############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate proportion of each cell type within each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

########## Create stacked bar plot for cell counts ##########
p1 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 24, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 24),  # Modify Y-axis text size
        axis.title.y = element_text(size = 26), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 24))  # Modify legend text size

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 8, dpi = 800)

########## Create stacked bar plot for cell proportions ##########
p2 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 24, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 24),  # Modify Y-axis text size
        axis.title.y = element_text(size = 26), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 24))  # Modify legend text size

file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 8, dpi = 800)


