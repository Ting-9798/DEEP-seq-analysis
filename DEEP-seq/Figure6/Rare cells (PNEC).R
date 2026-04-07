# Clean environment
rm(list = ls())

# Load required packages
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)

####################### Seurat analysis #####################
# Set output directory
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


setwd("/out(0+3+7+14)tdt+(1)/decontamination/PNEC/")
outdir <- "/out(0+3+7+14)tdt+(1)/decontamination/PNEC/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

data <- readRDS("/out(0+3+7+14)tdt+(1)/decontamination/celltype.rds")

# Select epithelial cells belonging to "Tumor", "Res" or "Sen" groups
ScRNA <- subset(data, idents = c("PNEC"))
table(ScRNA@meta.data$celltype)


#### 6. Normalization and PCA dimensionality reduction ####
# Normalization
ScRNA<-ScaleData(ScRNA)

# Run PCA
ScRNA<-RunPCA(ScRNA,npcs = 30)

pdf(paste(OUTPUT,"Dimplot.pdf"),width = 9,height = 6)
p1 <- DimPlot(object = ScRNA, reduction = "pca", pt.size = .1, group.by = "treatment",cols = col)
CombinePlots(plots=list(p1))
dev.off()

pdf(paste(OUTPUT,"vlnplot.pdf"),width = 9,height = 6)
p2 <- VlnPlot(object = ScRNA, features = "PC_1", group.by = "treatment", pt.size = 0,cols = col)
CombinePlots(plots=list(p2))
dev.off()

# PCA visualization
pdf(paste(OUTPUT, "DimHeatmap.pdf"),width = 9,height = 6)
DimHeatmap(ScRNA, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# Evaluate PC dimensions
pdf(paste0(OUTPUT,"PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(ScRNA)
dev.off()

save(ScRNA, file = "ScRNA(before_clustering).RData")



#### 7. Cell clustering and annotation ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         "#66CCCC",'#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#CCFFCC","#00CC66","#99FFFF", 
         "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#FF3300","#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA(before_clustering).RData")
# Cell clustering
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:20) %>% 
  RunTSNE(dims = 1:20) %>%
  FindNeighbors(dims = 1:20)

ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))

#Idents(ScRNA) <- "integrated_snn_res.1"
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident## Select your desired resolution based on the clustering tree
table(Idents(ScRNA))

# Ensure "treatment" factor levels are ordered as Non-infected and Infected
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Visualize clustering, displayed in order of Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 8, height = 4)
DimPlot(ScRNA, reduction = "umap", pt.size=2,label = TRUE, repel = TRUE, split.by = "treatment", label.size = 5,cols = col)+
  theme(
    strip.text = element_text(size = 22, face = "bold"),  # Increase subplot title font size
    axis.text.x = element_text(size = 16),  # X-axis label size
    axis.text.y = element_text(size = 16),  # Y-axis label size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
    axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18)    # Increase legend text size
  )
dev.off()

# Generate UMAP plot separately
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 5.5, height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=2,label = TRUE, repel = TRUE, label.size = 6,cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
        axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
        legend.title = element_text(size = 22),  # Increase legend title size
        legend.text = element_text(size = 22))
dev.off()

pdf(paste(OUTPUT, "cluster_umap1.pdf"), width = 6, height = 4)
DimPlot(ScRNA, reduction = "umap", pt.size=2,label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
        axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
        legend.title = element_text(size = 22),  # Increase legend title size
        legend.text = element_text(size = 22))
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

saveRDS(ScRNA, "ScRNA(after_clustering).rds")




########## Check expression proportions of "tdTomato" and "Epcam" ###########
##### Merge plotting ######
genes <- c("tdTomato", "Epcam","Uchl1","Resp18")
subset_data <- ScRNA

for (gene in genes) {
  # Get gene expression data
  gene_expr <- FetchData(subset_data, vars = gene)
  subset_data[[paste0(gene, "_expr")]] <- gene_expr[[gene]]
  
  # Calculate proportion of cells with non-zero expression
  expressed_cells <- sum(subset_data[[paste0(gene, "_expr")]] > 0)
  total_cells <- nrow(subset_data@meta.data)
  expression_ratio <- expressed_cells / total_cells * 100
  
  # Extract expression vector from metadata
  expr_vec <- subset_data@meta.data[[paste0(gene, "_expr")]]
  
  # Order cells (low expression at bottom)
  cells_ordered <- subset_data@meta.data[order(expr_vec, decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)
  
  # Set title
  plot_title <- paste0(gene, " (", round(expression_ratio, 2), "%) ")
  
  # PDF output path
  pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge", ".pdf")
  svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge", ".svg")
  
  # PDF plot
  pdf(pdf_path, width = 4, height = 4)
  print(
    FeaturePlot(
      subset_data,
      features = gene,
      reduction = "umap",
      pt.size=2,
      cells = cell_names_ordered,
      ncol = 1,
      cols = c('#E5D2DD',  "#FF3366")
    ) +
      ggtitle(plot_title) +
      theme(
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      ) +
      NoAxes()
  )
  dev.off()
  
  # SVG plot
  svg(svg_path, width = 4, height = 4)
  print(
    FeaturePlot(
      subset_data,
      features = gene,
      reduction = "umap",
      pt.size=2,
      cells = cell_names_ordered,
      ncol = 1,
      cols = c('#E5D2DD',  "#FF3366")
    ) +
      ggtitle(plot_title) +
      theme(
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      ) +
      NoAxes()
  )
  dev.off()
}


# Groups and genes for plotting
treatment_groups <- c("0d", "3d","7d", "14d")
genes <- c("tdTomato", "Epcam","Uchl1","Resp18")

# Loop through each treatment and gene for plotting
for (treat in treatment_groups) {
  # Subset cells for corresponding treatment
  subset_data <- subset(ScRNA, subset = treatment == treat)
  
  for (gene in genes) {
    # Get gene expression data
    gene_expr <- FetchData(subset_data, vars = gene)
    subset_data[[paste0(gene, "_expr")]] <- gene_expr[[gene]]
    
    # Calculate proportion of cells with non-zero expression
    expressed_cells <- sum(subset_data[[paste0(gene, "_expr")]] > 0)
    total_cells <- nrow(subset_data@meta.data)
    expression_ratio <- expressed_cells / total_cells * 100
    
    # Extract expression vector from metadata
    expr_vec <- subset_data@meta.data[[paste0(gene, "_expr")]]
    
    # Order cells (low expression at bottom)
    cells_ordered <- subset_data@meta.data[order(expr_vec, decreasing = FALSE), ]
    cell_names_ordered <- rownames(cells_ordered)
    
    # Set title
    plot_title <- paste0(gene, " (", round(expression_ratio, 2), "%) - ", treat)
    
    # PDF output path
    pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, ".pdf")
    svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, ".svg")
    
    # PDF plot
    pdf(pdf_path, width = 4, height = 4)
    print(
      FeaturePlot(
        subset_data,
        features = gene,
        reduction = "umap",
        pt.size=2,
        cells = cell_names_ordered,
        ncol = 1,
        cols = c('#E5D2DD',  "#FF3366")
      ) +
        ggtitle(plot_title) +
        theme(
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, size = 1)
        ) +
        NoAxes()
    )
    dev.off()
    
    # SVG plot
    svg(svg_path, width = 4, height = 4)
    print(
      FeaturePlot(
        subset_data,
        features = gene,
        reduction = "umap",
        pt.size=2,
        cells = cell_names_ordered,
        ncol = 1,
        cols = c('#E5D2DD',  "#FF3366")
      ) +
        ggtitle(plot_title) +
        theme(
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, size = 1)
        ) +
        NoAxes()
    )
    dev.off()
  }
}





# Display expression trends of marker genes for the first 6 clusters
output <- paste(outdir,'cell_localization', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA(after_clustering).rds")
ScRNA <- readRDS(file_path)

# Find marker genes
ScRNA.markers <- FindAllMarkers(ScRNA, only.pos = TRUE,   ### only.pos = TRUE: only find upregulated genes
                                min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ScRNA.markers,paste0(output,"./ScRNA.all.markers.csv"))

dim.use <-1:30
top5 <- ScRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top5,file=paste0(output,"/top20_marker_genes_tsne_",max(dim.use),"PC.csv"))

pdf(paste0(output,"/Heatmap_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 25,height = 20)
DoHeatmap(ScRNA, features = top5$gene,size = 5)+
  scale_fill_gradientn(colors = c("#437eb8", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18),  # Increase Y-axis title text size
        legend.text = element_text(size = 14),  # Adjust legend text size
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 30))  # Adjust grouping label color and size
dev.off()

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 70,height = 10)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18))  # Increase Y-axis title text size
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=70*dpi,h=10*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 20),  # Increase X-axis title text size
        axis.title.y = element_text(size = 20))  # Increase Y-axis title text size
dev.off()






######### Expression trends of marker genes in cells #######
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# Select differential clusters for display
output <- paste(outdir,'cluster(Club+AT2+Cilliated)', sep='/')
dir.create(output)

#file_path <- file.path(outdir, "celltype.rds")
file_path <- file.path(outdir, "ScRNA(after_clustering).rds")
ScRNA <- readRDS(file_path)


# Set T cell activation related genes
cellmarker <- c(
  "Ccn1","Ccn2","Resp18","Cbr2","Hes1","Krt8" ,"Foxj1","Scgb3a2","Lamp3","Ager"    #PNEC (pulmonary neuroendocrine cells)
  # "Cbr2","Notch2","Hes1","Myb","Calca","Mki67",
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# Create lists to store all RidgePlot, VlnPlot and FeaturePlot
ridge_plots <- list()
vln_plots <- list()
feature_plots <- list()

# Loop through immune cell marker gene list
for (gene in cellmarker) {
  # RidgePlot
  ridge_plots[[gene]] <- RidgePlot(ScRNA, features = gene, ncol = 1, cols = col) +
    theme(legend.position = "none")  
  
  # VlnPlot
  vln_plots[[gene]] <- VlnPlot(ScRNA, features = gene, ncol = 1, pt.size = 0, cols = col) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") 
  
  # Get gene expression data, ensure return as numeric vector
  gene_expr <- as.numeric(FetchData(ScRNA, vars = gene, assay = "RNA")[[gene]])  # Convert to numeric vector
  
  # Check if gene expression was successfully extracted
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # Add gene expression data to metadata
  ScRNA@meta.data[[paste0(gene, "_expr")]] <- gene_expr  # Add expression information to metadata
  
  # Use numeric vector for ordering, avoid directly indexing ScRNA object
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # Extract ordered cell names
  
  # FeaturePlot drawn in ordered cell sequence, using continuous color gradient
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", pt.size=3,
                                       cells = cell_names_ordered,  # Specify cell order
                                       ncol = 1) + 
    scale_color_gradientn(colors = c("#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366")) +  # Set continuous color gradient
    #  "#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366"
    theme(legend.position = "right", 
          plot.title = element_text(size = 22, face = "bold"),  # Increase title text size and make bold
          legend.text = element_text(size = 18)) +  # Increase legend text size
    NoAxes()  # Remove axes
  
}


# Save RidgePlot
#pdf(paste0(out, "cellmarker_RidgePlot.pdf"), width = 25, height = 12)
#print(cowplot::plot_grid(plotlist = ridge_plots, ncol = 4))
#dev.off()

# Save FeaturePlot
pdf(paste0(output, "/spacial_FeaturePlot_umap.pdf"), width = 15, height = 12)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

svg(paste0(output, "/spacial_FeaturePlot_umap.svg"), width = 15, height = 12)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

pdf(paste0(output, "/spacial_VlnPlot_umap.pdf"), width = 12, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()
svg(paste0(output, "/spacial_VlnPlot_umap.svg"), width = 12, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()




################ Batch differential analysis ################
library(scRNAtoolVis)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggrepel)
library(org.Mm.eg.db) # Mouse database
#library(org.Hs.eg.db) # Human database
library(clusterProfiler)
library(enrichplot)
library(DOSE)


col <- c('#437eb8','#FF6666',"#FFFFCC",'#FFCC99','#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFCCCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")


# Create output directory
output <- file.path(outdir, "differential_analysis")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# Read data
file_path <- file.path(outdir, "ScRNA(after_clustering).rds")
scRNAsub <- readRDS(file_path)

logFCfilter <- 0.25
adjPvalFilter <- 0.05

genes <- c("Gpnmb","Spp1","Ctsd","Nfkb1","Ifngr1","Camk4","Zeb1","Rora","Cd14","Il1b",
           "Il12rb2","Cxcl2","Fth1","Icos","Stat4","Cd63","Thbs1","Tyrobp")

timepoints <- c("7d","14d")

all_markers <- list()

for (tp in timepoints) {
  
  
  comp_name <- paste0(tp, " vs 3d")
  
  markers <- FindMarkers(object = scRNAsub,
                         ident.1 = tp,
                         ident.2 = "3d",
                         group.by = "treatment",
                         logfc.threshold = 0,
                         min.pct = 0.25,
                         test.use = "wilcox")
  
  markers$gene <- rownames(markers)
  markers <- markers %>%
    mutate(Significance = ifelse(p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter, 
                                 ifelse(avg_log2FC > 0, "Up", "Down"), "Normal")) %>%
    mutate(Group = comp_name) 
  
  write.table(markers, file = file.path(output, paste0("sig.markers_", comp_name, ".txt")),
              sep = "\t", row.names = TRUE, quote = FALSE)
  saveRDS(markers, file = file.path(output, paste0("ScRNA.sig.markers_", comp_name, ".rds")))
  
  up_df <- markers %>% filter(Significance == "Up")
  down_df <- markers %>% filter(Significance == "Down")
  write.csv(up_df, file = file.path(output, paste0("upregulated_genes_", comp_name, ".csv")), row.names = TRUE)
  write.csv(down_df, file = file.path(output, paste0("downregulated_genes_", comp_name, ".csv")), row.names = TRUE)
  
  
  # Calculate number of upregulated and downregulated genes
  upregulated_genes <- sum(markers$Significance == "Up")
  downregulated_genes <- sum(markers$Significance == "Down")
  total_diff_genes <- upregulated_genes + downregulated_genes
  
  # Save data frames for upregulated and downregulated genes separately
  upregulated_genes_df <- markers %>%
    filter(Significance == "Up")
  downregulated_genes_df <- markers %>%
    filter(Significance == "Down")
  
  # Select top 10 upregulated and downregulated genes to display labels
  top_genes_upregulated <- upregulated_genes_df %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(p_val_adj) 
  top_genes_downregulated <- downregulated_genes_df %>%
    filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
    arrange(p_val_adj) %>%
    head(10)
  
  # ---- Volcano Plot ----
  interested_genes <- markers %>% filter(gene %in% genes)
  
  p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Significance), size = 2, shape = 18) +
    scale_color_manual(values = c("Up" = "#ca0020", "Down" = "#0099CC", "Normal" = "#FFCCCC")) +
    geom_hline(yintercept = -log10(adjPvalFilter), linetype = "dashed") +
    geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed") +
    geom_text_repel(data = top_genes_upregulated, aes(label = top_genes_upregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
    #geom_text_repel(data = interested_genes, aes(label = gene),size = 5, fontface = "bold", box.padding = 0.6, max.overlaps = 50) +
    theme_classic() +
    labs(title = comp_name, x = "log2 Fold Change", y = "-log10 Adjusted P-value", color = "Significance") +
    scale_x_continuous(limits = c(-1.5, 2.5), breaks = seq(-1.5, 2.5, by = 1)) +
    scale_y_continuous(limits = c(0, max(-log10(markers$p_val_adj), na.rm = TRUE) + 1))+
    #scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 5)) +
    theme(plot.title = element_text(size = 22, face = "bold", hjust = 0), 
          legend.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 20, hjust = 0.5),
          axis.text = element_text(size = 18))
  
  
  ggsave(file.path(output, paste0(comp_name, "_volcano_plot.svg")), p, width = 8, height = 7)
  ggsave(file.path(output, paste0(comp_name, "_volcano_plot.pdf")), p, width = 8, height = 7)
  
  
  all_markers[[tp]] <- markers
}


# Merge three groups
merged_markers <- bind_rows(all_markers)

# Save multi-group input table
write.csv(merged_markers, file = file.path(output, "merged_markers_for_mutiVolcano.csv"), row.names = FALSE)


######### Draw multi-group differential volcano plot ###############
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)
library(ggsci)
library(patchwork)
library(tidyr)

# Load data
merged_markers <- read.csv(file.path(output, "merged_markers_for_mutiVolcano.csv"))

# Data preprocessing
data <- merged_markers %>%
  rename(logFC = avg_log2FC, adj.P.Val = p_val_adj) %>%
  mutate(change = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.25,
                         ifelse(logFC > 0, "Up", "Down"),
                         "No change")) %>%
  filter(change != "No change") %>%
  mutate(group = Group) %>%
  mutate(label = change)

# Set group order
group_order <- c("7d vs 3d", "14d vs 3d")
data <- data %>%
  mutate(group = factor(group, levels = group_order))

# Extract Top20 from Up and Down in each group separately
TopGene_up <- data %>%
  filter(change == "Up") %>%
  group_by(group) %>%
  distinct(gene, .keep_all = TRUE) %>%
  top_n(20, wt = logFC)

TopGene_down <- data %>%
  filter(change == "Down") %>%
  group_by(group) %>%
  distinct(gene, .keep_all = TRUE) %>%
  top_n(20, wt = -logFC)

TopGene <- bind_rows(TopGene_up, TopGene_down) %>%
  ungroup()

# Background bar plot data
dbar <- data %>%
  group_by(group) %>%
  summarise(logFC_min = min(logFC),
            logFC_max = max(logFC))

# Plot
p <- ggplot() +
  # Background bar plot
  geom_col(data = dbar, aes(x = group, y = logFC_min), fill = "#dcdcdc", alpha = 0.6, width = 0.7) +
  geom_col(data = dbar, aes(x = group, y = logFC_max), fill = "#dcdcdc", alpha = 0.6, width = 0.7) +
  
  # All points
  #geom_jitter(data = TopGene, aes(x = group, y = logFC, color = label), size = 0.85, width = 0.3) +
  
  # Top gene points
  geom_jitter(data = TopGene, aes(x = group, y = logFC, color = label), size = 1.5, width = 0.35) +
  
  # Middle tile
  geom_tile(data = TopGene, aes(x = group, y = 0, fill = group), height = 0.6, color = "black", alpha = 0.6, show.legend = FALSE) +
  
  # Top gene labels, label position offset by 0.1 towards logFC direction to prevent overlap
  geom_text_repel(data = TopGene,
                  aes(x = group, y = logFC - 0.1 * sign(logFC), label = gene),
                  size = 5, color = 'black',
                  force = 1.2,
                  max.overlaps = 50,
                  arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
  
  # Center group labels
  geom_text(data = TopGene, aes(x = group, y = 0, label = group), size = 5, color = "white") +
  
  ggsci::scale_fill_npg() +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5"))+
  labs(x="",y = "log2 Fold Change", color = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 18)
  )

# Save images
ggsave(filename = file.path(output, "multi_volcano_plot.pdf"), plot = p, width = 6, height = 5, bg = "white")
ggsave(filename = file.path(output, "multi_volcano_plot.svg"), plot = p, width = 6, height = 5, bg = "white")




# ======== Functional analysis: Batch GSEA, GO, KEGG analysis ========
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)

# Loop through each differential analysis result
for (tp in timepoints) {
  comp_name <- paste0(tp, " vs 3d")
  marker_file <- file.path(output, paste0("ScRNA.sig.markers_", comp_name, ".rds"))
  markers <- readRDS(marker_file)
  
  message(">>> Processing: ", comp_name)
  
  
  # ---- GO & KEGG enrichment analysis ----
  gene_up <- markers$gene[markers$Significance == "Up"]
  gene_down <- markers$gene[markers$Significance == "Down"]
  
  gene_up_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                               keys = gene_up, 
                                                               columns = 'ENTREZID', 
                                                               keytype = 'SYMBOL')[,2]))
  gene_down_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                                 keys = gene_down, 
                                                                 columns = 'ENTREZID', 
                                                                 keytype = 'SYMBOL')[,2]))
  
  # GO enrichment analysis
  go_up <- enrichGO(gene = gene_up_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.1)
  go_down <- enrichGO(gene = gene_down_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.1)
  
  # Remove species suffix from Description in GO analysis results
  go_up@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", go_up@result$Description)
  go_down@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", go_down@result$Description)
  
  # Convert geneID from ENTREZID to SYMBOL
  go_up@result$geneID <- sapply(strsplit(go_up@result$geneID, "/"), function(ids) {
    symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
    paste(symbols, collapse = "/")
  })
  
  go_down@result$geneID <- sapply(strsplit(go_down@result$geneID, "/"), function(ids) {
    symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
    paste(symbols, collapse = "/")
  })
  
  write.csv(as.data.frame(go_up), file = file.path(output, paste0(comp_name, "_GO_UP.csv")))
  write.csv(as.data.frame(go_down), file = file.path(output, paste0(comp_name, "_GO_DOWN.csv")))
  
  if (nrow(go_up) > 0) {
    gop_up <- dotplot(go_up) + ggtitle(paste(comp_name, "GO UP"))
    ggsave(file.path(output, paste0(comp_name, "_GO_UP.pdf")), gop_up, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_GO_UP.svg")), gop_up, width = 6, height = 6)
  }
  if (nrow(go_down) > 0) {
    gop_down <- dotplot(go_down) + ggtitle(paste(comp_name, "GO DOWN"))
    ggsave(file.path(output, paste0(comp_name, "_GO_DOWN.pdf")), gop_down, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_GO_DOWN.svg")), gop_down, width = 6, height = 6)
  }
  
  # KEGG enrichment analysis
  kegg_up <- enrichKEGG(gene = gene_up_entrez, organism = "mmu", pvalueCutoff = 0.1)
  kegg_down <- enrichKEGG(gene = gene_down_entrez, organism = "mmu", pvalueCutoff = 0.1)
  
  # Remove species suffix from Description in KEGG analysis results
  kegg_up@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_up@result$Description)
  kegg_down@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_down@result$Description)
  
  # Convert geneID from ENTREZID to SYMBOL
  kegg_up@result$geneID <- sapply(strsplit(kegg_up@result$geneID, "/"), function(ids) {
    symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
    paste(symbols, collapse = "/")
  })
  
  kegg_down@result$geneID <- sapply(strsplit(kegg_down@result$geneID, "/"), function(ids) {
    symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
    paste(symbols, collapse = "/")
  })
  
  write.csv(as.data.frame(kegg_up), file = file.path(output, paste0(comp_name, "_KEGG_UP.csv")))
  write.csv(as.data.frame(kegg_down), file = file.path(output, paste0(comp_name, "_KEGG_DOWN.csv")))
  
  if (nrow(kegg_up) > 0) {
    kp_up <- dotplot(kegg_up) + ggtitle(paste(comp_name, "KEGG UP"))
    ggsave(file.path(output, paste0(comp_name, "_KEGG_UP.pdf")), kp_up, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_KEGG_UP.svg")), kp_up, width = 6, height = 6)
  }
  if (nrow(kegg_down) > 0) {
    kp_down <- dotplot(kegg_down) + ggtitle(paste(comp_name, "KEGG DOWN"))
    ggsave(file.path(output, paste0(comp_name, "_KEGG_DOWN.pdf")), kp_down, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_KEGG_DOWN.svg")), kp_down, width = 6, height = 6)
  }
}






######### Draw pathways of interest ###############
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)
library(ggsci)
library(patchwork)
library(tidyr)

# Load data
go_up <- read.csv(file.path(output, "14d vs 3d_GO_UP.csv"))


# Prepare visualization data
# Extract GO analysis results
go_up_dt <- as.data.frame(go_up)
#go_down_dt <- as.data.frame(go_down)


# Extract GO pathways of interest
interested_go <- c(
  "tight junction assembly","tight junction organization","bicellular tight junction assembly","transepithelial transport",
  "epithelial fluid transport","cell-cell junction assembly","lung alveolus development","epithelial cell proliferation",
  "collagen biosynthetic process","collagen metabolic process","regulation of tissue remodeling",
  "regulation of bone remodeling","prostaglandin biosynthetic process","prostaglandin metabolic process",
  "icosanoid biosynthetic process","icosanoid metabolic process","endothelial cell proliferation",
  "vascular transport","chemokine production",
  "interleukin-6 production","interleukin-8 production","tumor necrosis factor-mediated signaling pathway",
  "response to interferon-alpha","response to type II interferon",
  "JAK-STAT cascade","NF-kappaB signaling","ERK1 and ERK2 cascade","p38MAPK cascade",
  
  "DNA damage response, signal transduction by p53 class mediator","signal transduction by p53 class mediator",
  "signal transduction in response to DNA damage"
)

# Filter pathways of interest from GO enrichment results
go_up_dt <- go_up_dt[go_up_dt$Description %in% interested_go, ]


# Set color classification
classification_colors <- c('#437eb8','#FF6666','#FFCC99','#FF9999', '#80c5d8',"#9999FF",
                           "#FFCCCC","#99CCFF","#FF3366","#CCCCFF","#CC0066","#FFFFCC",
                           "#66CCCC","#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                           "#6699CC","#CC99CC","#FF6699","#FF0000","#6666CC","#FF9966",
                           "#669999","#CC99FF","#FFCCFF",
                           '#437eb8','#FF6666','#FFCC99','#FF9999', '#80c5d8',"#9999FF",
                           "#FFCCCC","#99CCFF","#FF3366","#CCCCFF","#CC0066","#FFFFCC",
                           "#66CCCC","#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                           "#6699CC","#CC99CC","#FF6699","#FF0000","#6666CC","#FF9966",
                           "#669999","#CC99FF","#FFCCFF")

# Text wrapping function to limit character length
wrap_text <- function(text, width = 40) {
  sapply(text, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

# Bar plot for GO enrichment analysis
plot_GO_bar <- function(dt, title) {
  dt <- dt[order(-dt$pvalue, decreasing = TRUE),]  # Sort by p-value first
  #dt <- dt[order(dt$Count, decreasing = TRUE),]  # Sort by gene count first
  
  #dt <- head(dt,20)  # Select top 20 pathways
  dt$Description <- factor(wrap_text(dt$Description), levels = wrap_text(dt$Description))
  
  # Left plot: enrichment p-value
  p1 <- ggplot(dt, aes(x = Description, y = log10(p.adjust), fill = Description)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = classification_colors) +
    coord_flip() +
    ylab('-log10(P-value)') +
    xlab('') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          plot.title =  element_text(size = 16, face = "bold"),
          legend.position = "none",
          plot.margin = margin(10, 10, 10, 10),
          panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    ggtitle(title) 
  
  # Right plot: gene count
  p2 <- ggplot(dt, aes(x = Description, y = Count)) +
    geom_bar(stat = 'identity', fill = '#66CCCC') +
    coord_flip() +
    ylab('Gene Count') +
    xlab('') +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(size = 14, face = "bold"),
          plot.margin = margin(10, 10, 10, 10),
          panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add XY axes borders
  
  # Combine two plots
  p_combined <- p1 + p2 + plot_layout(widths = c(2, 1.5))
  return(p_combined)
}

# Draw bar plot for GO and KEGG enrichment analysis
go_up_plot <- plot_GO_bar(go_up_dt, "Upregulated Genes GO Enrichment")

# Save image
ggsave(file.path(output, 'go_enrich_up_bar.pdf'), plot = go_up_plot, width = 8, height = 6)





############ Batch draw pathway gene expression heatmaps ############
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# Load data
ScRNA <- readRDS("ScRNA(after_clustering).rds")

# Set group order
ScRNA$treatment <- factor(ScRNA$treatment, levels = c("3d", "7d", "14d"))

# Set gene lists for multiple pathways
pathway_genes <- list(
  ## p53 signaling pathway
  p53 = c("Trp53","Mdm2","Ccnd1","Ccnd2","Ccnd3","Ccne1","Ccne2",
          "Ccng1","Ccng2","Cdk1","Cdk2","Cdk4","Cdk6","Cdkn1a","Gadd45a","Gadd45g",
          "Ei24","Ppm1d","Bax","Bbc3","Bid","Apaf1","Ddit4","Pig3",
          "Casp3","Casp8"),
  
  ## Wnt signaling pathway
  Wnt = c("Axin1","Axin2","Ctnnb1","Gsk3β","Wnt3","Wnt5a","Wnt7b",
          "Lrp5","Lrp6","WISP1","Tnc","C-MYC","Cyclin D","VEGF"),
  
  ## Hippo signaling pathway
  Hippo = c("Mst1","Mst2","Stk3","Stk4","Lats1","Lats2","Yap1","Taz","Ajuba",
            "Wwtr1","Ctgf","Serpine1","Tead1","Tead2","Tead3","Tead4","Mob1","Sav1"),
  
  ## JAK-STAT signaling pathway
  JAK_STAT = c("Jak2","Jak3","Tyk2","Stat1","Stat3","Stat5","Stat6","Il11","Socs1","Socs3"),
  
  ## PI3K-Akt signaling pathway
  PI3K_Akt = c("Pik3ca","Pik3cb","Pik3r1","Akt1","Akt2","Akt3","Mtor","Rictor","Rheb",
               "Tsc1","Tsc2","Pten","Sgk1","Fgf2","Igf1","Igf1r","Cd274","Src","Dab2",
               "Hmox1","Nqo1","Nfe2l2"),
  
  spacial = c("Piezo1","Piezo2","Osr2","Egr3"),
  
  Hedgehog = c("Shh","Ptch1","Ptch2","Smo","Sufu","Kif7","Gli3","Gas1","Boc","Cdon","Hhip"),
  
  Tgfb = c("Tgfb1","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2","Tgfbr3","Acvr1","Acvr2a","Acvr2b","Bmpr1a","Bmpr1b","Bmpr2","Smad2","Smad3","Smad4","Smad1","Smad5","Smad9","Smad6","Smad7","Map3k7"),
  
  EGF = c("Tgfa","Hbegf","Areg","Btc","Nrg1","Nrg2","Nrg3","Nrg4","Egfr","Erbb2","Erbb3","Erbb4","Shc1","Grb2","Sos1","Sos2","Gab1","Kras","Hras","Nras","Raf1","Braf","Map2k1","Map2k2","Mapk1","Mapk3","Pik3ca","Pik3cb","Pik3r1","Akt1","Akt2","Akt3","Mtor","Plcg1","Jak2","Stat3","Stat5","Pten","Cbl"),
  
  FGF = c("Fgf1","Fgf2","Fgf7","Fgf9","Fgf18","Fgf21","Fgfr1","Fgfr2","Fgfr3","Frs2","Grb2","Sos1","Gab1","Pik3ca","Pik3cb","Pik3r1","Akt1","Akt2","Akt3","Plcg1","Mapk1","Mapk3","Klb","Kl","Sdc1","Sdc2","Sdc3","Sdc4"),
  
  Rb = c("Rb1","Rbl1","Rbl2","E2f1","E2f2","E2f3","E2f4","E2f5","E2f6","E2f7","Ccnd1","Ccnd2","Ccnd3","Cdk4","Cdk6","Ccne1","Ccne2","Cdk2","Cdkn2b","Cdkn1a","Cdkn1b")
  
  
  
)

# Loop to draw heatmap for each pathway
for (pathway_name in names(pathway_genes)) {
  
  genes <- pathway_genes[[pathway_name]]
  
  # Extract expression matrix and calculate average expression for each treatment group
  avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment", assays = "RNA")$RNA
  
  avg_expr_mat <- as.matrix(avg_expr)
  rownames(avg_expr_mat) <- rownames(avg_expr)
  
  # Z-score transformation
  mat_scaled <- t(scale(t(avg_expr_mat)))
  
  # Group information
  group <- colnames(avg_expr_mat)  # treatment
  group_anno <- HeatmapAnnotation(
    Treatment = factor(group, levels = c("0d", "3d", "7d", "14d")),
    col = list(Treatment = c( "3d" = "#66CCCC", "7d" = "#FF9933", "14d" = "#CC0066"))   # "0d" = "#99CCFF",
  )
  
  # Heatmap color gradient
  heatmap_col <- colorRampPalette(c('#3399CC', "white", "#FF3366"))(100)
  
  # Save heatmap
  pdf(file.path(output, paste0("Heatmap_", pathway_name, "_byTreatment.pdf")), width = 6, height = 5)
  Heatmap(mat_scaled,
          name = "Z-score",
          top_annotation = group_anno,
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 12),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          col = heatmap_col,
          heatmap_legend_param = list(
            title = "Z-score",
            title_gp = gpar(fontsize = 12),
            labels_gp = gpar(fontsize = 10)
          )) %>% draw()
  dev.off()
}


##################### Pseudotime analysis monocle2 ###############################
library(Seurat)
library(monocle)
library(igraph)
#devtools::install_version("igraph", version = "2.0.8", repos = "http://cran.us.r-project.org")

#BiocManager::install("monocle")
#install.packages("igraph")
dpi=300

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

output <- paste(outdir,'monocle2', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA(after_clustering).rds")
data1 <- readRDS(file_path)
summary(data1$treatment)

# Remove 0d
data <- subset(data1, subset = treatment != "0d")
summary(data$treatment)

#data <- subset(data, subset = treatment == "Tumor")
#View(data@meta.data)

## Extract original expression matrix and sparsify: UMI count
expr_matrix<-as(as.matrix(data@assays$RNA@data), 'sparseMatrix')
print(head(expr_matrix[,1:4]))

## Extract phenotype information, i.e., cell information
p_data<-data@meta.data
rownames(p_data)<-colnames(expr_matrix)
head(data@active.ident)

# Keep only M1 and M2 cells
#p_data <- p_data[p_data$celltype %in% c("M1", "M2"), ]
#expr_matrix <- expr_matrix[, rownames(p_data)]

## Extract gene information
f_data<-data.frame(gene_id=rownames(expr_matrix),gene_short_name=rownames(expr_matrix))
rownames(f_data)<-rownames(expr_matrix)
print(head(f_data))

## Construct monocle2 object
fd<-new("AnnotatedDataFrame", data = f_data)
pd<-new("AnnotatedDataFrame", data =p_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())

# Estimate library size and dispersion - normalization
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## Ordering gene selection and visualization
## Calculate number of cells expressing each gene
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

## Select genes expressed in 5 or more cells to speed up calculation
express_genes <- row.names(subset(fData(cds),num_cells_expressed>=5)) #subset: take subset
str(express_genes)
head(fData(cds))

## Perform differential analysis on remaining genes, differential analysis based on different celltypes, cores=2 means using 2 CPUs
diff <- differentialGeneTest(cds[express_genes,],fullModelFormulaStr="~celltype",cores=2) 
head(diff)

## Differentially expressed genes as genes for trajectory construction, selection criteria is qval<0.01, decreasing=F means sort by increasing value
deg <- subset(diff)  #qval < 0.01
head(deg)

write.csv(deg,file=paste0(output,"/monocle.DEG.csv"),row.names=FALSE)


## Trajectory construction genes (ordering genes) visualization
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
pdf(paste0(output,"/ordergenes.pdf"))
plot_ordering_genes(cds)
dev.off()
ggsave(paste0(output,"/ordergenes.png"),plot_ordering_genes(cds),width = dpi*6, height = dpi*6, units = "px",type='cairo')

## Dimensionality reduction
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)


######## Trajectory visualization
## Pseudotime indicates pseudotime value, State indicates cell state, celltype indicates cell type
## Trajectory visualization by different types
types <- c("Pseudotime", "State", "celltype", "treatment")

# Custom colors
custom_colors <- c(
  "Tumor" = '#FF6666',       # Red
  "Normal" = '#E5D2DD',     # Green
  "DCs-cancer" = '#FF6666',  # Blue
  "DCs-no_cancer" = '#E5D2DD' # Orange
)


for (type in types) {
  if (type == "Pseudotime") {
    # For continuous data (Pseudotime), use gradient color
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # Set gradient color
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # For discrete data (State), use color vector col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # Use col as colors
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # For discrete data (celltype and treatment), use custom colors custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # Use custom_colors as colors
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  # Save as PDF format
  ggsave(filename = paste(output, paste0("monocle_", type, ".pdf", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 10, height = 5)
  
  # Save as SVG format
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 10, height = 5)
}


#saveRDS(cds,  file = file.path(output, "monocle2.rds"))
saveRDS(cds,"monocle2.rds")


##### Calculate cell proportions after pseudotime #####

# Extract state information for each cell
state_data <- pData(cds)$State
celltype_data <- pData(cds)$celltype

# Calculate counts of different cell types under each state
cell_counts_state <- as.data.frame(table(celltype_data, state_data))
colnames(cell_counts_state) <- c("CellType", "State", "Counts")

# Calculate proportion of different cell types under each state
cell_counts_state$Ratio <- ave(cell_counts_state$Counts, cell_counts_state$State, FUN = function(x) x / sum(x))

# Draw bar plot of cell type proportions under each state
p <- ggplot(cell_counts_state, aes(x = State, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='State', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.title = element_blank(),         
        legend.text = element_text(size = 12))  

# Save image
ggsave(paste0(output, "/state_celltype_proportion.pdf"), plot = p, width = 7)
ggsave(paste0(output, "/state_celltype_proportion.svg"), plot = p, width = 7)



file_path <- file.path(outdir, "monocle2.rds")
cds <- readRDS(file_path)

## Confirm root node
cds <- orderCells(cds,root_state=6)
### Re-execute visualization, pseudotime direction changes
types <- c("Pseudotime","State","celltype","treatment")

# Custom colors
custom_colors <- c(
  "Tumor" = '#FF6666',       # Red
  "Normal" = '#E5D2DD',     # Green
  "BC-cancer" = '#FF6666',  # Blue
  "BC-no_cancer" = '#E5D2DD' # Orange
)


for (type in types) {
  if (type == "Pseudotime") {
    # For continuous data (Pseudotime), use gradient color
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # Set gradient color
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # For discrete data (State), use color vector col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # Use col as colors
      theme(legend.text = element_text(size = 14),  # Adjust legend text size
            legend.title = element_text(size = 16),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # For discrete data (celltype and treatment), use custom colors custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # Use custom_colors as colors
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  # Save as PDF format
  ggsave(filename = paste(output, paste0("monocle_", type, ".pdf", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 6, height = 4)
  
  # Save as SVG format
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 6, height = 4)
}


## Color by cell state (split)
# Generate split state trajectory plot
plot_cell_traj_facet <- plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap("~State", nrow = 1) +
  scale_color_manual(values =col)+
  theme(legend.text = element_text(size = 14),  # Adjust legend text size
        legend.title = element_text(size = 14))  # Adjust legend title size

ggsave(filename = paste0(output, "/monocle_state_facet.pdf"), plot = plot_cell_traj_facet, width=10, height=6)
ggsave(filename = paste0(output, "/monocle_state_facet.svg"), plot = plot_cell_traj_facet, width=10, height=6)

## Heatmap based on expression trends
topgene <- ordergene[1:50]
topgene <- c(
  "Foxj1","Tppp3","Tubb4b","Tubb1",   #Cilliated cells,
  "Scgb1a1","Scgb3a2","Chad",     #club cell
  "Cd14","Cd74","H2-K1",    # Activated Club
  "H2-Ab1","Cst3",     # MHC-II+ Club
  'Ager', 'Hopx', 'Rtkn2', 'Aqp5',"Cav1 ","Spock2",   #AT1
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  "Lrg1","Lcn2","Retnla","Il33","Car8","Ank3","Cftr",         # Activated AT2
  "Birc5","Top2a",       # Proliferating AT2s
  "Cldn4","Sfn","Clu","Krt19","Krt8"    # PATs [pre-AT1 transitional state] Epithelial transitional states 
  
)
plot_pseu_heatmap <- plot_pseudotime_heatmap(cds[topgene,],num_clusters = 5,cores = 1,
                                             show_rownames = T,return_heatmap=T,hmcols = colorRampPalette(c("#1f77b4", "#ffffff", "#FF3366"))(100))
pdf(paste0(output,"/monocle_pheatmap.pdf"))
print(plot_pseu_heatmap)
dev.off()
ggsave(paste0(output,"/monocle_pheatmap.svg"),plot_pseu_heatmap,width = 6, height = 7)



## Expression change plots of key driver genes
# Select first 4 genes
keygenes <- ordergene[1:8] 
#cds_subset <- cds[keygenes,]

print(ordergene)

# Select genes of interest (user needs to provide gene names or indices of interest)
interested_genes <- c( "PLOD1", "SLC2A5", "TNFRSF14","TNFRSF9", "ERRFI1", "ISG15", "AURKAIP1", "ENO1")

keygenes <- intersect(interested_genes, rownames(cds))  # Ensure genes exist in the dataset

# Check selected genes
if (length(keygenes) == 0) {
  stop("No genes of interest found, please check gene names.")
} else {
  message("Found the following genes of interest: ", paste(keygenes, collapse = ", "))
}

# Subset data
cds_subset <- cds[keygenes,]

# Custom colors
custom_colors <- c(
  "Tumor" = '#FF6666',       # Red
  "Normal" = '#E5D2DD',     # Green
  "BC-cancer" = '#FF6666',  # Blue
  "BC-no_cancer" = '#E5D2DD' # Orange
)

for (type in types) {
  if (type == "Pseudotime") {
    # Use gradient color for continuous data
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # Set gradient color
      theme(legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16),
            legend.position = "top",
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            strip.text = element_text(size = 14)) +
      facet_wrap(~ gene_short_name, ncol = 2)
  } else if (type %in%c("State","celltype")) {
    # Use custom colors for discrete data
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_manual(values = col) +  # Apply discrete colors
      theme(legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16),
            legend.position = "top",
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            strip.text = element_text(size = 14)) +
      facet_wrap(~ gene_short_name, ncol = 2)+
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # Use custom colors for discrete data
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_manual(values = col) +  # Apply discrete colors
      theme(legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16),
            legend.position = "top",
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            strip.text = element_text(size = 14)) +
      facet_wrap(~ gene_short_name, ncol = 2)+
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  
  pdf(file = paste0(output, "/keygene_", type, ".pdf"),width = 8, height = 8)
  print(plot_cell_pseu)
  dev.off()
  
  ggsave(filename = paste0(output, "/keygene_", type, ".svg"), plot = plot_cell_pseu, width = 8, height = 6)
}



