####################### Seurat Analysis #####################
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

# Construct full path
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)

# Generate violin plot to display QC metrics
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 5)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4, group.by = "treatment", pt.size = 0, cols = col)
dev.off()



#### 5. Expression normalization ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Calculate genes with significant expression variation using FindVariableFeatures
ScRNA <- FindVariableFeatures(ScRNA, selection.method = "vst",
                              nfeatures = 2000)

# Display highly variable genes
pdf(paste(OUTPUT,"variable gene.pdf"),width = 9,height = 6)
top10 <- head(VariableFeatures(ScRNA), 10) 
plot1 <- VariableFeaturePlot(ScRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=3)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()


#### 5. Expression normalization ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Calculate genes with significant expression variation using FindVariableFeatures
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
# Scale data
ScRNA <- ScaleData(ScRNA)

# Run PCA
ScRNA <- RunPCA(ScRNA, npcs = 30)

pdf(paste(OUTPUT,"Dimplot.pdf"),width = 9,height = 6)
DimPlot(object = ScRNA, reduction = "pca", pt.size = .1, group.by = "treatment", cols = col)
dev.off()

pdf(paste(OUTPUT,"vlnplot.pdf"),width = 9,height = 6)
VlnPlot(object = ScRNA, features = "PC_1", group.by = "treatment", pt.size = 0, cols = col)
dev.off()

# PCA visualization
pdf(paste(OUTPUT, "DimHeatmap.pdf"),width = 9,height = 6)
DimHeatmap(ScRNA, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# Evaluate PC dimensions
pdf(paste0(OUTPUT,"PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(ScRNA)
dev.off()



save(ScRNA, file = "ScRNA_before_clustering_after_batch_correction.RData")



#### 7. Cell clustering and annotation ####

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#6A4C93',
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA_before_clustering_after_batch_correction.RData")
# Cell clustering
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:20, spread = 0.5) %>% 
  #RunTSNE(dims = 1:20) %>%
  FindNeighbors(dims = 1:20)

ScRNA <- FindClusters(ScRNA, resolution = seq(from = 0.1, 
                                              to = 1.0, 
                                              by = 0.1))
#pdf(paste(OUTPUT, "clustree.pdf"),width=10,height=9)
#library(clustree)
#clustree(ScRNA)
#dev.off()

#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident ## Select the resolution based on the clustering tree
table(Idents(ScRNA))

#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))

# Ensure the factor levels of "treatment" are ordered as Non-infected and Infected
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Display clusters split by Non-infected and Infected order
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 6*length(unique(ScRNA$treatment)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size = 2, label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Display clusters split by sample
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 6*length(unique(ScRNA$orig.ident)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size = 2, label = TRUE, repel = TRUE, split.by = "orig.ident", cols = col)
dev.off()

# Generate UMAP plot only
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size = 2, label = TRUE, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

DimPlot(ScRNA, reduction = "umap", pt.size = 2, label = FALSE, repel = TRUE, cols = col, group.by ="treatment") + ggtitle(NULL) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())
dev.off()

# Display tumor and normal together
pdf(paste(OUTPUT, "cluster-diff_umap.pdf"), width = 6, height = 6)
DimPlot(ScRNA, repel = TRUE, pt.size = 2, 
        reduction = "umap",
        group.by = "treatment") +
  scale_color_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
        legend.position = c(.01, .1)) +
  labs(title = "Sample Origin")
dev.off()

saveRDS(ScRNA, "ScRNA_after_clustering.rds")



### Calculate doublet rate - DoubletFinder ####
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder',force = TRUE)
library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)


output <- paste(outdir, "Doublet", sep = '/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA_after_clustering.rds")
keloid <- readRDS(file_path)


## (2) pK Identification ----------------------------------------------------------
# This is a process for testing the optimal parameter and runs slowly
sweep.res.list_keloid <- paramSweep(keloid, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_keloid)
sweep.stats_keloid <- summarizeSweep(sweep.res.list_keloid, GT = FALSE)
bcmvn_keloid <- find.pK(sweep.stats_keloid) # The point with the best parameter can be identified
## Therefore, the optimal parameter is:
mpK <- as.numeric(as.character(bcmvn_keloid$pK[which.max(bcmvn_keloid$BCmetric)]))



## (3) Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- keloid@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(keloid)*8*1e-6 # Calculated as the doublet rate increases by 0.8% for every additional 1000 cells
DoubletRate = 0.042104
# Estimate the homotypic doublet proportion. Artificial doublets are generated based on the parameter in modelHomotypic(). Here they are mixed from seurat_clusters.


#nExp_poi <- round(0.008 *nrow(keloid@meta.data)) 
nExp_poi <- round(DoubletRate*length(keloid$seurat_clusters))  # It is better to provide celltype rather than seurat_clusters.
# Calculate doublet proportion
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## (4) Finally, use the selected parameters to identify Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
keloid <- doubletFinder(keloid, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, sct = F)
keloid <- doubletFinder(keloid, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, sct = F)  # reuse.pANN = FALSE,
# Perform doublet identification using nExp = nExp_poi and nExp = nExp_poi.adj, respectively, to determine which cells are Doublet-High Confidence later.


## Plot results ---------------------------------------------------------------------------
View(keloid@meta.data)

keloid@meta.data[,"DF_hi.lo"] <- keloid@meta.data$DF.classifications_0.25_0.3_25
keloid@meta.data$DF_hi.lo[which(keloid@meta.data$DF_hi.lo == "Doublet" & keloid@meta.data$DF.classifications_0.25_0.3_19 == "Singlet")] <- "Doublet-Low Confidence"
keloid@meta.data$DF_hi.lo[which(keloid@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidence"
table(keloid@meta.data$DF_hi.lo)
# Doublet-High Confidence  Doublet-Low Confidence                  Singlet 
# 198                      24                                      5041 

## Result display: the classification results are stored in pbmc@meta.data
pdf(paste(output, "doubletfinder.pdf", sep = '/'), width = 4, height = 3)
DimPlot(keloid, reduction = "umap", group.by ="DF.classifications_0.25_0.3_19", cols = c("#DC050C", "#1965B0")) +
  ggtitle("DoubletFinder") 
dev.off()





########## Check the expression ratio of "tdTomato" and "Epcam" ###########
##### Combined plotting ######
genes <- c("EPCAM","KRT7","KRT19")
subset_data <- ScRNA

for (gene in genes) {
  # Get gene expression data
  gene_expr <- FetchData(subset_data, vars = gene)
  subset_data[[paste0(gene, "_expr")]] <- gene_expr[[gene]]
  
  # Calculate the proportion of non-zero expressing cells
  expressed_cells <- sum(subset_data[[paste0(gene, "_expr")]] > 0)
  total_cells <- nrow(subset_data@meta.data)
  expression_ratio <- expressed_cells / total_cells * 100
  
  # Extract expression vector from metadata
  expr_vec <- subset_data@meta.data[[paste0(gene, "_expr")]]
  
  # Sort cells (low expression plotted underneath)
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
      pt.size = 2,
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
      pt.size = 2,
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



########### Manual cell annotation ########

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


setwd("/out/")
outdir <- "/out/"


output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA_after_clustering.rds")
scedata <- readRDS(file_path)

# Define marker genes for different cell types
cellmarker <- c(
  "SERPINB3", "CEACAM5", "ENO2", "KRT19","GRP", 
  "EPCAM","KRT7","KRT19",   # Cancer cells
  "SCGB1A1", "SCGB3A2",       # Club cells
  #"KCNE1", "FOXJ1", "TPPP3", "TUBB4B", "TUBB", "TP73", "CCDC7", # Ciliated cells
  "EPCAM", "KRT18",  "KRT19", "MUC1", "KRT8" ,        # Epithelial Cells / epithelial ovarian cancer
  #"STEAP4","CEACAM6","SCGB1A1","MUC5B","MUC5AC","SPDEF","FOXJ1",    # Secretory epithelial cells
  "DNAH5", "DNAH9", "TEKT1",  # Ciliated epithelial cells
  'AGER', 'CAV1', 'CLIC5' ,'HOPX', 'SEMA3E', 'COL4A3',  # AT1
  "LAMP3", "ABCA3", "SLC34A2", "LPCAT1", "SFTPC", "SFTPA1", "SFTPB", "SFTPD"  # AT2 cells
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Use DotPlot to visualize marker gene expression in immune cells
library(ggplot2)
plot <- DotPlot(scedata, features = unique(cellmarker)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 18),  # Increase X-axis text size
    axis.text.y = element_text(size = 18),  # Increase Y-axis text size
    legend.title  = element_text(size = 18),
    legend.text = element_text(size = 16)   # Increase legend text size
  ) +
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066','#336699','#66CC66','#FFCC33'))

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 10, height = 4)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 10, height = 4)




library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)


col <- c(
  # UMAP
  "#CCCCCC","#3C8487", "#D0F199", "#79BC98",  "#094867",'#E59CC4',"#6666CC",
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue -> green -> light green gradient
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq light colors
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)

file_path <- file.path(outdir, "ScRNA_after_clustering.rds")
scedata <- readRDS(file_path)

# Annotate cell types for each cluster
scedata <- RenameIdents(scedata, c(
  "0"="Unknown", 
  "1"="Unknown", 
  "2"="Lung cancer", 
  "3"="Unknown", 
  "4"="Lung cancer", 
  "5"="Unknown")
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
# Plot cell type UMAP
pdf(paste(output, "ann_umap.pdf",sep = '/'), width = 6, height = 4)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        axis.title.x = element_text(size = 12, face = "bold"),  # Increase X-axis title size
        axis.title.y = element_text(size = 12, face = "bold"),  # Increase Y-axis title size
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

dev.off() 

#        legend.position = c(0.99, 0.12),  # Move legend to the lower-right corner
#        legend.justification = c("right", "bottom"))

# Plot cell type UMAP
svg(paste(output, "ann_umap.svg",sep = '/'), width = 6, height = 4)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        axis.title.x = element_text(size = 12, face = "bold"),  # Increase X-axis title size
        axis.title.y = element_text(size = 12, face = "bold"),  # Increase Y-axis title size
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

#        legend.position = c(0.99, 0.12),  # Move legend to the lower-right corner
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width = 5*length(unique(scedata$treatment)),height = 5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # Increase facet title size
    axis.text.x = element_text(size = 16),  # X-axis label size
    axis.text.y = element_text(size = 16),  # Y-axis label size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
    axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18)    # Increase legend text size
  )
dev.off()

svg(paste(output, "ann-diff-umap.svg",sep = '/'),width = 5*length(unique(scedata$treatment)),height = 5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # Increase facet title size
    axis.text.x = element_text(size = 16),  # X-axis label size
    axis.text.y = element_text(size = 16),  # Y-axis label size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
    axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18)    # Increase legend text size
  )
dev.off()


###### Add count and proportion labels ######
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)

# Calculate count and proportion for each cell type
cell_counts <- table(scedata$celltype)
cell_prop <- prop.table(cell_counts)

cell_stats <- data.frame(
  celltype = names(cell_counts),
  count = as.numeric(cell_counts),
  prop = round(100 * as.numeric(cell_prop), 1)
)

# Extract UMAP coordinates
umap_data <- Embeddings(scedata, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(celltype = scedata$celltype)

# Calculate center coordinates for each group to place labels
centers <- umap_data %>%
  group_by(celltype) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# Merge count and proportion information
centers <- centers %>%
  left_join(cell_stats, by = "celltype") %>%
  mutate(label = paste0("n = ", count, " cells (", prop, "%)"))

# Plot PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 6, height = 4)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 2,
        label = FALSE, cols = col) +
  geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = label),
            color = "black", size = 5, fontface = "bold") +
  theme_dr(xlength = 0.15,
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_blank())
dev.off()

# Plot SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 6, height = 4)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 2,
        label = FALSE, cols = col) +
  geom_text(data = centers, aes(x = UMAP_1, y = UMAP_2, label = label),
            color = "black", size = 5, fontface = "bold") +
  theme_dr(xlength = 0.15,
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_blank())
dev.off()





# Define marker genes for different cell types
cellmarker <- c(
  "EPCAM","MUC1","ERBB2","MET"      # Cancer cells
  
)


# Create lists to store all RidgePlot, VlnPlot, and FeaturePlot objects
ridge_plots <- list()
vln_plots <- list()
feature_plots <- list()

# Iterate over the marker gene list
for (gene in cellmarker) {
  # RidgePlot
  ridge_plots[[gene]] <- RidgePlot(scedata, features = gene, ncol = 1, cols = col) +
    theme(legend.position = "none")  
  
  # VlnPlot
  vln_plots[[gene]] <- VlnPlot(scedata, features = gene, ncol = 1, pt.size = 0, cols = col) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") 
  
  # Get gene expression data and ensure it is returned as a numeric vector
  gene_expr <- as.numeric(FetchData(scedata, vars = gene, assay = "RNA")[[gene]])  # Convert to numeric vector
  
  # Check whether gene expression was successfully extracted
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # Add gene expression data to metadata
  scedata@meta.data[[paste0(gene, "_expr")]] <- gene_expr  # Add expression information to metadata
  
  # Sort using numeric vector to avoid directly indexing the ScRNA object
  cells_ordered <- scedata@meta.data[order(scedata@meta.data[[paste0(gene, "_expr")]], 
                                           decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # Extract sorted cell names
  
  # Plot FeaturePlot according to sorted cell order and use a continuous color gradient
  feature_plots[[gene]] <- FeaturePlot(scedata, features = gene, reduction = "umap", 
                                       cells = cell_names_ordered,  # Specify cell order
                                       ncol = 1) + 
    scale_color_gradientn(colors = c(  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
                                       "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
                                       "#F5E24B")) +  # Set continuous color gradient
    theme(legend.position = "right", 
          plot.title = element_text(size = 22, face = "bold"),  # Increase title size and make bold
          legend.text = element_text(size = 18)) +  # Increase legend text size
    NoAxes()  # Remove axes
  
}


# Save RidgePlot figures
#pdf(paste0(out, "cellmarker_RidgePlot.pdf"), width = 25, height = 12)
#print(cowplot::plot_grid(plotlist = ridge_plots, ncol = 4))
#dev.off()

# Save FeaturePlot figures
pdf(paste0(output, "/marker_FeaturePlot_umap.pdf"), width = 10, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 2))
dev.off()

svg(paste0(output, "/marker_FeaturePlot_umap.svg"), width = 10, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 2))
dev.off()

pdf(paste0(output, "/marker_VlnPlot_umap.pdf"), width = 10, height = 3)
print(cowplot::plot_grid(plotlist = vln_plots, ncol = 3))
dev.off()
svg(paste0(output, "/marker_VlnPlot_umap.svg"), width = 10, height = 3)
print(cowplot::plot_grid(plotlist = vln_plots, ncol = 3))
dev.off()






# Load required packages
if (!require("ggridges")) install.packages("ggridges")
if (!require("cowplot")) install.packages("cowplot")
if (!require("patchwork")) install.packages("patchwork")
if (!require("grid")) install.packages("grid")
library(ggridges)
library(cowplot)
library(patchwork)
library(grid)

# Define marker genes for different cell types
cellmarker <- c("EPCAM", "MUC1", "ERBB2", "MET")  # Cancer cells

# Create a list to store all combined plots
combined_plots <- list()

# Iterate over the marker gene list
for (gene in cellmarker) {
  # Get gene expression data
  gene_expr <- as.numeric(FetchData(scedata, vars = gene, assay = "RNA")[[gene]])
  
  # Check whether gene expression was successfully extracted
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # Get UMAP coordinates
  umap_data <- as.data.frame(scedata@reductions$umap@cell.embeddings)
  colnames(umap_data) <- c("UMAP_1", "UMAP_2")
  umap_data$expression <- gene_expr
  
  # Create UMAP scatter plot - add border and axis titles
  p_main <- FeaturePlot(scedata, features = gene, reduction = "umap", 
                        pt.size = 0.5, order = TRUE) + 
    scale_color_gradientn(colors = c("#ADD8C0", "#094867"),
                          name = "Expression") +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      #legend.key.height = unit(1.5, "cm"),
      #legend.key.width = unit(0.5, "cm"),
      # Add panel border
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      # Axis settings
      axis.line = element_line(color = "black", size = 0.8),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.text = element_text(size = 12, color = "black")
    ) +
    labs(
      title = gene,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    # Add axis ticks
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05))
  
  # Create expression level categories
  umap_data$expr_level <- "Low"
  if (any(umap_data$expression > 0)) {
    expr_positive <- umap_data$expression[umap_data$expression > 0]
    median_expr <- quantile(expr_positive, 0.5)
    
    umap_data$expr_level[umap_data$expression > 0 & 
                           umap_data$expression <= median_expr] <- "Medium"
    umap_data$expr_level[umap_data$expression > median_expr] <- "High"
  }
  
  # Set factor order to ensure correct legend order
  umap_data$expr_level <- factor(umap_data$expr_level, 
                                 levels = c("Low", "Medium", "High"))
  
  # Method 1: Use density-shaded area plot to show expression distribution along the UMAP1 axis (with legend)
  p_top <- ggplot(umap_data, aes(x = UMAP_1, fill = expr_level)) +
    # Add density shading for cells with different expression levels
    geom_density(data = subset(umap_data, expr_level == "Low"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.5, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "Medium"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "High"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    # Add smoothed expression trend line
    geom_smooth(aes(y = expression, color = "Expression Trend"), 
                method = "loess", span = 0.3, 
                fill = "#20958B", alpha = 0.2, size = 0.8) +
    scale_fill_manual(
      name = "Expression Level",
      values = c("Low" = "#E0E0E0", "Medium" = "#6ECC5A", "High" = "#F5E24B"),
      labels = c("Low (expr = 0)", "Medium", "High"),
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        nrow = 1,
        keywidth = unit(0.8, "cm"),
        keyheight = unit(0.4, "cm")
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Expression Trend" = "#20958B")
    ) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(2, 0, 2, 0),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL)
  
  # Use density-shaded area plot to show expression distribution along the UMAP2 axis (with legend)
  p_right <- ggplot(umap_data, aes(x = UMAP_2, fill = expr_level)) +
    # Add density shading for cells with different expression levels
    geom_density(data = subset(umap_data, expr_level == "Low"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.5, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "Medium"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "High"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    # Add smoothed expression trend line
    geom_smooth(aes(y = expression, color = "Expression Trend"), 
                method = "loess", span = 0.3, 
                fill = "#20958B", alpha = 0.2, size = 0.8) +
    scale_fill_manual(
      name = "Expression Level",
      values = c("Low" = "#E0E0E0", "Medium" = "#6ECC5A", "High" = "#F5E24B"),
      labels = c("Low (expr = 0)", "Medium", "High"),
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        nrow = 1,
        keywidth = unit(0.8, "cm"),
        keyheight = unit(0.4, "cm")
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Expression Trend" = "#20958B")
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    coord_flip()
  
  # Create an empty plot as a placeholder
  p_empty <- ggplot() + 
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  # Combine plots using patchwork
  combined_plots[[gene]] <- p_top + p_empty + p_main + p_right +
    plot_layout(ncol = 2, nrow = 2, 
                widths = c(4, 1), 
                heights = c(1, 4))
}



# Provide a more compact version by combining all genes into one figure
if (length(cellmarker) > 0) {
  # Create a multi-page PDF, one gene per page
  pdf(paste0(output, "/marker_FeaturePlot_individual.pdf"), width = 5, height = 6)
  for (gene in cellmarker) {
    print(combined_plots[[gene]])
  }
  dev.off()
  
  # Create a combined figure (all genes on one page)
  combined_grid <- cowplot::plot_grid(plotlist = combined_plots, ncol = 2, nrow = ceiling(length(cellmarker)/2))
  
  pdf(paste0(output, "/marker_FeaturePlot_combined.pdf"), width = 10, height = 6 * ceiling(length(cellmarker)/2))
  print(combined_grid)
  dev.off()
  
}










####### Calculate cell proportions ###########

col <- c(
  # UMAP
  "#3C8487", "#D0F199", "#79BC98",  "#094867",'#E59CC4',"#6666CC",
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue -> green -> light green gradient
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

######## Calculate cell counts of different cell groups across all samples
cell_counts <- as.data.frame(table(Idents(scedata)))
colnames(cell_counts) <- c("CellType", "Counts")

# Sort by cell count from high to low
cell_counts <- cell_counts[order(-cell_counts$Counts), ]
# Save cell counts of all cell groups to the specified directory
write.csv(cell_counts, paste(output, "cell_counts.csv", sep='/'), row.names = FALSE)

# Select the top 11 cell groups and save to the specified directory
cell_counts_top9 <- head(cell_counts, 11)
write.csv(cell_counts_top9, paste(output, "cell_counts_top9.csv", sep='/'), row.names = FALSE)

# Load required package
library(ggplot2)
# Plot bar chart
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

# Save figures in PDF format
ggsave(paste(output, "cell_type_distribution.pdf", sep='/'), plot = p, width = 7, height = 6, dpi = 800)
ggsave(paste(output, "cell_type_distribution.svg", sep='/'), plot = p, width = 7, height = 6, dpi = 800)


# Calculate the number of different cell groups in each group
# Calculate the count of each cell type grouped by sample
cell_counts_group <- as.data.frame(table(scedata$orig.ident, Idents(scedata)))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# Add group information (assuming the grouping variable is `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # Ensure the uniqueness of group information
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# Calculate the proportion of each cell type in each sample
cell_counts_group <- cell_counts_group %>%
  group_by(Sample) %>%
  mutate(Ratio = Counts / sum(Counts))

p <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # Modify X-axis labels
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size
# Add cell count text labels
p <- p + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)


p <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # Modify X-axis labels
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size
# Add cell proportion text labels
p <- p + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)



############ Grouping ############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate the proportion of each cell type in each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

########## Plot stacked bar chart of cell counts ##########
p1 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 24, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 24),  # Modify Y-axis text size
        axis.title.y = element_text(size = 26), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 24))  # Modify legend text size

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)

########## Plot stacked bar chart of cell proportions ##########
p2 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 24, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 24),  # Modify Y-axis text size
        axis.title.y = element_text(size = 26), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 24))  # Modify legend text size

file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)


