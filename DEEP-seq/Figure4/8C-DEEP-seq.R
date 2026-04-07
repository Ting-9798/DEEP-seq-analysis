# Clear environment
rm(list = ls())

# Load required packages
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(grid)
library(tidydr)
#install.packages("tidydr")


setwd("/data/")
# Define folder path
data_dir <- "/data/"

folders <- c("A","B")

# Initialize Seurat object list
seurat.list <- list()

# Loop through each subfolder to read expression matrix files and perform Seurat analysis
for (folder in folders) {
  folder_path <- file.path(data_dir, folder)
  
  # Read 10X data format
  sample_data <- Read10X(data.dir = folder_path)
  sample_name <- basename(folder)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 10, min.features = 1000)
  
  # Add mitochondrial gene percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Calculate ribosomal protein percentage
  seurat_obj[["percent.rps"]] <- PercentageFeatureSet(seurat_obj, pattern = c("^RPS","^RPL"))
  
  # Quality control filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 25)
  
  # Add treatment column
  seurat_obj$treatment <- sample_name 
  
  
  # Add processed Seurat object to list
  seurat.list[[sample_name]] <- seurat_obj
}


combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)
combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
View(combined_seurat@meta.data)

# Save merged Seurat object
saveRDS(combined_seurat, file = "/out/combined_seurat.rds")


####################### Seurat analysis #####################
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

# Concatenate full path
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
View(ScRNA@meta.data)

# Generate violin plots to display QC metrics
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 6)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"), ncol = 4, group.by = "treatment", pt.size = 0,cols = col)
dev.off()


#### 5. Expression normalization ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Calculate highly variable genes
ScRNA <- FindVariableFeatures(ScRNA, selection.method = "vst",
                              nfeatures = 2000)

# Display variable genes
pdf(paste(OUTPUT,"variable gene.pdf"),width = 9,height = 6)
top10 <- head(VariableFeatures(ScRNA), 10) 
plot1 <- VariableFeaturePlot(ScRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=3)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()

#### 6. Normalization and PCA dimensionality reduction ####
# Scaling
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

#pdf(paste(OUTPUT,  "Dimplot-corret.pdf"),width = 12,height = 6)
#DimPlot(object = ScRNA, reduction = "harmony",
#        pt.size = 0.1, group.by = "treatment")
#dev.off()

#pdf(paste(OUTPUT, "vlnplot-corret.pdf"),width = 12,height = 6)
#VlnPlot(object = ScRNA, features = "harmony_1", 
#        group.by = "treatment", pt.size =0)
#dev.off()


save(ScRNA, file = "ScRNA_before_clustering_after_batch_correction.RData")



#### 7. Cell clustering and annotation ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#FEDD81", "#FF9A84",  "#094867",'#E59CC4',"#6666CC",
  "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue -> green -> light green gradient
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq light colors
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


load("ScRNA_before_clustering_after_batch_correction.RData")


ScRNA <- ScRNA %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)

ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))

#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident## Select the desired resolution based on the clustering tree
table(Idents(ScRNA))

# Ensure treatment factor levels are ordered as Non-infected and Infected
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Display clusters split by Non-infected and Infected order
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, label.size = 5,repel = TRUE, split.by = "treatment", cols = col)+
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

# Display clusters split by sample
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 40, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, label.size = 5,repel = TRUE, split.by = "orig.ident", cols = col)+
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

# Generate UMAP plot only
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE,label.size = 5, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size=14),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20))

DimPlot(ScRNA, reduction = "umap", label = FALSE,label.size = 5, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size=14),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20))
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



# Generate UMAP plot only
pdf(paste(OUTPUT, "cluster_umap_11.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
dev.off()



saveRDS(ScRNA, "ScRNA_after_clustering.rds")


# Select the first 6 clusters to display marker gene expression trends
output <- paste(outdir,'Cell_localization', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA_after_clustering.rds")
ScRNA <- readRDS(file_path)

# Identify marker genes
ScRNA.markers <- FindAllMarkers(ScRNA, only.pos = TRUE,   ### only.pos = TRUE: only identify upregulated genes
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
        strip.text.x = element_text(size = 30))  # Modify group label color and size
dev.off()

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 100,height = 10)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): rotate X-axis text
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18))  # Increase Y-axis title text size
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=100*dpi,h=10*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): rotate X-axis text
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 20),  # Increase X-axis title text size
        axis.title.y = element_text(size = 20))  # Increase Y-axis title text size
dev.off()




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


# Display nFeature_RNA on UMAP and save as PDF
pdf(file = paste(output, "nFeature_RNA_UMAP.pdf", sep='/'),width = 7, height = 6)
FeaturePlot(scedata, features = "nFeature_RNA", reduction = "umap", cols = c("lightgrey", "#FF3366"))
dev.off()

# Display nCount_RNA on UMAP and save as PDF
pdf(file = paste(output, "nCount_RNA_UMAP.pdf", sep='/'),width = 7, height = 6)
FeaturePlot(scedata, features = "nCount_RNA", reduction = "umap", cols = c("lightgrey", '#6A4C93'))
dev.off()



# Get TPRX1 gene expression data
TPRX1_expression <- scedata[["RNA"]]@data["TPRX1", ]

# Create a new column named 'celltype' based on TPRX1 expression
scedata$celltype <- ifelse(TPRX1_expression > 0, "8CLC", "non-8CLC")
#View(scedata@meta.data)

# Plot cell type UMAP
pdf(paste(output, "ann_umap_8c.pdf", sep='/'), width = 6, height = 5)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=1,label=TRUE,label.size = 5,repel = TRUE,cols=col,
        label.box = TRUE)+
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

saveRDS(scedata,  "celltype_8C.rds")

###### Check TPRX1 expression ########
# Get TPRX1 gene expression data
tdTomato_expr <- FetchData(scedata, vars = "TPRX1")

# Add expression values to cell metadata
scedata$tdTomato_expr <- tdTomato_expr$TPRX1

# Calculate TPRX1 expression ratio (percentage of non-zero expressing cells)
expressed_cells <- sum(scedata$tdTomato_expr > 0)
total_cells <- nrow(scedata@meta.data)
expression_ratio <- expressed_cells / total_cells * 100

# Sort cells by expression level
cells_ordered <- scedata@meta.data[order(scedata$tdTomato_expr, decreasing = FALSE), ]

# Extract sorted cell names
cell_names_ordered <- rownames(cells_ordered)

# Set title with expression ratio annotation
plot_title <- paste0("TPRX1 Expression (", round(expression_ratio, 2), "%)")

# Plot FeaturePlot using sorted cell order
pdf(paste0(output, "/TPRX1_FeaturePlot_umap.pdf"), width = 4, height = 4)
FeaturePlot(
  scedata, 
  features = "TPRX1",
  reduction = "umap", 
  cells = cell_names_ordered,  # Plot according to cell order
  ncol = 1,
  cols = c('#E5D2DD',  "#FF3366")
) +
  ggtitle(plot_title) +  # Add title
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + 
  NoAxes()  # Remove axes
dev.off()

svg(paste0(output, "/TPRX1_FeaturePlot_umap.svg"), width = 4, height = 4)
FeaturePlot(
  scedata, 
  features = "TPRX1",
  reduction = "umap", 
  cells = cell_names_ordered,  # Plot according to cell order
  ncol = 1,
  cols = c('#E5D2DD',  "#FF3366")
) +
  ggtitle(plot_title) +  # Add title
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + 
  NoAxes()  # Remove axes
dev.off()



####### Calculate cell proportions ###########
col <- c('#E5D2DD','#FF6666',"#66CCCC","#A4CDE1","#CCFFCC","#FF3366",'#58A4C3',"#FFFFCC",'#E5D2DD',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype_8C.rds")
scedata <- readRDS(file_path)
View(scedata@meta.data)

table(scedata$celltype)

# Calculate cell counts for different cell groups in each sample
# Group by sample and calculate the count of each cell type
cell_counts_group <- as.data.frame(table(scedata$orig.ident,scedata$celltype))
colnames(cell_counts_group) <- c("Sample", "CellType","Counts")

cell_counts_group$CellType <- factor(cell_counts_group$CellType, levels = c( "non-8CLC","8CLC"))

# Sort by CellType
cell_counts_group <- cell_counts_group %>%
  arrange(CellType)

# Add group information (assuming grouping variable is `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # Ensure uniqueness of group information
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# Calculate the proportion of each cell type in each sample
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

file_path <- paste0(output, "/genecount_8c.pdf")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount_8c.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)


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

file_path <- paste0(output, "/geneRatio_8c.pdf")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio_8c.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)





cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)


cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Use DotPlot to visualize immune cell marker gene expression
library(ggplot2)
plot <- DotPlot(scedata, features = unique(cellmarker))+
  theme_bw()+theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 16),  # Increase X-axis text size
    axis.text.y = element_text(size = 16),  # Increase Y-axis text size
    legend.title  = element_text(size = 18),
    legend.text = element_text(size = 16)   # Increase legend text size
  ) +
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 14, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 14, height = 5)




######### Plot cell annotation heatmap ############
# change annotation color
library("scales")
library(ggsci)
library(scRNAtoolVis)

# Set annotation colors
mycol1 <- pal_simpsons()(18)


#file_path <- file.path(outdir, "celltype.rds")
#scedata <- readRDS(file_path)

pdf(file = paste(output, "ann_Heatmap_11.pdf",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # Custom high-value color
dev.off()

svg(file = paste(output, "ann_Heatmap_11.svg",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # Custom high-value color
dev.off()



######### Plot marker gene expression violin plot

library(ggplot2)
library(reshape2)
library(dplyr)

cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)

cellmarker <- c(
  "ZSCAN5B","ZNF280A","TPRX1","DUXA","DNMT3L","KLF17", "DPPA3","KHDC1L","KLF4","SOX2", "NANOG","CD24"
)

# Filter marker genes present in the dataset
existing_markers <- cellmarker[cellmarker %in% rownames(scedata[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# Extract expression data and convert to data format suitable for plotting
vln.df <- as.data.frame(scedata[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

# Continue with the original steps
anno <- scedata@meta.data[, c("CB", "seurat_clusters")]
vln.df <- inner_join(vln.df, anno, by = "CB")
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)


# Plot violin plot with swapped X and Y axes
plot <- vln.df %>%
  ggplot(aes(exp, seurat_clusters)) +
  geom_violin(aes(fill = seurat_clusters), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # Adjust facet_grid with genes as columns
  scale_fill_manual(values = col) +
  scale_x_continuous("") + scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 28),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    strip.text = element_text(size = 18, face = "bold"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Save violin plot
ggsave(filename = paste(output, "marker_ViolinPlot.pdf", sep = '/'), plot = plot, width = 18,height=5,limitsize = FALSE)
ggsave(filename = paste(output, "marker_ViolinPlot.svg", sep = '/'), plot = plot, width = 18,height=5,limitsize = FALSE)



library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)

col <- c('#E5D2DD','#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FF3366",'#58A4C3',"#FFFFCC",'#E5D2DD',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

file_path <- file.path(outdir, "ScRNA_after_clustering.rds")
scedata <- readRDS(file_path)

# Annotate cell types for each cluster
scedata <- RenameIdents(scedata, c(
  "0"="non-8CLC",
  "1"="non-8CLC", 
  "2"="non-8CLC",
  "3"= "non-8CLC",
  "4"="8CLC",
  "5"="non-8CLC",
  "6"="non-8CLC",
  "7"="non-8CLC",
  "8"= "non-8CLC",
  "9"="non-8CLC",
  "10"="non-8CLC",
  "11"="non-8CLC")
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
pdf(paste(output, "ann_umap.pdf",sep = '/'), width = 6, height = 5)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=1,label=TRUE,label.size = 6,repel = TRUE,cols=col,
        label.box = TRUE)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size = 16),
        legend.position = "none",
        #legend.title = element_text(size = 18), 
        #legend.text = element_text(size = 18),
        plot.title = element_blank())

dev.off() 

#        legend.position = c(0.99, 0.12),  # Move legend to the lower right corner
#        legend.justification = c("right", "bottom"))


# Plot cell type UMAP
svg(paste(output, "ann_umap.svg",sep = '/'), width = 6, height = 5)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=1,label=TRUE,label.size = 5,repel = TRUE,cols=col,
        label.box = TRUE)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size = 16),
        legend.position = "none",
        #legend.title = element_text(size = 18), 
        #legend.text = element_text(size = 18),
        plot.title = element_blank())
#        legend.position = c(0.99, 0.12),  # Move legend to the lower right corner
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=6*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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

svg(paste(output, "ann-diff-umap.svg",sep = '/'),width=6*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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


########### Plot annotated cluster dot plot
library(ggh4x)

cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)

cellmarker <- c(
  "ZSCAN5B","ZNF280A","TPRX1","DUXA","DNMT3L","KLF17", "DPPA3","KHDC1L","KLF4","SOX2", "NANOG","CD24"
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Generate DotPlot
p <- DotPlot(scedata, features = unique(cellmarker))

# Extract data
dat <- p$data

# Get cluster annotations and merge
anno <- distinct(data.frame(id = scedata$celltype, celltype = scedata$seurat_clusters))
colnames(anno) <- c("id", "celltype")
df <- left_join(dat, anno, by = "id")

# Define cluster order
cluster.order <- c(0,2,1,3,4,5,6,7,8,9,10,11,12)
df$celltype <- factor(df$celltype, levels = cluster.order)

# Plot dot plot
p <- ggplot(df, aes(features.plot, interaction(celltype, id), size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank(),
    axis.text.y = element_text(color = 'black', size = 18),
    axis.text.x = element_text(color = 'black', size = 16, angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_fill_gradientn(colours = c('#5749a0',  '#00bbb1', '#bef0b0', '#fdf4af', '#f9b64b', '#ec840e', '#ca443d', '#FF6666')) +
  guides(y = "axis_nested") +
  theme(
    ggh4x.axis.nesttext.y = element_text(colour = c('#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0', '#24796C', '#DAA51B', '#2F8AC4', '#764E9F', '#ED645A', '#CC3A8E')),
    ggh4x.axis.nestline.y = element_line(size = 3)
  )

# Save figure
ggsave(plot = p, filename = paste0(output, "/ann_DotPlot.pdf"), width = 10, height = 5)
ggsave(plot = p, filename = paste0(output, "/ann_DotPlot.svg"), width = 10, height = 5)



######### Plot cell annotation heatmap ############
# change annotation color
library("scales")
library(ggsci)
library(scRNAtoolVis)

# Set annotation colors
mycol1 <- pal_simpsons()(18)


#file_path <- file.path(outdir, "celltype.rds")
#scedata <- readRDS(file_path)

pdf(file = paste(output, "ann_Heatmap.pdf",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # Custom high-value color
dev.off()

svg(file = paste(output, "ann_Heatmap.svg",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # Custom high-value color
dev.off()




######### Plot marker gene expression violin plot

library(ggplot2)
library(reshape2)
library(dplyr)


cellmarker <- c(
  "OCT4","SOX2", "NANOG","DPPA3","KLF17","TFAP2C","ZSCAN4","DUXA","DUXB","TPRX1","ZNF280A",
  "MBD3L2","FAM151A","TRIM43","GATA6"
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Filter marker genes present in the dataset
existing_markers <- cellmarker[cellmarker %in% rownames(ScRNA[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# Extract expression data and convert to data format suitable for plotting
vln.df <- as.data.frame(ScRNA[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

# Continue with the original steps
anno <- ScRNA@meta.data[, c("CB", "celltype")]
vln.df <- inner_join(vln.df, anno, by = "CB")
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)


# Plot violin plot with swapped X and Y axes
plot <- vln.df %>%
  ggplot(aes(exp, celltype)) +
  geom_violin(aes(fill = celltype), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # Adjust facet_grid with genes as columns
  scale_fill_manual(values = col) +
  scale_x_continuous("") + scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 28),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    strip.text = element_text(size = 18, face = "bold"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Save violin plot
ggsave(filename = paste(output, "marker_ViolinPlot_ann.pdf", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)
ggsave(filename = paste(output, "marker_ViolinPlot_ann.svg", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)



####### Calculate cell proportions ###########
col <- c('#E5D2DD','#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FF3366",'#58A4C3',"#FFFFCC",'#E5D2DD',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

table(scedata$seurat_clusters)

######## Calculate cell counts for all clusters across all samples
cell_counts <- as.data.frame(table(Idents(scedata)))
colnames(cell_counts) <- c("CellType", "Counts")

# Sort by cell count in descending order
cell_counts <- cell_counts[order(-cell_counts$Counts), ]
# Save cell counts of all clusters to the specified directory
write.csv(cell_counts, paste(output, "cell_counts.csv", sep='/'), row.names = FALSE)

# Select the top 11 clusters and save to the specified directory
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

# Save images in PDF and SVG format
ggsave(paste(output, "cell_type_distribution.pdf", sep='/'), plot = p, width = 7, height = 6, dpi = 800)
ggsave(paste(output, "cell_type_distribution.svg", sep='/'), plot = p, width = 7, height = 6, dpi = 800)


# Calculate cell counts for different cell groups in each sample
# Group by sample and calculate the count of each cell type
cell_counts_group <- as.data.frame(table(scedata$orig.ident, Idents(scedata)))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# Add group information (assuming grouping variable is `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # Ensure uniqueness of group information
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# Calculate the proportion of each cell type in each sample
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
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)


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
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)



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
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)

########## Plot stacked bar chart of cell proportions ##########
p2 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size

file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)




############# 8C #################

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
col <- c('#E5D2DD','#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#4F6272',"#A4CDE1",'#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


output <- paste(outdir,'8C', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data <- readRDS(file_path)

# Select epithelial cells belonging to the "8CLC" group
Cells.sub <- subset(data@meta.data, celltype == c("8CLC"))
summary(Cells.sub$celltype)
scedata <- subset(data, cells=row.names(Cells.sub))


# Define exp column based on Cbr2 expression level
Cbr2_threshold <- 0  # Assume expression threshold is 0; can be adjusted as needed
scedata@meta.data$exp <- ifelse(scedata@assays$RNA@data["TPRX1", ] > Cbr2_threshold, "TPRX1+", "TPRX1-")

# Check definition results
#View(scedata@meta.data)

library(ggsci)
# Plot cell type UMAP
pdf(file = paste(output, "ann_umap.pdf",sep='/'), width = 5, height = 4)
DimPlot(object=scedata,group.by = "exp",reduction='umap',pt.size=1,label=FALSE,label.size = 6,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.20, 
           ylength = 0.20,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size = 16),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        plot.title = element_blank())
dev.off() 

# Plot cell type UMAP (SVG format)
svg(file = paste(output, "ann_umap.svg",sep='/'), width = 5, height = 4)
DimPlot(object=scedata,group.by = "exp",reduction='umap',pt.size=1,label=FALSE,label.size = 6,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.20, 
           ylength = 0.20,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size = 16),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        plot.title = element_blank())
dev.off() 

# Save updated data object
saveRDS(scedata,"celltype(8C).rds")




####### Calculate cell proportions ###########

output <- paste(outdir,'8C', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype(8C).rds")
scedata <- readRDS(file_path)

table(scedata$exp)

# Calculate cell counts for different cell groups in each sample
# Group by sample and calculate the count of each cell type
cell_counts_group <- as.data.frame(table(scedata$orig.ident, scedata$exp))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# Add group information (assuming grouping variable is `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # Ensure uniqueness of group information
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# Calculate the proportion of each cell type in each sample
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
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


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
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)




############# Compare WF and SNJ 8C ############
setwd("/out/8C/WF_SNJ/")
outdir <- "/out/8C/WF_SNJ/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

file_path <- file.path("/out/celltype(8C).rds")
WF <- readRDS(file_path)

file_path <- file.path("/out(SNJ1)/celltype(8C).rds")
SNJ <- readRDS(file_path)


# Add identifiers for each sample
WF$treatment <- "WF"
SNJ$treatment <- "SNJ"

# Merge the two Seurat objects
merged_seurat <- merge(WF, SNJ, add.cell.ids = c("WF", "SNJ"))

# Optionally add sample information to metadata
merged_seurat$orig.ident <- paste0(merged_seurat$treatment, "_", merged_seurat$orig.ident)
View(merged_seurat@meta.data)

# Save merged object
saveRDS(merged_seurat, file = file.path("merged_WF_SNJ_8C.rds"))



file_path <- file.path("/out/8C/WF_SNJ/merged_WF_SNJ_8C.rds")
ScRNA <- readRDS(file_path)


#### 5. Expression normalization ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Calculate highly variable genes
ScRNA <- FindVariableFeatures(ScRNA, selection.method = "vst",
                              nfeatures = 2000)

# Display variable genes
pdf(paste(OUTPUT,"variable gene.pdf"),width = 9,height = 6)
top10 <- head(VariableFeatures(ScRNA), 10) 
plot1 <- VariableFeaturePlot(ScRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=3)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()


#### 6. Normalization and PCA dimensionality reduction ####
# Scaling
all.genes <- rownames(ScRNA)
ScRNA<-ScaleData(ScRNA,features = all.genes)

# Run PCA
ScRNA<-RunPCA(ScRNA,npcs = 60)

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
ScRNA<-RunHarmony(ScRNA,group.by.vars = c("treatment"),npcs = 60, 
                  plot_convergence = TRUE, verbose = TRUE)

pdf(paste(OUTPUT,  "Dimplot-corret.pdf"),width = 12,height = 6)
DimPlot(object = ScRNA, reduction = "harmony",
        pt.size = 0.1, group.by = "treatment")
dev.off()

pdf(paste(OUTPUT, "vlnplot-corret.pdf"),width = 12,height = 6)
VlnPlot(object = ScRNA, features = "harmony_1", 
        group.by = "treatment", pt.size =0)
dev.off()


save(ScRNA, file = "ScRNA_before_clustering_after_batch_correction.RData")



#### 7. Cell clustering and annotation ####

col <- c('#58A4C3','#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")

load("ScRNA_before_clustering_after_batch_correction.RData")
# Cell clustering
#ScRNA <- ScRNA %>% 
#  RunUMAP(dims = 1:20) %>% 
#  RunTSNE(dims = 1:20) %>%
#  FindNeighbors(dims = 1:20)


# Cell clustering
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:60,reduction = "harmony",
          n.neighbors = 15,     # Default 30; reducing it separates clusters more
          min.dist = 0.1,       # Default 0.3; reducing it makes clusters tighter and more separated
          spread = 2,
          seed.use = 42) %>%       # Increasing spread enlarges distances between clusters
  RunTSNE(dims = 1:60, reduction = "harmony",         # Use more principal components
          perplexity = 50,      # Default is 30; can be increased appropriately
          theta = 0.3,          # More precise calculation
          seed.use = 42) %>%
  FindNeighbors(dims = 1:50)


ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))

#library(clustree)
#pdf(paste(OUTPUT, "clustree.pdf"),width=10,height=9)
#clustree(ScRNA)
#dev.off()


Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident## Select the desired resolution based on clustering results
table(Idents(ScRNA))

# Display clusters split by Non-infected and Infected order
pdf(paste(OUTPUT, "split.by_cluster_tsne.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "tsne", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Generate UMAP plot only
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())
dev.off()

pdf(paste(OUTPUT, "cluster_umap_1.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_blank())
dev.off()

svg(paste(OUTPUT, "cluster_umap_1.svg"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_blank())
dev.off()

pdf(paste(OUTPUT, "cluster-diff_umap.pdf"),width=6,height=6)
DimPlot(ScRNA, repel = TRUE,
        reduction = "umap",
        group.by ="treatment")+
  scale_color_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")
dev.off()



saveRDS(ScRNA, "ScRNA_after_clustering.rds")




########## Check the expression ratio of "TPRX1" and "Epcam" ###########
##### Combined plotting ######
genes <- c("TPRX1")
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
  
  # Sort cells (low expression at bottom layer)
  cells_ordered <- subset_data@meta.data[order(expr_vec, decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)
  
  # Set title
  plot_title <- paste0(gene, " (", round(expression_ratio, 2), "%) ")
  
  # Output paths
  pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge(filter)", ".pdf")
  svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge(filter)", ".svg")
  
  # PDF plot
  pdf(pdf_path, width = 4, height = 4)
  print(
    FeaturePlot(
      subset_data,
      features = gene,
      reduction = "umap",
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



col <- c("#CCCCCC",'#58A4C3','#FF6666','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")

output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA_after_clustering.rds")
scedata <- readRDS(file_path)


# Define exp column based on Cbr2 expression level
Cbr2_threshold <- 0  # Assume expression threshold is 0; can be adjusted as needed
scedata@meta.data$exp <- ifelse(scedata@assays$RNA@data["TPRX1", ] > Cbr2_threshold, "TPRX1+", "TPRX1-")

# Check definition results
#View(scedata@meta.data)

library(ggsci)
# Plot cell type UMAP
pdf(file = paste(output, "ann_umap.pdf",sep='/'), width = 6.5, height = 4)
DimPlot(object=scedata,group.by = "exp",reduction='umap',pt.size=1,label=FALSE,label.size = 6,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.20, 
           ylength = 0.20,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size = 16),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        plot.title = element_blank())
dev.off() 

# Plot cell type UMAP (SVG format)
svg(file = paste(output, "ann_umap.svg",sep='/'), width = 6.5, height = 4)
DimPlot(object=scedata,group.by = "exp",reduction='umap',pt.size=1,label=FALSE,label.size = 6,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.20, 
           ylength = 0.20,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size = 16),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        plot.title = element_blank())
dev.off() 



##### Add total cell counts for each group #####
# Count the number of cells in each treatment
cell_counts <- scedata@meta.data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(n = n()) %>%
  mutate(label = paste0(treatment, " (", n, " cells)"))

# Construct named vector for replacing facet labels
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# Plot PDF
pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),
    width = 4*length(unique(scedata$treatment)), height = 4)

DimPlot(scedata, reduction = "umap",group.by = "exp",  split.by = "treatment",
        pt.size = 1, label = FALSE, cols = col) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # Bold facet titles
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
dev.off()

# Same for SVG output
svg(paste(output, "ann-diff-umap.svg",sep = '/'),
    width = 4*length(unique(scedata$treatment)), height = 4)

DimPlot(scedata, reduction = "umap", group.by = "exp", split.by = "treatment",
        pt.size = 1, label = FALSE, cols = col) +
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


# Save updated data object
saveRDS(scedata,"celltype(8C).rds")




####### Calculate cell proportions ###########

output <- paste(outdir,'celltype/', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype(8C).rds")
scedata <- readRDS(file_path)

table(scedata$exp)

# Calculate cell counts for different cell groups in each sample
# Group by sample and calculate the count of each cell type
cell_counts_group <- as.data.frame(table(scedata$orig.ident, scedata$exp))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# Add group information (assuming grouping variable is `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # Ensure uniqueness of group information
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# Calculate the proportion of each cell type in each sample
cell_counts_group <- cell_counts_group %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(Ratio = Counts / sum(Counts))

p1 <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) + 
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
p1 <- p1 + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


p2 <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) + 
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
p2 <- p2 + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 1.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 1.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)





############ Grouping ############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, scedata$exp))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate the proportion of each cell type in each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(Ratio = Counts / sum(Counts))

########## Plot stacked bar chart of cell counts ##########
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
# Add cell count text labels
p1 <- p1 + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)

########## Plot stacked bar chart of cell proportions ##########
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
# Add cell count text labels
p2 <- p2 + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)


file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)





######### Plot marker gene expression violin plot

col <- c('#FF6666','#58A4C3','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")


library(ggplot2)
library(reshape2)
library(dplyr)


cellmarker <- c(
  "OCT4","SOX2", "NANOG","DPPA3","KLF17","TFAP2C","ZSCAN4","DUXA","DUXB","TPRX1","ZNF280A",
  "MBD3L2","FAM151A","TRIM43","GATA6"
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Filter marker genes present in the dataset
existing_markers <- cellmarker[cellmarker %in% rownames(scedata[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# Extract expression data and convert to data format suitable for plotting
vln.df <- as.data.frame(scedata[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "expression")

# Remove prefixes "WF-" and "SNJ-" from CB column
vln.df$CB <- gsub("^(WF_|SNJ_)", "", vln.df$CB)

# You can use the following command to confirm whether the prefixes were removed
head(vln.df$CB)

# Continue with the original steps
anno <- scedata@meta.data[, c("CB", "exp","treatment")]
vln.df <- inner_join(vln.df, anno, by = "CB")
# Set gene and treatment order
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)
vln.df$treatment <- factor(vln.df$treatment, levels = c("WF", "SNJ"))  # 👈 Force order

# Plot violin plot with swapped X and Y axes
plot <- vln.df %>%
  ggplot(aes(expression, treatment)) +
  geom_violin(aes(fill = treatment), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # Adjust facet_grid with genes as columns
  scale_fill_manual(values = col) +
  scale_x_continuous("") + scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 28),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    strip.text = element_text(size = 18, face = "bold"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Save violin plot
ggsave(filename = paste(output, "marker_ViolinPlot_ann.pdf", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)
ggsave(filename = paste(output, "marker_ViolinPlot_ann.svg", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)



