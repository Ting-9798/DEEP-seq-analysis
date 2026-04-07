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
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 10, min.features = 1000)
  
  # Add mitochondrial gene proportion
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Calculate ribosomal protein proportion
  seurat_obj[["percent.rps"]] <- PercentageFeatureSet(seurat_obj, pattern = c("^RPS","^RPL"))
  
  #seurat_obj <- subset(seurat_obj, subset = percent.rps < 10)
  
  # Quality control filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 25)
  
  
  # Data normalization and highly variable gene finding
  #seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  
  # Add treatment column
  seurat_obj$treatment <- sample_name 
  
  
  # Add processed Seurat object to list
  seurat.list[[sample_name]] <- seurat_obj
}

# Integrate data
#anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)
#combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:20)

combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)
combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
#View(combined_seurat@meta.data)

# Save merged Seurat object
saveRDS(combined_seurat, file = "/out(SNJ1)/combined_seurat.rds")


#######################Seurat Analysis#####################
# Set output directory
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("/out(SNJ1)//")
outdir <- "/out(SNJ1)/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# Splice the full path
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
#View(ScRNA@meta.data)

# Generate violin plot showing QC metrics
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 6)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"), ncol = 4, group.by = "treatment", pt.size = 0,cols = col)
dev.off()



#### 5. Expression normalization ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Calculate highly variable genes using FindVariableFeatures
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



save(ScRNA, file = "ScRNA（after batch correction before clustering）.RData")


#### 7. Cell clustering and annotation ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487",  "#094867","#FEDD81", "#FF9A84", '#E59CC4',"#6666CC",
  "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue → green → light green gradient
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq light colors
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


load("ScRNA（after batch correction before clustering）.RData")
# Cell clustering
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:20) %>% 
  RunTSNE(dims = 1:20) %>%
  FindNeighbors(dims = 1:20)



ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))


#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident## Select your resolution based on clustering tree
table(Idents(ScRNA))

# Ensure "treatment" factor levels are ordered as Non-infected and Infected
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Display clustering, ordered by Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Display clustering, ordered by Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 40, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "orig.ident", cols = col)
dev.off()

# Generate separate umap plot
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size=14),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20))


DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
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



# Generate separate umap plot
pdf(paste(OUTPUT, "cluster_umap_11.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size=14),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20))
dev.off()



saveRDS(ScRNA, "ScRNA（after clustering）.rds")




# Display expression trends of marker genes for top 6 clusters
output <- paste(outdir,'Cell_annotation', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
ScRNA <- readRDS(file_path)

# Find Marker genes
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




########### Manual cell annotation ########

col <- c('#3399CC',"#CCCCCC",'#E5D2DD','#FF9999',"#66CCCC","#A4CDE1","#FFCCCC","#CCFFCC","#FFFFCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", "#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("/out(SNJ1)//")
outdir <- "/out(SNJ1)/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
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

# Create a new column named 'celltype', classify based on TPRX1 gene expression
scedata$celltype <- ifelse(TPRX1_expression > 0, "8CLC", "non-8CLC")
#View(scedata@meta.data)

# Plot cell type umap plot
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

######View tdt expression########
# Get tdTomato gene expression data
tdTomato_expr <- FetchData(scedata, vars = "TPRX1")

# Add expression information to cell metadata
scedata$tdTomato_expr <- tdTomato_expr$TPRX1

# Calculate tdTomato expression ratio (proportion of cells with non-zero expression)
expressed_cells <- sum(scedata$tdTomato_expr > 0)
total_cells <- nrow(scedata@meta.data)
expression_ratio <- expressed_cells / total_cells * 100

# Sort cell data by expression level
cells_ordered <- scedata@meta.data[order(scedata$tdTomato_expr, decreasing = FALSE), ]

# Extract ordered cell names
cell_names_ordered <- rownames(cells_ordered)

# Set title, annotating expression ratio
plot_title <- paste0("TPRX1 Expression (", round(expression_ratio, 2), "%)")

# Draw FeaturePlot in ordered sequence
pdf(paste0(output, "/TPRX1_FeaturePlot_umap.pdf"), width = 4, height = 4)
FeaturePlot(
  scedata, 
  features = "TPRX1",
  reduction = "umap", 
  cells = cell_names_ordered,  # Plot in cell order
  ncol = 1,
  cols = c("#CCCCCC",'#3399CC')
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
  cells = cell_names_ordered,  # Plot in cell order
  ncol = 1,
  cols = c("#CCCCCC",'#3399CC')
) +
  ggtitle(plot_title) +  # Add title
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + 
  NoAxes()  # Remove axes
dev.off()



####### Calculate cell proportion ###########
col <- c("#CCCCCC",'#3399CC','#E5D2DD','#FF9999',"#66CCCC","#A4CDE1","#FFCCCC","#CCFFCC","#FFFFCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", "#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype_8C.rds")
scedata <- readRDS(file_path)
View(scedata@meta.data)

table(scedata$celltype)

# Calculate cell counts for different cell groups across each group
# Calculate count of each cell type by sample group
cell_counts_group <- as.data.frame(table(scedata$orig.ident,scedata$celltype))
colnames(cell_counts_group) <- c("Sample", "CellType","Counts")

cell_counts_group$CellType <- factor(cell_counts_group$CellType, levels = c( "non-8CLC","8CLC"))

# Sort by "CellType"
cell_counts_group <- cell_counts_group %>%
  arrange(CellType)

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

file_path <- paste0(output, "/genecount_8c.pdf")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount_8c.svg")
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

file_path <- paste0(output, "/geneRatio_8c.pdf")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio_8c.svg")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 8, dpi = 800)



cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)


cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Visualize immune cell marker gene expression using DotPlot
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

# Save DotPlot figure
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 14, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 14, height = 5)



#########Draw cell annotation heatmap############
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
               markerGene = cellmarker)   # Custom high-value colors
dev.off()

svg(file = paste(output, "ann_Heatmap_11.svg",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # Custom high-value colors
dev.off()



#########Draw marker gene expression violin box plot

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


# Filter marker genes that exist in the dataset
existing_markers <- cellmarker[cellmarker %in% rownames(scedata[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# Extract expression data and convert to data format suitable for plotting
vln.df <- as.data.frame(scedata[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

# Continue with original steps
anno <- scedata@meta.data[, c("CB", "seurat_clusters")]
vln.df <- inner_join(vln.df, anno, by = "CB")
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)


# Draw Violin Plot, swap X and Y axes
plot <- vln.df %>%
  ggplot(aes(exp, seurat_clusters)) +
  geom_violin(aes(fill = seurat_clusters), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # Adjust facet_grid, use genes as columns
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

# Save Violin Plot
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

col <- c("#CCCCCC",'#3399CC','#E5D2DD','#FF9999',"#66CCCC","#A4CDE1","#FFCCCC","#CCFFCC","#FFFFCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", "#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
scedata <- readRDS(file_path)

# Annotate cell types for Clusters
scedata <- RenameIdents(scedata, c(
  "0"="non-8CLC",
  "1"="non-8CLC", 
  "2"="non-8CLC",
  "3"= "non-8CLC",
  "4"="non-8CLC",
  "5"="8CLC",
  "6"="non-8CLC",
  "7"="non-8CLC",
  "8"= "non-8CLC",
  "9"="non-8CLC",
  "10"="non-8CLC",
  "11"="non-8CLC")
)

# Add cell type to meta data
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
# Plot cell type umap plot
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

#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom"))


# Plot cell type umap plot
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
#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=6*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # Increase subplot title font size
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
    strip.text = element_text(size = 18, face = "bold"),  # Increase subplot title font size
    axis.text.x = element_text(size = 16),  # X-axis label size
    axis.text.y = element_text(size = 16),  # Y-axis label size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
    axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18)    # Increase legend text size
  )
dev.off()



#########Draw marker gene expression violin box plot

library(ggplot2)
library(reshape2)
library(dplyr)

# Filter marker genes that exist in the dataset
existing_markers <- cellmarker[cellmarker %in% rownames(ScRNA[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# Extract expression data and convert to data format suitable for plotting
vln.df <- as.data.frame(ScRNA[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

# Continue with original steps
anno <- ScRNA@meta.data[, c("CB", "celltype")]
vln.df <- inner_join(vln.df, anno, by = "CB")
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)


# Draw Violin Plot, swap X and Y axes
plot <- vln.df %>%
  ggplot(aes(exp, celltype)) +
  geom_violin(aes(fill = celltype), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # Adjust facet_grid, use genes as columns
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

# Save Violin Plot
ggsave(filename = paste(output, "marker_ViolinPlot_ann.pdf", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)
ggsave(filename = paste(output, "marker_ViolinPlot_ann.svg", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)



####### Calculate cell proportion ###########
col <- c("#CCCCCC",'#3399CC','#E5D2DD','#FF9999',"#66CCCC","#A4CDE1","#FFCCCC","#CCFFCC","#FFFFCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", "#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

table(scedata$seurat_clusters)

######## Calculate cell counts for different cell groups across all samples
cell_counts <- as.data.frame(table(Idents(scedata)))
colnames(cell_counts) <- c("CellType", "Counts")

# Sort by cell count descending
cell_counts <- cell_counts[order(-cell_counts$Counts), ]
# Save cell counts for all cell groups to specified directory
write.csv(cell_counts, paste(output, "cell_counts.csv", sep='/'), row.names = FALSE)

# Select top 11 cell types and save to specified directory
cell_counts_top9 <- head(cell_counts, 11)
write.csv(cell_counts_top9, paste(output, "cell_counts_top9.csv", sep='/'), row.names = FALSE)

# Load required packages
library(ggplot2)
# Draw bar plot
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

# Save image as png format
ggsave(paste(output, "cell_type_distribution.pdf", sep='/'), plot = p, width = 7, height = 6, dpi = 800)
ggsave(paste(output, "cell_type_distribution.svg", sep='/'), plot = p, width = 7, height = 6, dpi = 800)


# Calculate cell counts for different cell groups across each group
# Calculate count of each cell type by sample group
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



############ Grouping ############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate proportion of each cell type within each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

########## Draw cell count stacked bar plot ##########
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

########## Draw cell proportion stacked bar plot ##########
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

#######################Seurat Analysis#####################
# Set output directory
col <- c("#CCCCCC",'#3399CC','#E5D2DD','#FF9999',"#66CCCC","#A4CDE1","#FFCCCC","#CCFFCC","#FFFFCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", "#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


output <- paste(outdir,'8C', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data <- readRDS(file_path)

# Select epithelial cells and cells belonging to "Tumor", "Res", or "Sen" groups
Cells.sub <- subset(data@meta.data, celltype == c("8CLC"))
summary(Cells.sub$celltype)
scedata <- subset(data, cells=row.names(Cells.sub))


# Define exp column based on Cbr2 expression
Cbr2_threshold <- 0  # Assuming expression threshold is 0, can be adjusted based on actual needs
scedata@meta.data$exp <- ifelse(scedata@assays$RNA@data["TPRX1", ] > Cbr2_threshold, "TPRX1+", "TPRX1-")

# Check definition results
#View(scedata@meta.data)

library(ggsci)
# Plot cell type umap plot
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

# Plot cell type umap plot (SVG format)
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





####### Calculate cell proportion ###########

output <- paste(outdir,'8C', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype(8C).rds")
scedata <- readRDS(file_path)

table(scedata$exp)

# Calculate cell counts for different cell groups across each group
# Calculate count of each cell type by sample group
cell_counts_group <- as.data.frame(table(scedata$orig.ident, scedata$exp))
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
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


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
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 4*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


