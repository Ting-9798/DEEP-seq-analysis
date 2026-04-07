# Clean environment
rm(list = ls())

library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
#install.packages("tidydr")

# Set working directory
setwd("/data/")
# Define folder path
data_dir <- "/data/"

outdir <- "/data/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

folders <- c("WF-7/","WF-14/","WF-PNEC-0/","WF-PNEC-3/")

# Red blood cell marker gene list
HB.genes <- c("Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-y","Hbb-bt","Hbb-bs", "Hbb-bh1", "Hbb-bh2")

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
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # Match red blood cell genes and calculate percentage
  HB_m <- match(HB.genes, rownames(seurat_obj@assays$RNA))
  HB.genes.filtered <- rownames(seurat_obj@assays$RNA)[HB_m]
  HB.genes.filtered <- HB.genes.filtered[!is.na(HB.genes.filtered)]
  
  seurat_obj[["percent.HB"]] <- PercentageFeatureSet(seurat_obj, features = HB.genes.filtered)
  
  # Plot red blood cell gene percentage distribution
  pdf(paste0(OUTPUT,sample_name, "_percent_HB_genes.pdf"), w=5, h=5)
  print(VlnPlot(seurat_obj, features = "percent.HB", pt.size = 0) + 
          ggtitle(paste(sample_name, "Hemoglobin gene percentage")) +
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  # Remove cells with high red blood cell gene percentage
  #seurat_obj <- subset(seurat_obj, subset = percent.HB < 1)
  
  
  # Quality control filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
  
  # Data normalization and highly variable gene finding
  #seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # Assign orig.ident grouping based on folder name
  if (str_detect(sample_name, "WF-7")) {
    seurat_obj$orig.ident <- "7d"
    seurat_obj$treatment <- "7d"
  } else if (str_detect(sample_name, "WF-14")) {
    seurat_obj$orig.ident <- "14d"
    seurat_obj$treatment <- "14d"
  } else if (str_detect(sample_name, "WF-PNEC-0")) {
    seurat_obj$orig.ident <- "0d"
    seurat_obj$treatment <- "0d"
  } else if (str_detect(sample_name, "WF-PNEC-3")) {
    seurat_obj$orig.ident <- "3d"
    seurat_obj$treatment <- "3d"
  }
  
  # Add treatment column
  #seurat_obj$treatment <- sample_name
  
  # Add processed Seurat object to list
  seurat.list[[sample_name]] <- seurat_obj
}


combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)
combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
View(combined_seurat@meta.data)
#summary(combined_seurat@meta.data$treatment)

# Set working directory and save integrated Seurat object
saveRDS(combined_seurat, file = "/out(0+3+7+14)tdt+(1)/combined_seurat.rds")



#######################Seurat analysis#####################
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# Set output directory
setwd("/out(0+3+7+14)tdt+(1)/")
outdir <- "/out(0+3+7+14)tdt+(1)/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("0d", "3d","7d", "14d"))

View(ScRNA@meta.data)
summary(ScRNA$`treatment`)

# Define exp column based on Cbr2 expression level
Cbr2_threshold <- 0  # Assuming expression threshold is 0, can be adjusted based on actual needs
ScRNA@meta.data$exp <- ifelse(ScRNA@assays$RNA@data["tdTomato", ] > Cbr2_threshold, "tdTomato+", "tdTomato-")

# Keep only tdTomato+ cells
ScRNA <- subset(ScRNA, subset = exp == "tdTomato+")
#View(ScRNA@meta.data)
summary(ScRNA$`treatment`)


# Generate violin plot showing QC metrics
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 9, height = 6)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "treatment", pt.size = 0,cols = col)
dev.off()



#### 5. Expression normalization ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Calculate significantly variable genes FindVariableFeatures
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
# Normalization
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


save(ScRNA, file = "ScRNA(before_clustering).RData")



#### 7. Cell clustering and annotation ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")

load("ScRNA(before_clustering).RData")

# Cell clustering
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:60,reduction = "harmony",
          n.neighbors = 15,     # Default 30, decreasing separates clusters
          min.dist = 0.1,       # Default 0.3, decreasing makes clusters more compact and separated
          spread = 2,
          seed.use = 42) %>%       # Increasing spread increases distance between clusters
  RunTSNE(dims = 1:60, reduction = "harmony",         # Use more principal components
          perplexity = 50,      # Default is 30, can be increased appropriately
          theta = 0.3,          # More accurate calculation
          seed.use = 42) %>%
  FindNeighbors(dims = 1:50)


ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))

Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident## Select your resolution based on clustering tree
table(Idents(ScRNA))

# Visualize clustering, displayed in Non-infected and Infected order
pdf(paste(OUTPUT, "split.by_cluster_tsne.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "tsne", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Generate umap plot individually
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

saveRDS(ScRNA, "ScRNA(post_clustering).rds")


##########Check expression proportions of "tdTomato" and "Epcam"###########
#####Combined plotting######
genes <- c("tdTomato", "Epcam")
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


# Groups and genes to plot
treatment_groups <- c("0d", "3d","7d", "14d")
genes <- c("tdTomato", "Epcam")

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
    pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, "_filter.pdf")
    svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, "_filter.svg")
    
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
}



# Visualize expression trends of marker genes for first 6 clusters
output <- paste(outdir,'cell_annotation', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA(post_clustering).rds")
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
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18))  # Increase Y-axis title text size
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=80*dpi,h=7*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): tilt X-axis text
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

# Set output directory
setwd("/out(0+3+7+14)tdt+(1)/")
outdir <- "/out(0+3+7+14)tdt+(1)/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA(post_clustering).rds")
scedata <- readRDS(file_path)
scedata$`treatment` <- factor(scedata$`treatment`, levels = c("0d", "3d","7d", "14d"))

# Define marker genes for different cell types
cellmarker <- c(
  "Epcam", "Krt18", "Cd24", "Krt19",           # Epithelial Cells /EOC epithelial ovarian cancer
  "Pecam1", "Cdh5","Gpihbp1",  # Endothelial Cells
  # Endothelial
  "Pecam1","Kdr","Klf2","Car4",
  # Lymphatic Endothelium
  "Prox1","Lyve1","Pdpn",
  # Fibroblast
  "Col1a1","Col1a2","Pdgfra","Dcn","Acta2",
  # Smooth muscle / Pericytes
  "Acta2","Myh11","Tagln","Pdgfrb","Rgs5",
  #"Lsr","Cxcl10", "Clec4g", "Igfbp7", "Adamts1", "Plpp3", "Iigp1", "Kdr", "Nrp1", "Cyp4b1", "Socs3",  # Endothelial Cell
  #"Ctsd","Ctss",    # Monocyte
  #"Ccl22", "Gm2a", "H2-Ab1",  # Dendritic cells DC
  #"Col1a1", "Col1a2", "Dcn",  # Fibroblasts / Stromal
  "Bgn", "Mgp",  "Sparc",        #Stromal "Pdgfra",
  #"Myh11", "Cnn1", "Smtn",               # Smooth Muscle Cells
  
  "Pcna", "Top2a", "Cdk1",            # Proliferating cells
  #"Birc5",    # Proliferative endothelial cell
  "Pdgfrb", "Rgs5",   # Pericytes
  "Trp63", "Krt5", "Krt14",    # Basal cells
  #"Muc5ac", "Spdef", "Clca1",  # Goblet cells
  "Sulf1",                 # Tuft cells
  "Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells
  "Tppp3","Tubb4b","Tubb1",   #Cilliated
  "Scgb1a1","Scgb3a2",     #club cell
  "Ascl1","Mash1","Calca","Calcb","Uchl1","Syp","Resp18","Pcsk1","Scg5","Chgb",    #PNEC (pulmonary neuroendocrine cells)
  #"Igfbp5","Fmnl2","Pcsk2","Cacna2d1","Chgb","Sez6l2", ### Neuroendocrine cells
  'Ager', 'Hopx', 'Rtkn2', 'Aqp5',   #AT1
  'Lamp3',  'Slc34a2', 'Lpcat1',     #AT2
  "Ndrg1","Cldn4","Krt8","Sprr1a","AW112010"   # DATPs (Damage-associated transition progenitors)
  
  
  #"Cd79a","Cd79b", "Ms4a1",       # B cell
  #"Cd68", "Tspo", "Cd163", "Sepp1","Lgals3","Cd11b", "Apoe",   # Macrophage
  #"Ccl3", "Tmem176a", "Tmem176b", "Arg1",  #Alveolar Macrophages
  #"Tmem176a",  # T helper 17(Th17) cell
  #"Il2rb","Nkg7",       # Natural Killer cell,NK
  #"Trbc2","Icos", # T cells
  #"Cd53","Ptprc","Coro1a",   #Immune cells
  #"Thy1","Il2rb","Il17rb",    #ILC2
  #"B4galt1", "Plcb1", "Brca1", "Mcm6", "Hells", "Cdt1", "Dtl", "Ung", "Rmi2", # NK T cells
  #"Cd8a","Cd8b", # CD8+ T
  #"Gzma", "Gzmb",  "Ifng", "Ccl5", "Il2", "Tbx21", # Cytotoxic T cells
  #"Icos","Grap2", "Csf2", "Gata3" , "Pdcd1" ,               #T helper cell(Th cell)
  #"S100a8", "S100a9",    # Neutrophils
  #"Hba-a1", "Hbb-bs", "Hba-a2", "Hbb-bt",    #Erythrocytes
  #"Tpsab1", "Tpsb2",                                   # Mast cells
  
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Visualize expression of immune cell marker genes using DotPlot
library(ggplot2)
plot <- DotPlot(scedata, features = unique(cellmarker))+
  theme_bw()+theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 12),  # Increase X-axis text size
    axis.text.y = element_text(size = 12),  # Increase Y-axis text size
    legend.title  = element_text(size = 14),
    legend.text = element_text(size = 12)   # Increase legend text size
  ) +
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 12, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 12, height = 8)




############################################################################################
########################### Remove Endothelial Cells and Stromal #############################
# Read existing file
file_path <- file.path(outdir, "ScRNA(post_clustering).rds")
scedata1 <- readRDS(file_path)

# Remove Endothelial Cells and Stromal
meta_filtered <- subset(scedata1@meta.data, !seurat_clusters %in% c("3","6","9","16","18","21","23","25","27","28"))
scedata <- subset(scedata1, cells = rownames(meta_filtered))

# Plot UMAP (by celltype)
pdf(file.path(output, "ann_umap_filtered.pdf"), width = 5.5, height = 5)
DimPlot(scedata, group.by = "seurat_clusters", reduction = 'umap', pt.size = 0.1, label = TRUE, label.size = 5, repel = TRUE, cols = col, label.box = TRUE) +
  theme_dr(xlength = 0.15, ylength = 0.15, arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_blank())
dev.off()

# Plot grouped UMAP (by treatment)
pdf(file.path(output, "ann-diff-umap_filtered.pdf"), width = 6 * length(unique(scedata$treatment)), height = 5)
DimPlot(scedata, reduction = "umap", split.by = "treatment", pt.size = 0.1, label = FALSE, cols = col) +
  theme(strip.text = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))
dev.off()



# Set output directory
setwd("/out(0+3+7+14)tdt+(1)/filtered/")
outdir <- "/out(0+3+7+14)tdt+(1)/filtered/"

# Save filtered data
saveRDS(scedata, file.path(outdir, "celltype_filtered.rds"))

##############Re-run dimensionality reduction and clustering####################

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

file_path <- file.path(outdir, "celltype_filtered.rds")
ScRNA <- readRDS(file_path)
View(ScRNA@meta.data)


#### 6. Normalization and PCA dimensionality reduction ####
# Normalization
ScRNA<-ScaleData(ScRNA)

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


save(ScRNA, file = "ScRNA(before_clustering).RData")



#### 7. Cell clustering and annotation ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")

load("ScRNA(before_clustering).RData")
# Cell clustering
#ScRNA <- ScRNA %>% 
#  RunUMAP(dims = 1:20) %>% 
#  RunTSNE(dims = 1:20) %>%
#  FindNeighbors(dims = 1:20)


# Cell clustering
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:30,reduction = "harmony") %>%       # Increasing spread increases distance between clusters
  RunTSNE(dims = 1:30, reduction = "harmony") %>%
  FindNeighbors(dims = 1:30)


ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))

#library(clustree)
#pdf(paste(OUTPUT, "clustree.pdf"),width=10,height=9)
#clustree(ScRNA)
#dev.off()


Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident## Select your resolution based on clustering tree
table(Idents(ScRNA))

# Visualize clustering, displayed in Non-infected and Infected order
pdf(paste(OUTPUT, "split.by_cluster_tsne.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "tsne", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Generate umap plot individually
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

saveRDS(ScRNA, "ScRNA(post_clustering).rds")


##########Check expression proportions of "tdTomato" and "Epcam"###########
#####Combined plotting######
genes <- c("tdTomato", "Epcam")
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


# Groups and genes to plot
treatment_groups <- c("0d", "3d","7d", "14d")
genes <- c("tdTomato", "Epcam")

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
}


# Visualize expression trends of marker genes for first 6 clusters
output <- paste(outdir,'cell_annotation', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA(post_clustering).rds")
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
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18))  # Increase Y-axis title text size
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=80*dpi,h=7*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  # RotatedAxis(): tilt X-axis text
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

# Set output directory
setwd("/out(0+3+7+14)tdt+(1)/filtered/")
outdir <- "/out(0+3+7+14)tdt+(1)/filtered/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA(post_clustering).rds")
scedata <- readRDS(file_path)

scedata$`treatment` <- factor(scedata$`treatment`, levels = c("0d", "3d","7d", "14d"))

# Define marker genes for different cell types
cellmarker <- c(
  "Epcam", "Krt18", "Cd24", "Krt19",           # Epithelial Cells /EOC epithelial ovarian cancer
  "Pcna", "Top2a", "Cdk1",            # Proliferating cells
  #"Birc5",    # Proliferative endothelial cell
  "Trp63", "Krt5", "Krt14",    # Basal cells
  #"Muc5ac","Muc5b", "Spdef", "Clca1","Gsto1","Tff3",  # Goblet cells
  #"Sulf1",                 # Tuft cells
  "Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells
  "Foxj1","Tppp3","Tubb4b",  #Cilliated
  "Cldn10","Hp","Aldh1a7","Aox3",   # Clara cells	
  "Upk3a","Cyp2f2","N1icd","Scgb1a1","Scgb3a2","Chad",     # club cells   "Upk3a","Cyp2f2","N1icd",
  "Ascl1","Mash1","Calca","Calcb","Uchl1","Syp","Resp18","Pcsk1","Scg5","Chgb",    #PNEC (pulmonary neuroendocrine cells)
  #"Igfbp5","Fmnl2","Pcsk2","Cacna2d1","Chgb","Sez6l2", ### Neuroendocrine cells
  'Ager', 'Hopx', 'Rtkn2', 'Aqp5',"Cav1 ","Spock2",   #AT1
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  #"Lrg1","Lcn2","Retnla","Il33","Car8","Ank3","Cftr",         # Activated AT2
  #"Birc5","Top2a",       # Proliferating AT2s
  "Cldn4","Sfn","Clu","Krt19","Krt8"     # PATs [pre-AT1 transitional state] Epithelial transitional states 
  #"Ndrg1","Cldn4","Krt8","Sprr1a","AW112010"   # DATPs (Damage-associated transition progenitors)
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Visualize expression of immune cell marker genes using DotPlot
library(ggplot2)
plot <- DotPlot(scedata, features = unique(cellmarker))+
  theme_bw()+theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 12),  # Increase X-axis text size
    axis.text.y = element_text(size = 12),  # Increase Y-axis text size
    legend.title  = element_text(size = 14),
    legend.text = element_text(size = 12)   # Increase legend text size
  ) +
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 10, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 10, height = 8)



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
  # Dark blue→green→light green gradient
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq light colors
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


file_path <- file.path(outdir, "ScRNA(post_clustering).rds")
scedata <- readRDS(file_path)

# Annotate cell types for Clusters
scedata <- RenameIdents(scedata, c(
  "0"="AT2",
  "1"="Club cells",
  "2"="Club cells",
  "3"="Club cells",
  "4"="Ciliated",
  "5"="Club cells",
  "6"="Club cells",
  "7"="Club cells",
  "8"="Proliferating cells",
  "9"="Ciliated",
  "10"="Club cells",
  "11"="AT2",
  "12"="Ciliated",
  "13"="Club cells",
  "14"="Ciliated",
  "15"="Basal cells",
  "16"="AT1",
  "17"="AT2",
  "18"="PNEC",
  "19"="AT1",
  "20"="Proliferating cells",
  "21"="Ciliated"
)
)

# Add cell type to meta data
scedata$celltype <- scedata@active.ident
head(scedata@meta.data)

saveRDS(scedata, "celltype.rds")


library(ggsci)
# Plot umap by cell type
pdf(paste(output, "ann_umap.pdf",sep = '/'), width = 5.5, height = 5)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=0.1,label=TRUE,label.size = 5,repel = TRUE,cols=col,
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

#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom"))


# Plot umap by cell type
svg(paste(output, "ann_umap.svg",sep = '/'), width = 5.5, height = 5)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=0.1,label=TRUE,label.size = 5,repel = TRUE,cols=col,
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
#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=6*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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





library(ggsci)
# Plot umap by cell type
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


# Plot umap by cell type
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



########## Add total cell count #########
library(ggsci)

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

# Plot umap by cell type and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 5, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + # Center title
  ggtitle(title_label)
dev.off()

# Plot umap by cell type and save as SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 5, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + # Center title
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

# Set cell type colors, assuming col is a predefined color vector matching cell types
celltype_colors <- col  # Assuming col is a predefined color vector

# Extract UMAP coordinates
umap_coords <- Embeddings(scedata, "umap")  # Extract UMAP coordinates
umap_data <- as.data.frame(umap_coords)
umap_data$celltype <- scedata$celltype  # Add cell type information

# Create data frame with center coordinates for each celltype
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# Calculate center coordinates for each celltype
celltype_centers <- umap_df %>%
  group_by(celltype) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Calculate boundary ranges for each cell type to determine label placement direction
celltype_ranges <- umap_df %>%
  group_by(celltype) %>%
  summarise(
    min_x = min(UMAP1), max_x = max(UMAP1),
    min_y = min(UMAP2), max_y = max(UMAP2),
    width = max_x - min_x, height = max_y - min_y
  )

# Merge center coordinates and boundary information
celltype_labels <- left_join(celltype_centers, celltype_ranges, by = "celltype")

# Determine label placement direction: based on cell cluster shape
celltype_labels <- celltype_labels %>%
  mutate(
    # Determine label placement based on aspect ratio of cell cluster
    direction_x = ifelse(width > height, 1, 0),
    direction_y = ifelse(height > width, 1, 0),
    # Set nudge parameters for labels: push labels outward from cell cluster
    nudge_x = ifelse(UMAP1 > mean(UMAP1), 2, -2),  # Push right clusters left, left clusters right
    nudge_y = ifelse(UMAP2 > mean(UMAP2), 2, -2)   # Push top clusters down, bottom clusters up
  )

# Plot umap by cell type and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1) +
  ggrepel::geom_text_repel(
    data = celltype_labels, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 7, 
    fontface = "bold",
    box.padding = 1.5,        # Increase border padding to move labels away from points
    point.padding = 0.8,      # Increase point padding
    nudge_x = celltype_labels$nudge_x,  # Horizontal offset
    nudge_y = celltype_labels$nudge_y,  # Vertical offset
    min.segment.length = 1,   # Increase minimum segment length
    segment.size = 0.5,       # Line thickness
    segment.color = "grey40", # Line color
    segment.alpha = 0.7,      # Line transparency
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

# Plot umap by cell type and save as SVG
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




##### Add text label borders ####
library(ggsci)
library(Seurat)
library(dplyr)
library(ggrepel)

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

# Set cell type colors, assuming col is a predefined color vector matching cell types
celltype_colors <- col  # Assuming col is a predefined color vector

# Extract UMAP coordinates
umap_coords <- Embeddings(scedata, "umap")  # Extract UMAP coordinates
umap_data <- as.data.frame(umap_coords)
umap_data$celltype <- scedata$celltype  # Add cell type information
umap_data$label <- rownames(umap_data)  # Use row names as labels (can be replaced with other labels)

# Create data frame with center coordinates for each celltype
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# Calculate center coordinates for each celltype
celltype_centers <- umap_df %>%
  group_by(celltype) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot umap by cell type and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype, label = label)) +
  geom_point(size = 0.1) +
  ggrepel::geom_label_repel(
    data = celltype_centers, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),  # Use corresponding colors
    size = 7, fontface = "bold",   # Set label to bold
    box.padding = 0.5,             
    point.padding = 0.5,           
    segment.color = "grey30",
    label.size = 0.5,             # Set text box border width
    label.r = 0.3,                # Set rounded corners
    label.border = "black",       # Set text box border to black
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +  # Set cell type colors
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20)) + # Center title
  ggtitle(title_label)  # Add total cell count as title
dev.off()

# Plot umap by cell type and save as SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype, label = label)) +
  geom_point(size = 0.1) +
  ggrepel::geom_label_repel(
    data = celltype_centers, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),  # Use corresponding colors
    size = 7, fontface = "bold",   # Set label to bold
    box.padding = 0.5,             
    point.padding = 0.5,           
    segment.color = "grey30",
    label.size = 0.5,             # Set text box border width
    label.r = 0.3,                # Set rounded corners
    label.border = "black",       # Set text box border to black
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +  # Set cell type colors
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20)) + # Center title
  ggtitle(title_label)  # Add total cell count as title
dev.off()







##### Add cell counts for each group #####
# Count cell numbers for each treatment
cell_counts <- scedata@meta.data %>%
  group_by(treatment) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(treatment, " ( n = ", n, " cells)"))

# Construct named vector for facet label replacement
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# Plot PDF
pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),
    width = 3.5*length(unique(scedata$treatment)), height = 10.5)

DimPlot(scedata, reduction = "umap", split.by = "treatment",
        pt.size = 0.1, label = FALSE, cols = col) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # Subplot title bold
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
    width = 3.5*length(unique(scedata$treatment)), height = 10.5)

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





###########Plot dot plot for annotation
library(ggh4x)

# Define marker genes for different cell types
cellmarker <- c(
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  "Cldn10","Hp","Aldh1a7","Aox3",   # Clara cells	
  "Scgb1a1","Scgb3a2","Chad",     # club cells
  "Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells
  "Foxj1","Tppp3","Tubb4b",  #Cilliated
  "Pcna", "Top2a", "Cdk1",            # Proliferating cells
  "Trp63", "Krt5", "Krt14",    # Basal cells
  'Ager', 'Hopx', 'Rtkn2', "Cav1 ","Spock2",   #AT1
  "Ascl1","Mash1","Calca","Calcb","Uchl1","Syp","Resp18","Pcsk1","Scg5","Chgb"    #PNEC (pulmonary neuroendocrine cells)
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Generate DotPlot
p <- DotPlot(scedata, features = unique(cellmarker))

# Get data
dat <- p$data

# Get cluster annotations and merge
anno <- distinct(data.frame(id = scedata$celltype, celltype = scedata$seurat_clusters))
colnames(anno) <- c("id", "celltype")
df <- left_join(dat, anno, by = "id")

# Define cluster order
cluster.order <- c(0,2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
df$celltype <- factor(df$celltype, levels = cluster.order)

# Plot dot plot
p <- ggplot(df, aes(features.plot, interaction(celltype, id), size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 16),    # Increase legend text size
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

# Save images
ggsave(plot = p, filename = paste(output, "ann_DotPlot.pdf",sep = '/'), width = 12, height = 6)
ggsave(plot = p, filename = paste(output, "ann_DotPlot.svg",sep = '/'), width = 12, height = 6)


#########Plot heatmap for cell annotation############
# change annotation color
library("scales")
library(ggsci)
library(scRNAtoolVis)

# Set annotation colors
mycol1 <- pal_simpsons()(18)

col <- c('#FF6666','#E5D2DD',"#BC8F8F",'#FFCC99','#FF9999','#4F6272','#58A4C3',"#CC0066",
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#FFFFCC")


col <- setNames(c(
  "#FF6666", "#E5D2DD", "#6A4C93", "#FFCC99", "#FF9999",
  "#58A4C3", "#CC0066", "#F3B1A0", "#66CCCC", "#CCFFCC",
  "#6699CC", "#CCCCFF", "#BC8F8F", "#FF3366", "#99FFFF",
  "#00CC66", "#FF9933", "#FFFFCC", "#9999FF"
), c(
  "Monocytes", "Dendritic cells", "T cells", "Fibroblasts", 
  "Ciliated epithelial cells", "Endothelial Cells", "Macrophages", 
  "B cells", "AT2", "NK", "Smooth Muscle Cells", 
  "Epithelial Cells", "pDCs", "Pericytes", "Neutrophils", 
  "Plasma Cells", "Lymphatic endothelial cells", "Mast cells", 
  "Club cells"
))

all(scedata$celltype %in% names(col))

#file_path <- file.path(outdir, "celltype.rds")
#scedata <- readRDS(file_path)

pdf(file = paste(output, "ann_Heatmap.pdf",sep = '/'), width = 7, height = 6)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # Custom high value color
dev.off()

svg(file = paste(output, "ann_Heatmap.svg",sep = '/'), width = 7, height = 6)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # Custom high value color
dev.off()




#########Expression trends of marker genes in cells#######
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# Select differential cells for display
output <- paste(outdir,'cluster', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
#file_path <- file.path(outdir, "ScRNA(post_clustering).rds")
ScRNA <- readRDS(file_path)

ScRNA_3d <- subset(ScRNA, subset = treatment == "3d")

# Save results
saveRDS(ScRNA_3d, file = file.path(outdir, "celltype_3d.rds"))


# Define marker genes for different cell types
cellmarker <- c(
  
  "Scgb1a1",    # club cells
  'Lamp3',    #AT2
  "Dnah5",
  "Pcna",            # Proliferating cells
  "Krt5",  # Basal cells
  'Ager'   #AT1
)

cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# Create lists to store RidgePlot, VlnPlot, and FeaturePlot
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
  
  # Use numeric vector for sorting, avoid directly indexing ScRNA object
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # Extract sorted cell names
  
  # FeaturePlot plotted in sorted order using continuous color gradient
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", 
                                       cells = cell_names_ordered,  # Specify cell order
                                       ncol = 1) + 
    scale_color_gradientn(colors = c("#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366")) +  # Set continuous color gradient
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
pdf(paste0(output, "/marker_FeaturePlot_umap.pdf"), width = 15, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

svg(paste0(output, "/marker_FeaturePlot_umap.svg"), width = 15, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

pdf(paste0(output, "/marker_VlnPlot_umap.pdf"), width = 18, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()
svg(paste0(output, "/marker_VlnPlot_umap.svg"), width = 18, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()



######## PNEC #########
# Set T cell activation related genes
cellmarker <- c(
  "Ascl1","Mash1","Calca","Calcb","Syp","Resp18","Pcsk1","Scg5","Chgb"      #PNEC (pulmonary neuroendocrine cells)
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# Create lists to store RidgePlot, VlnPlot, and FeaturePlot
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
  
  # Use numeric vector for sorting, avoid directly indexing ScRNA object
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # Extract sorted cell names
  
  # FeaturePlot plotted in sorted order using continuous color gradient
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", 
                                       cells = cell_names_ordered,  # Specify cell order
                                       ncol = 1) + 
    scale_color_gradientn(colors = c('#E5D2DD',  "#FF3366")) +  # Set continuous color gradient
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
pdf(paste0(output, "/spacial_FeaturePlot_umap.pdf"), width = 20, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 4))
dev.off()

svg(paste0(output, "/spacial_FeaturePlot_umap.svg"), width = 20, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 4))
dev.off()

pdf(paste0(output, "/spacial_VlnPlot_umap.pdf"), width = 20, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =4))
dev.off()
svg(paste0(output, "/spacial_VlnPlot_umap.svg"), width = 20, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =4))
dev.off()




######## PNEC #########
# Set T cell activation related genes
cellmarker <- c(
  "Trp63","Calca","Mki67"             #"tdTomato","Epcam"   
)

cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# Create lists to store RidgePlot, VlnPlot, and FeaturePlot
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
  
  # Use numeric vector for sorting, avoid directly indexing ScRNA object
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # Extract sorted cell names
  
  # FeaturePlot plotted in sorted order using continuous color gradient
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", 
                                       cells = cell_names_ordered,  # Specify cell order
                                       ncol = 1) + 
    scale_color_gradientn(colors = c("#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366")) +  # Set continuous color gradient
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
pdf(paste0(output, "/spacial_FeaturePlot_umap11.pdf"), width = 15, height = 4)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

svg(paste0(output, "/spacial_FeaturePlot_umap11.svg"), width = 15, height = 4)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

pdf(paste0(output, "/spacial_VlnPlot_umap11.pdf"), width = 12, height = 3)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()
svg(paste0(output, "/spacial_VlnPlot_umap11.svg"), width = 12, height = 3)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()





###### Positive cells #########
cellmarker <- c("tdTomato","Epcam")
cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# Create data frame to store average expression
avg_expr_data <- data.frame(Gene = character(), Avg_Expression = numeric(), stringsAsFactors = FALSE)

# Loop to calculate average expression and plot
for (gene in cellmarker) {
  
  # Get gene expression data
  gene_expr <- FetchData(ScRNA, vars = gene)
  
  # Ensure gene_expr is a vector, not a data frame
  gene_expr_vec <- gene_expr[[gene]]
  
  # Calculate average expression for this gene
  avg_expression <- mean(gene_expr_vec)
  
  # Calculate number of expressing cells based on fixed threshold 1.5
  expressed_cells <- sum(gene_expr_vec > 0)
  total_cells <- nrow(ScRNA@meta.data)
  expression_ratio <- expressed_cells / total_cells * 100
  
  # Record average expression
  avg_expr_data <- rbind(avg_expr_data, data.frame(Gene = gene, Avg_Expression = avg_expression))
  
  # Sort by expression level
  ordered_cells <- order(gene_expr_vec, decreasing = FALSE)
  cell_names_ordered <- rownames(ScRNA@meta.data)[ordered_cells]
  
  # Set title with expression ratio annotation
  plot_title <- paste0(gene, " Expression (", round(expression_ratio, 2), "%)")
  
  # Add gene expression information to Seurat object
  ScRNA[[paste0(gene, "_expr")]] <- gene_expr_vec
  
  # Plot FeaturePlot
  p <- FeaturePlot(
    ScRNA, 
    features = gene,
    reduction = "umap", 
    cells = cell_names_ordered,  
    ncol = 1,
    cols = c("#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366")
  ) +
    ggtitle(plot_title) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) + 
    NoAxes() 
  
  # Save PDF
  ggsave(paste0(output, "/", gene, "_FeaturePlot_umap.pdf"), plot = p, width = 4, height = 4)
  
  # Save SVG
  ggsave(paste0(output, "/", gene, "_FeaturePlot_umap.svg"), plot = p, width = 4, height = 4)
}

# Save average expression to TXT file
write.table(avg_expr_data, file = paste0(output, "/Tcell_activation_genes_avg_expression.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




####### Calculate cell proportions ###########

col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#094867",'#E59CC4',"#6666CC",
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue→green→light green gradient
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

######## Calculate cell counts for different cell populations across all samples
cell_counts <- as.data.frame(table(Idents(scedata)))
colnames(cell_counts) <- c("CellType", "Counts")

# Sort by counts in descending order
cell_counts <- cell_counts[order(-cell_counts$Counts), ]
# Save cell counts for all cell populations to specified directory
write.csv(cell_counts, paste(output, "cell_counts.csv", sep='/'), row.names = FALSE)

# Select top 11 cell types and save to specified directory
cell_counts_top9 <- head(cell_counts, 11)
write.csv(cell_counts_top9, paste(output, "cell_counts_top9.csv", sep='/'), row.names = FALSE)

# Load required packages
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

# Save image as png format
ggsave(paste(output, "cell_type_distribution.pdf", sep='/'), plot = p, width = 7, height = 6, dpi = 800)
ggsave(paste(output, "cell_type_distribution.svg", sep='/'), plot = p, width = 7, height = 6, dpi = 800)


# Calculate cell counts for different cell populations by group
# Calculate counts for each cell type by sample group
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

# Sort
cell_counts_group$Sample <- factor(cell_counts_group$Sample, levels = c("0d", "3d","7d", "14d"))

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
# Add count text labels
p <- p + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


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
# Add proportion text labels
p <- p + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)



############Grouping############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate proportion of each cell type within each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

# Sort
cell_counts_treatment$Treatment <- factor(cell_counts_treatment$Treatment, levels = c("0d", "3d","7d", "14d"))

########## Plot stacked bar chart for cell counts ##########
p1 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 22),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)

########## Plot stacked bar chart for cell proportions ##########
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
ggsave(file_path, plot = p2, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)





#################Plot line area chart for cell proportions###############
# Replace with your provided color scheme

col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#094867",'#E59CC4',"#6666CC",
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # Dark blue→green→light green gradient
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq light colors
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


# Plot line area chart for counts
p_counts_area <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, group = CellType)) +
  stat_summary(geom = 'line', fun = 'mean', color = 'white', linewidth = 1) +
  geom_area(aes(fill = CellType)) +
  scale_fill_manual(values = col) +
  labs(x = NULL, y = "Counts") +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.01, 0.01), name = "Counts") +  # ✅ Use correct scale_y and add Y axis title
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  geom_vline(aes(xintercept = "3d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "7d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "14d"), linetype = "dashed", size = 1, colour = "white")

ggsave(paste0(output, "/genecount_treatment_area.pdf"), plot = p_counts_area, width = 6, height = 4, dpi = 800)


# cell_counts_treatment$Ratio should already exist
p_ratio_area <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, group = CellType)) +
  stat_summary(geom = 'line', fun = 'mean', color = 'white', linewidth = 1) +
  geom_area(aes(fill = CellType)) +
  scale_fill_manual(values = col) +
  labs(x = NULL, y = "Ratio") +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.01, 0.01), name = "Ratio") +  # ✅ Use correct scale_y and add Y axis title
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  geom_vline(aes(xintercept = "3d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "7d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "14d"), linetype = "dashed", size = 1, colour = "white")

ggsave(paste0(output, "/geneRatio_treatment_area.pdf"), plot = p_ratio_area, width = 6, height = 4, dpi = 800)







###############Plot each sample individually########################
library(ggplot2)
library(dplyr)

# Define samples to plot
samples_to_plot <- c("0d", "3d","7d", "14d")

# Create output directory (if it doesn't exist)
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

for (sample_id in samples_to_plot) {
  
  # Filter data for current sample
  cell_counts_group_filtered <- subset(cell_counts_group, Sample == sample_id)
  
  ### —— Cell count plot —— ###
  # Aggregate and create legend labels (cell counts)
  cell_counts_group_agg <- aggregate(Counts ~ CellType, cell_counts_group_filtered, sum)
  cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", cell_counts_group_agg$Counts, ")")
  legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)
  
  # Plot count bar chart
  p_counts <- ggplot(cell_counts_group_filtered, aes(x = Sample, y = Counts, fill = CellType)) + 
    geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
    theme_classic() +
    labs(x = '', y = 'Counts') +
    scale_fill_manual(values = col, labels = legend_labels_group) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5),
          axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          legend.title = element_blank(),
          legend.text = element_text(size = 20))
  
  # Save cell count plot
  ggsave(paste0(output, "/genecount_", sample_id, ".pdf"), plot = p_counts, width = 6, height = 8, dpi = 800)
  ggsave(paste0(output, "/genecount_", sample_id, ".svg"), plot = p_counts, width = 6, height = 8, dpi = 800)
  
  ### —— Cell proportion plot —— ###
  # Aggregate and create legend labels (proportions)
  cell_counts_group_agg <- aggregate(Ratio ~ CellType, cell_counts_group_filtered, sum)
  cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", round(cell_counts_group_agg$Ratio * 100, 2), "%)")
  legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)
  
  # Plot proportion bar chart
  p_ratio <- ggplot(cell_counts_group_filtered, aes(x = Sample, y = Ratio, fill = CellType)) + 
    geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
    theme_classic() +
    labs(x = '', y = 'Ratio') +
    scale_fill_manual(values = col, labels = legend_labels_group) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5),
          axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          legend.title = element_blank(),
          legend.text = element_text(size = 20))
  
  # Save cell proportion plot
  ggsave(paste0(output, "/geneRatio_", sample_id, ".pdf"), plot = p_ratio, width = 6.5, height = 8, dpi = 800)
  ggsave(paste0(output, "/geneRatio_", sample_id, ".svg"), plot = p_ratio, width = 6.5, height = 8, dpi = 800)
}






################Batch differential expression analysis################
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


# Set output directory
setwd("/out(0+3+7+14)tdt+(1)/filtered/")
outdir <- "/out(0+3+7+14)tdt+(1)/filtered/"

# Create output directory
output <- file.path(outdir, "differential_analysis")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# Read data
file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)

logFCfilter <- 0.25
adjPvalFilter <- 0.05

genes <- c("Gpnmb","Spp1","Ctsd","Nfkb1","Ifngr1","Camk4","Zeb1","Rora","Cd14","Il1b",
           "Il12rb2","Cxcl2","Fth1","Icos","Stat4","Cd63","Thbs1","Tyrobp")

timepoints <- c("3d", "7d", "14d")

for (tp in timepoints) {
  comp_name <- paste0(tp, "_vs_0d")
  
  markers <- FindMarkers(object = scRNAsub,
                         ident.1 = tp,
                         ident.2 = "0d",
                         group.by = "treatment",
                         logfc.threshold = 0,
                         min.pct = 0.25,
                         test.use = "wilcox")
  
  markers$gene <- rownames(markers)
  markers <- markers %>%
    mutate(Significance = ifelse(p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter, 
                                 ifelse(avg_log2FC > 0, "Up", "Down"), "Normal"))
  
  write.table(markers, file = file.path(output, paste0("sig.markers_", comp_name, ".txt")),
              sep = "\t", row.names = TRUE, quote = FALSE)
  saveRDS(markers, file = file.path(output, paste0("ScRNA.sig.markers_", comp_name, ".rds")))
  
  up_df <- markers %>% filter(Significance == "Up")
  down_df <- markers %>% filter(Significance == "Down")
  write.csv(up_df, file = file.path(output, paste0("upregulated_genes_", comp_name, ".csv")), row.names = TRUE)
  write.csv(down_df, file = file.path(output, paste0("downregulated_genes_", comp_name, ".csv")), row.names = TRUE)
  
  
  # Calculate numbers of upregulated and downregulated genes
  upregulated_genes <- sum(markers$Significance == "Up")
  downregulated_genes <- sum(markers$Significance == "Down")
  total_diff_genes <- upregulated_genes + downregulated_genes
  
  # Save data frames for upregulated and downregulated genes separately
  upregulated_genes_df <- markers %>%
    filter(Significance == "Up")
  downregulated_genes_df <- markers %>%
    filter(Significance == "Down")
  
  # Filter top 10 upregulated and downregulated genes for labeling
  top_genes_upregulated <- upregulated_genes_df %>%
    filter(p_val_adj < 0.01 & avg_log2FC > 0) %>%
    arrange(p_val_adj) %>%
    head(15)
  top_genes_downregulated <- downregulated_genes_df %>%
    filter(p_val_adj < 0.01 & avg_log2FC < 0) %>%
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
    scale_x_continuous(limits = c(-3, 2), breaks = seq(-3, 2, by = 1)) +
    scale_y_continuous(limits = c(0, max(-log10(markers$p_val_adj), na.rm = TRUE) + 1))+
    #scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 5)) +
    theme(plot.title = element_text(size = 22, face = "bold", hjust = 0), 
          legend.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 20, hjust = 0.5),
          axis.text = element_text(size = 18))
  
  
  ggsave(file.path(output, paste0(comp_name, "_volcano_plot.svg")), p, width = 8, height = 7)
  ggsave(file.path(output, paste0(comp_name, "_volcano_plot.pdf")), p, width = 8, height = 7)
  
  # ==== GSEA analysis ====
  deg <- markers[, c("avg_log2FC", "p_val_adj")]
  colnames(deg) <- c("log2FoldChange", "pvalue")
  gene <- bitr(rownames(deg), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  gene$logFC <- deg$log2FoldChange[match(gene$SYMBOL, rownames(deg))]
  geneList <- gene$logFC
  names(geneList) <- gene$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  kk_gse <- gseKEGG(geneList = geneList, organism = "mmu", nPerm = 1000, minGSSize = 10,
                    pvalueCutoff = 0.25, verbose = FALSE)
  kk_gse@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kk_gse@result$Description)
  kk_gse <- DOSE::setReadable(kk_gse, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID')
  
  write.csv(as.data.frame(kk_gse), file = file.path(output, paste0(comp_name, "_kk_gse.csv")))
  
  kk_gse_cut <- kk_gse[kk_gse$pvalue < 0.05 & kk_gse$p.adjust < 0.25 & abs(kk_gse$NES) > 1, ]
  kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0, ]
  kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0, ]
  
  if (nrow(kk_gse_cut_up) > 0) {
    gseap_up <- gseaplot2(kk_gse, kk_gse_cut_up$ID[1],
                          title = kk_gse_cut_up$Description[1],
                          color = c("#FF4500", "#32CD32"),
                          base_size = 30, rel_heights = c(1.5, 0.5, 1), subplots = 1:3)
    ggsave(file.path(output, paste0("GSEA_", comp_name, "_up.pdf")), gseap_up, width = 14, height = 10)
    ggsave(file.path(output, paste0("GSEA_", comp_name, "_up.svg")), gseap_up, width = 14, height = 10)
  }
  
  if (nrow(kk_gse_cut_down) > 0) {
    gseap_down <- gseaplot2(kk_gse, kk_gse_cut_down$ID[1],
                            title = "DOWN_GSEA",
                            color = c("#FF4500", "#32CD32"),
                            base_size = 30, rel_heights = c(1.5, 0.5, 1), subplots = 1:3)
    ggsave(file.path(output, paste0("GSEA_", comp_name, "_down.pdf")), gseap_down, width = 14, height = 10)
    ggsave(file.path(output, paste0("GSEA_", comp_name, "_down.svg")), gseap_down, width = 14, height = 10)
  }
  
  # ridgeplot
  ridgep <- ridgeplot(kk_gse, showCategory = 15, fill = "pvalue", core_enrichment = TRUE) +
    ggtitle(paste("Ridgeplot -", comp_name))
  ggsave(file.path(output, paste0("ridgeplot_GSEA_", comp_name, ".pdf")), ridgep, width = 10, height = 8)
  ggsave(file.path(output, paste0("ridgeplot_GSEA_", comp_name, ".svg")), ridgep, width = 10, height = 8)
  
  # ==== GO / KEGG enrichment ====
  gene_up <- up_df$gene
  gene_down <- down_df$gene
  
  # Convert SYMBOL to ENTREZID
  gene_up_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                               keys = gene_up, 
                                                               columns = 'ENTREZID', 
                                                               keytype = 'SYMBOL')[,2]))
  gene_down_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                                 keys = gene_down, 
                                                                 columns = 'ENTREZID', 
                                                                 keytype = 'SYMBOL')[,2]))
  
  # Perform GO enrichment analysis
  go_up <- enrichGO(gene = gene_up_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)
  go_down <- enrichGO(gene = gene_down_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)
  
  go_up@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", go_up@result$Description)
  go_down@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", go_down@result$Description)
  
  write.csv(as.data.frame(go_up), file = file.path(output, paste0("go_up_", comp_name, ".csv")))
  write.csv(as.data.frame(go_down), file = file.path(output, paste0("go_down_", comp_name, ".csv")))
  
  ggsave(file.path(output, paste0("go_up_dot_", comp_name, ".pdf")), dotplot(go_up) + ggtitle("GO Up"), width = 6, height = 5)
  ggsave(file.path(output, paste0("go_down_dot_", comp_name, ".pdf")), dotplot(go_down) + ggtitle("GO Down"), width = 6, height = 5)
  
  kegg_up <- enrichKEGG(gene = gene_up_entrez, organism = 'mmu')
  kegg_down <- enrichKEGG(gene = gene_down_entrez, organism = 'mmu')
  
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
  
  write.csv(as.data.frame(kegg_up), file = file.path(output, paste0("kegg_up_", comp_name, ".csv")))
  write.csv(as.data.frame(kegg_down), file = file.path(output, paste0("kegg_down_", comp_name, ".csv")))
  
  ggsave(file.path(output, paste0("kegg_up_dot_", comp_name, ".pdf")), dotplot(kegg_up) + ggtitle("KEGG Up"), width = 6, height = 5)
  ggsave(file.path(output, paste0("kegg_down_dot_", comp_name, ".pdf")), dotplot(kegg_down) + ggtitle("KEGG Down"), width = 6, height = 5)
  
}






###############Differential gene expression##################
library(Seurat)
library(tidyverse)
library(ggsci)

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# Load data
ScRNA <- readRDS("celltype.rds")

# Set T cell activation related genes for analysis
genes <- c("Gpnmb","Spp1","Ctsd","Nfkb1","Ifngr1","Camk4","Zeb1","Rora","Cd14","Il1b","Il12rb2","Cxcl2","Fth1","Icos","Stat4","Cd63","Thbs1","Tyrobp")


# Filter genes that actually exist in the expression matrix
genes <- genes[genes %in% rownames(ScRNA)]

# Extract treatment information
ScRNA$treatment <- as.factor(ScRNA@meta.data$treatment)
treatment_groups <- levels(ScRNA$treatment)

# Get expression matrix
expr_matrix <- GetAssayData(ScRNA, slot = "data")[genes, ]


avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment")$RNA
avg_expr_selected <- avg_expr[, treatment_groups]

# Save average expression table
avg_expr_df <- avg_expr_selected %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Gene")
write.table(avg_expr_df, file = paste0(output, "/differential_gene_expression.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)


# Plot dot plot for gene expression (grouped by treatment)
plot <- DotPlot(ScRNA, features = unique(genes), group.by = "treatment") + 
  RotatedAxis() +
  coord_flip() +
  scale_color_gradientn(colors = c('#CCCCCC', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18))

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.pdf", sep='/'), plot = plot, width = 5, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.svg", sep='/'), plot = plot, width = 5, height = 5)















#####################Pseudotime analysis monocle2 ###############################
library(Seurat)
library(monocle)
library(igraph)
#devtools::install_version("igraph", version = "2.0.8", repos = "http://cran.us.r-project.org")

#BiocManager::install("monocle")
#install.packages("igraph")
dpi=300

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272',"#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")

# Set output directory
setwd("/out(0+3+7+14)tdt+(1)/filtered/")
outdir <- "/out(0+3+7+14)tdt+(1)/filtered/"

output <- paste(outdir,'monocle2(notch)', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data1 <- readRDS(file_path)
summary(data1$treatment)

# Remove 0d
data <- subset(data1, subset = treatment != "0d")
summary(data$treatment)

saveRDS(data, "celltype(no-0d).rds")


## Extract raw expression matrix and sparsify: UMI count
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

## Build monocle2 object
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

# Set target gene list (Notch pathway related)
target_genes <- c("Rest", "Numbl", "Numb", "Dlk2", "Dlk1", "Dll3", "Rbpj", "Mamld1", "Maml3", "Maml2", "Maml1",
                  "Rfng", "Mfng", "Lfng", "Pofut2", "Pofut1", "Jag2", "Jag1", "Dll4", "Dll1", 
                  "Nrarp", "Heyl", "Hey2", "Hey1", "Hes7", "Hes6", "Hes5", "Hes3", "Hes2", "Hes1",
                  "Notch4", "Notch3", "Notch2", "Notch1")

# Filter target genes that actually exist in the expression matrix
ordering_genes <- intersect(target_genes, rownames(expr_matrix))

cds <- setOrderingFilter(cds, ordering_genes)

## Since there are few genes, use all genes for ordering gene selection
#cds <- setOrderingFilter(cds, rownames(expr_matrix))

# Visualize ordering genes
pdf(paste0(output,"/ordergenes.pdf"))
plot_ordering_genes(cds)
dev.off()
ggsave(paste0(output,"/ordergenes.png"),plot_ordering_genes(cds),width = dpi*6, height = dpi*6, units = "px",type='cairo')


## Dimensionality reduction
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)


######## Trajectory visualization
## Pseudotime represents pseudotime value, State represents cell state, celltype represents cell type
## Visualize trajectory by different types
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
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # Set gradient color
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # For discrete data (State), use color vector col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # Use col as colors
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # For discrete data (celltype and treatment), use custom colors custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # Use custom_colors as colors
      theme(legend.text = element_text(size = 16),  # Adjust legend text size
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  # Save as PDF
  ggsave(filename = paste(output, paste0("monocle_", type, ".pdf", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 10, height = 5)
  
  # Save as SVG
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 10, height = 5)
}


#saveRDS(cds,  file = file.path(output, "monocle2.rds"))
saveRDS(cds,"monocle2.rds")


#####Calculate cell proportions after pseudotime#####

# Extract state information for each cell
state_data <- pData(cds)$State
celltype_data <- pData(cds)$celltype

# Calculate counts of different cell types in each state
cell_counts_state <- as.data.frame(table(celltype_data, state_data))
colnames(cell_counts_state) <- c("CellType", "State", "Counts")

# Calculate proportion of different cell types in each state
cell_counts_state$Ratio <- ave(cell_counts_state$Counts, cell_counts_state$State, FUN = function(x) x / sum(x))

# Plot bar chart of cell type proportions by state
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

# Save images
ggsave(paste0(output, "/state_celltype_proportion.pdf"), plot = p, width = 7)
ggsave(paste0(output, "/state_celltype_proportion.svg"), plot = p, width = 7)








############## Merge multiple subclusters back to main clusters ###############

col <- c("#CC0066","#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272',"#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")

setwd("/out(0+3+7+14)tdt+(1)/filtered/")
outdir <- "/out(0+3+7+14)tdt+(1)/filtered/"

output <- paste(outdir,"merge_subclusters", sep='/')
dir.create(output)

# Read data
#DCs <- readRDS('D:/R/GS/HZ/20250324-脾/out1/DCs/celltype.rds')  # Subcluster-annotated subpopulation
PNEC <- readRDS('/out(0+3+7+14)tdt+(1)/filtered/PNEC/celltype(Notch2).rds')  # Subcluster-annotated subpopulation
ScRNA <- readRDS("/out(0+3+7+14)tdt+(1)/filtered/celltype.rds")  # Main population with all cells


# Set identifier for each subcluster
#Idents(DCs) <- "celltype"  # Subcluster information for DCs cells
Idents(PNEC) <- "exp"  # Subcluster information for T cells
Idents(ScRNA) <- "celltype"  # Main population information for all cells

# Extract subclusters

all <- subset(ScRNA)
table(all@meta.data$celltype)

# Merge multiple subcluster identifiers into main population
Idents(all, cells = colnames(PNEC)) <- Idents(PNEC)

# Create new meta information column to store merged cell population information
all$celltype_merged <- Idents(all)

# Update Idents
Idents(all) <- "celltype_merged"

# View merged cluster statistics
table(all@meta.data$celltype_merged)

# Save merged data as rds file
saveRDS(all, file = "celltype_merged.rds")

# Plot UMAP by cell type and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 10, height = 6)
DimPlot(object = all, group.by = "celltype_merged", reduction = 'umap', pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        axis.title.x = element_text(size = 16, face = "bold"),  # Increase X-axis title size
        axis.title.y = element_text(size = 16, face = "bold"),  # Increase Y-axis title size
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18),
        plot.title = element_blank())
dev.off()

# Plot UMAP by cell type and save as SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 10, height = 6)
DimPlot(object = all, group.by = "celltype_merged", reduction = 'umap', pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        axis.title.x = element_text(size = 16, face = "bold"),  # Increase X-axis title size
        axis.title.y = element_text(size = 16, face = "bold"),  # Increase Y-axis title size
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18),
        plot.title = element_blank())
dev.off()


# Plot umap by cell type
pdf(file = file.path(output, "ann_umap1.pdf"), width = 6, height = 4)
DimPlot(object=all,group.by = "treatment",reduction='umap',pt.size=0.1,label=FALSE,label.size = 6,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

dev.off() 

#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom"))


# Plot umap by cell type
svg(file = file.path(output, "ann_umap1.svg"), width = 6, height = 4)
DimPlot(object=all,group.by = "treatment",reduction='umap',pt.size=0.1,label=FALSE,label.size = 6,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())
#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom")) +
dev.off()


# Plot UMAP grouped by treatment and save as PDF
pdf(paste(output, "ann-diff-umap.pdf", sep = '/'), width = 13, height = 5)
DimPlot(all, reduction = "umap", split.by = "treatment", pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col)+
  theme(
    # Increase label text size
    axis.text.x = element_text(size = 16),  # X-axis label size
    axis.text.y = element_text(size = 16),  # Y-axis label size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
    axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18)    # Increase legend text size
  )
dev.off()

# Plot UMAP grouped by treatment and save as SVG
svg(paste(output, "ann-diff-umap.svg", sep = '/'), width = 13, height = 5)
DimPlot(all, reduction = "umap", split.by = "treatment", pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col)+
  theme(
    # Increase label text size
    axis.text.x = element_text(size = 16),  # X-axis label size
    axis.text.y = element_text(size = 16),  # Y-axis label size
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase X-axis title size
    axis.title.y = element_text(size = 18, face = "bold"),  # Increase Y-axis title size
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # Increase title size
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18)    # Increase legend text size
  )
dev.off()









#############################Cell-cell communication comparison between different time points (0d, 3d, 7d, 14d)###########################
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(NMF)
library(ComplexHeatmap)
library(cowplot)
library(gridExtra)
library(grid)

options(stringsAsFactors = FALSE)

# Define colors
col <- c('#437eb8','#FF6666','#FFCC99','#FF9999',"#FFCCCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFFFCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272',"#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")


# Read data
file_path <- file.path(outdir, "celltype.rds")
seuratdata <- readRDS(file_path)
head(seuratdata@meta.data)


# Set output directory
output <- paste(outdir,"cellchat", sep='/')
dir.create(output, recursive = TRUE)
setwd(output)


# -----------------------------
#seuratdata <- subset(seuratdata,subset = celltype != "Proliferating cells")

# Clean factor levels (very important)
seuratdata@meta.data$celltype <- droplevels(seuratdata@meta.data$celltype)

# Check if removed
table(seuratdata@meta.data$celltype)

# Check time points in treatment
print("Available treatment levels:")
print(unique(seuratdata@meta.data$treatment))

# Extract data for each time point
time_points <- c("3d", "7d", "14d")
seurat_subset <- list()

for (tp in time_points) {
  seurat_subset[[tp]] <- subset(seuratdata, subset = treatment %in% tp)
  print(paste("Number of cells in", tp, ":", ncol(seurat_subset[[tp]])))
}

# Create list of CellChat objects
cellchat_list <- list()

for (tp in time_points) {
  cat("Processing", tp, "...\n")
  
  # Create CellChat object
  cellchat_obj <- createCellChat(
    object = seurat_subset[[tp]]@assays$RNA@data, 
    meta = seurat_subset[[tp]]@meta.data, 
    group.by = "celltype"
  )
  
  # Check and remove unused factor levels
  cellchat_obj@idents <- droplevels(cellchat_obj@idents)
  print(paste("Cell types in", tp, ":"))
  print(levels(cellchat_obj@idents))
  
  # Set database
  cellchat_obj@DB <- CellChatDB.mouse
  
  # Analysis pipeline
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- computeCommunProb(cellchat_obj, raw.use = TRUE, population.size = TRUE)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
  
  # Save individual object
  saveRDS(cellchat_obj, file = paste0("cellchat_", tp, ".rds"))
  
  cellchat_list[[tp]] <- cellchat_obj
}

# Merge CellChat objects from all time points
cellchat <- mergeCellChat(
  cellchat_list, 
  add.names = names(cellchat_list), 
  cell.prefix = TRUE
)
saveRDS(cellchat, "cellchat.rds")

# If reading saved objects
# cellchat_list <- list()
# for (tp in time_points) {
#   cellchat_list[[tp]] <- readRDS(paste0("cellchat_", tp, ".rds"))
# }
# cellchat_merged <- readRDS("cellchat_merged.rds")

# Comparative analysis - interaction counts and strength
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = 1:3, 
                           measure = "count", color.use = col) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  ggtitle("Number of Interactions")

gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = 1:3, 
                           measure = "weight", color.use = col) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  ggtitle("Interaction Strength")

p <- gg1 + gg2
ggsave("Overview_interaction_number_strength.pdf", p, width = 10, height = 5)



################### Ligand-receptor comparative analysis##########################
## Display differences for all ligand-receptor pairs
levels(cellchat@idents$joint)
#levels(cellchat_obj@idents)
#levels(cellchat_normal@idents)


##### DCs #####
p <- netVisual_bubble(cellchat, sources.use = c(7), targets.use = c(1,2,3,4,5,6),comparison = c(1, 2,3),angle.x = 45)
# Increase legend text and title size
p <- p + theme(
  legend.title = element_text(size = 18),  # Legend title size
  legend.text = element_text(size = 16),   # Legend text size
  plot.title = element_text(size = 18),    # Title size
  axis.title = element_text(size = 16),    # Axis title size
  axis.text = element_text(size = 16)      # Axis tick text size
)
ggsave("Compare_LR_bubble(PNEC).pdf", p, width = 12, height = 60,limitsize = FALSE)




sig.interactions <- subsetCommunication(cellchat)
unique(sig.interactions$interaction_name)



pairLR.use <- as.data.frame(c(
  "LAMA5_ITGA6_ITGB1", 
  "FN1_ITGA5_ITGB1",
  "FN1_ITGAV_ITGB1",
  "HBEGF_EGFR",
  "HBEGF_EGFR_ERBB2",
  "AREG_EGFR",
  "AREG_EGFR_ERBB2",
  "APP_CD74",
  "DLL1_NOTCH1",
  "DLL1_NOTCH2",
  "JAG1_NOTCH1",
  "WNT4_FZD1_LRP5",
  "WNT4_FZD1_LRP6",
  "WNT7B_FZD1_LRP5",
  "WNT7B_FZD1_LRP6",
  "SEMA3A_NRP1_PLXNA1",
  "SEMA3A_NRP1_PLXNA2",
  "SEMA3A_NRP1_PLXNA4",
  "PTN_SDC1","PTN_SDC4",
  "THBS1_SDC1","THBS1_SDC4"
))


colnames(pairLR.use) <- 'interaction_name'

# Use netVisual_bubble to plot specific ligand-receptor pair bubble plot
p <- netVisual_bubble(cellchat, pairLR.use = pairLR.use, sources.use = c(7), targets.use = c(1,2,3,4,5,6),comparison = c(1, 2,3),angle.x = 45)
p <- p + theme(
  legend.title = element_text(size = 18),  # Legend title size
  legend.text = element_text(size = 16),   # Legend text size
  plot.title = element_text(size = 18),    # Title size
  axis.title = element_text(size = 16),    # Axis title size
  axis.text = element_text(size = 16)      # Axis tick text size
)

# Save as PDF
ggsave("Compare_LR_bubble_Selected_pairs(PNEC).pdf", p, width = 10, height = 8)



# Network analysis and heatmap analysis
png("Diff_interaction_count.png", width = 800, height = 800)
par(cex = 2) 
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
dev.off()

png("Diff_interaction_weight.png", width = 800, height = 800)
par(cex = 2) 
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

p1_img <- ggdraw() + draw_image("Diff_interaction_count.png")
p2_img <- ggdraw() + draw_image("Diff_interaction_weight.png")
combined_plot <- plot_grid(p1_img, p2_img, ncol = 2)
ggsave("Diff_combined_interaction.pdf", combined_plot, width = 8, height = 4)


# Heatmap analysis
png("Diff_heatmap_count.png", width = 900, height = 800, res = 150)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex = 1.5)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

png("Diff_heatmap_weight.png", width = 900, height = 800, res = 150)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex = 1.5)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

p1_img <- ggdraw() + draw_image("Diff_heatmap_count.png")
p2_img <- ggdraw() + draw_image("Diff_heatmap_weight.png")
combined_plot <- plot_grid(p1_img, p2_img, ncol = 2)
ggsave("Diff_combined_heatmap.pdf", combined_plot, width = 12, height = 5)



## Network plots comparing interaction counts between groups
# Get maximum weight
weight.max <- getMaxWeight(cellchat_list, attribute = c("org.idents", "count"))
# Loop to plot network plots comparing interaction counts between groups and save as PDF
for (i in 1:length(cellchat_list)) {
  # Define PDF filename
  pdf_filename <- paste0("circle_compare_Number_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Number of interactions - ", names(cellchat_list)[i])
  # Set graphics parameters
  par(oma = c(0, 0, 0, 0))
  par(cex = 2)
  # Plot interaction network
  netVisual_circle(cellchat_list[[i]]@net$count, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 10)
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off() 
}

## Network plots comparing interaction strength between groups
weight.max.weight <- getMaxWeight(cellchat_list, attribute =c("org.idents", "weight"))
for (i in 1:length(cellchat_list)) {
  # Define PDF filename
  pdf_filename <- paste0("circle_compare_Weight_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Weight of interactions - ", names(cellchat_list)[i])
  # Set graphics parameters
  par(oma = c(0, 0, 0, 0))
  par(cex = 2)
  # Plot interaction network
  netVisual_circle(cellchat_list[[i]]@net$weight, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.weight.max = weight.max.weight[2], 
                   edge.width.max = 30)
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off() 
}


#########Plot chord diagrams##########
library(circlize)

# Get maximum weight
weight.max <- getMaxWeight(cellchat_list, attribute = c("org.idents", "count"))

# Loop to plot chord diagrams comparing interaction counts between groups and save as PDF
for (i in 1:length(cellchat_list)) {
  # Define PDF filename
  pdf_filename <- paste0("chord_compare_Number_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Number of interactions - ", names(cellchat_list)[i])
  
  # Get interaction count data
  interaction_data <- cellchat_list[[i]]@net$count
  
  # Set chord diagram parameters
  chordDiagram(interaction_data, 
               preAllocate = 1, 
               direction.type = c("diffHeight"), 
               link.arr.type = "big.arrow", 
               annotationTrack = c("name", "grid"),  # Keep only one definition
               link.border = NA,
               transparency = 0.5)
  
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off() 
}



# Get maximum weight
weight.max.weight <- getMaxWeight(cellchat_list, attribute =c("org.idents", "weight"))

# Loop to plot chord diagrams comparing interaction strength between groups and save as PDF
for (i in 1:length(cellchat_list)) {
  # Define PDF filename
  pdf_filename <- paste0("chord_compare_Weight_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Weight of interactions - ", names(cellchat_list)[i])
  
  # Get interaction strength data
  interaction_weight_data <- cellchat_list[[i]]@net$weight
  
  # Set chord diagram parameters
  chordDiagram(interaction_weight_data, 
               preAllocate = 1, 
               direction.type = c("diffHeight"), 
               link.arr.type = "big.arrow", 
               annotationTrack = c("name", "grid"),
               link.border = NA,
               transparency = 0.5)
  
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off() 
}



## Identification and visualization of conserved and specific signaling pathways
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strength.pdf", p, width = 10, height = 10)


saveRDS(cellchat, "cellchat_compare.rds")

## Comparison of cell signaling patterns
pathway.union <- Reduce(union, list(cellchat_list[[1]]@netP$pathways, 
                                    cellchat_list[[2]]@netP$pathways))

ht1 = netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cellchat_list)[1], width = 6, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cellchat_list)[2], width = 6, height = 22)

ht1_grob <- grid.grabExpr(draw(ht1))
ht2_grob <- grid.grabExpr(draw(ht2))

combined_plot <- grid.arrange(ht1_grob, ht2_grob,ncol = 2, 
                              widths = unit.c(unit(12, "cm"), unit(11, "cm"))) 
pdf("combined_heatmap_all.pdf", width = 10, height = 12)
grid.draw(combined_plot)
dev.off()

# outgoing
ht1 = netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cellchat_list)[1], width = 6, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cellchat_list)[2], width = 6, height = 22)

ht1_grob <- grid.grabExpr(draw(ht1))
ht2_grob <- grid.grabExpr(draw(ht2))

combined_plot <- grid.arrange(ht1_grob, ht2_grob,ncol = 2, 
                              widths = unit.c(unit(12, "cm"), unit(11, "cm")))
pdf("combined_heatmap_outgoing.pdf", width = 10, height = 12)
grid.draw(combined_plot)
dev.off()

# incoming
ht1 = netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cellchat_list)[1], width = 6, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cellchat_list)[2], width = 6, height = 22)

ht1_grob <- grid.grabExpr(draw(ht1))
ht2_grob <- grid.grabExpr(draw(ht2))

combined_plot <- grid.arrange(ht1_grob, ht2_grob, ncol = 2, 
                              widths = unit.c(unit(12, "cm"), unit(11, "cm")))

pdf("combined_heatmap_incoming.pdf", width = 10, height = 12)
grid.draw(combined_plot)
dev.off()


# Specific signaling pathway comparison
df.net <- subsetCommunication(cellchat_obj)
table(df.net$pathway_name)  ### View specific signaling pathways
#levels(df.net$source)


# Set list of signaling pathways of interest
pathways_list <- c( "TGFb", "BMP", "EGF", "FGF", "WNT", "GAS", "THBS", "NOTCH", "VEGF", "COLLAGEN", "LAMININ")


# Loop through signaling pathway list and generate plots
# Loop through signaling pathway list and generate plots
for (pathways.show in pathways_list) {
  # Check if pathways.show exists in pathways of cellchat_list
  available_pathways <- unique(unlist(lapply(cellchat_list, function(x) x@netP$pathways)))
  if (!(pathways.show %in% available_pathways)) {
    message(paste("Skipping pathway:", pathways.show, "- not found in any cellchat_list"))
    next
  }
  
  # Ensure weight.max uses correct parameters
  weight.max <- tryCatch({
    getMaxWeight(cellchat_list, slot.name = "netP", attribute = pathways.show)
  }, error = function(e) {
    message(paste("Error in getMaxWeight for pathway:", pathways.show))
    return(NULL)
  })
  
  weight.max <- getMaxWeight(cellchat_list, slot.name = c("netP"), attribute = pathways.show)
  
  # Loop to generate gene expression heatmaps, contribution analysis plots, heatmaps
  for (i in 1:length(cellchat_list)) {
    SampleOutput <- paste0("compare_", pathways.show, "_", names(cellchat_list)[i])
    
    # Generate and save individual gene expression heatmap (PNG format)
    png(file = paste0(SampleOutput, "_GeneExpression_heatmap.png"), width = 1200, height = 1400, res = 300)
    p3 <- plotGeneExpression(cellchat_list[[i]], signaling = pathways.show)
    print(p3)
    dev.off()
    
    # Generate and save individual contribution analysis plot (PNG format)
    png(file = paste0(SampleOutput, "_contribution.png"), width = 1200, height = 800, res = 300)
    p4 <- netAnalysis_contribution(cellchat_list[[i]], signaling = pathways.show)
    print(p4)
    dev.off()
    
    # Generate and save individual heatmap visualization (PNG format)
    png(file = paste0(SampleOutput, "_netVisual_heatmap.png"), width = 1200, height = 1200, res = 300)
    p5 <- netVisual_heatmap(cellchat_list[[i]], signaling = pathways.show, color.heatmap = "Reds")
    print(p5)
    dev.off()
  }
  
  # Combine gene expression heatmaps
  gene_expr_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_GeneExpression_heatmap.png"))
  gene_expr_imgs <- lapply(gene_expr_filenames, function(x) ggdraw() + draw_image(x))
  combined_gene_expr_plot <- plot_grid(plotlist = gene_expr_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_GeneExpression_heatmap.pdf"), combined_gene_expr_plot, width = 10, height = 8)
  
  # Combine contribution analysis plots
  contribution_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_contribution.png"))
  contribution_imgs <- lapply(contribution_filenames, function(x) ggdraw() + draw_image(x))
  combined_contribution_plot <- plot_grid(plotlist = contribution_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_contribution.pdf"), combined_contribution_plot, width = 8, height = 8)
  
  # Combine heatmap visualizations
  heatmap_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_netVisual_heatmap.png"))
  heatmap_imgs <- lapply(heatmap_filenames, function(x) ggdraw() + draw_image(x))
  combined_heatmap_plot <- plot_grid(plotlist = heatmap_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_netVisual_heatmap.pdf"), combined_heatmap_plot, width = 8, height = 8)
  
  
  # Loop to plot interaction network diagrams and save as PDF
  for (i in 1:length(cellchat_list)) {
    # Define PDF filename
    pdf_filename <- paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_net.pdf")
    pdf(pdf_filename, width = 10, height = 10)  # Set PDF file size
    
    plot_title <- paste0(pathways.show, " - ", names(cellchat_list)[i])  # Set plot title
    # Set graphics parameters
    par(oma = c(0, 0, 0, 0))
    par(cex = 2)
    
    # Plot interaction network diagram
    netVisual_aggregate(cellchat_list[[i]], signaling = pathways.show, layout = "circle", 
                        edge.weight.max = weight.max[1], edge.width.max = 30,
                        vertex.label.cex = 1.2)
    
    title(main = plot_title, line = 1, cex.main = 1.5)  # Set title
    dev.off()  # Close PDF device
  }
  
  # Loop to generate chord diagrams
  for (i in 1:length(cellchat_list)) {
    png_filename <- paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_chord.png")
    png(png_filename, width = 800, height = 800)
    par(oma = c(0, 0, 0, 0))  
    par(cex = 2)  
    netVisual_aggregate(cellchat_list[[i]], signaling = pathways.show, layout = "chord", 
                        pt.title = 3, title.space = 0.05, signaling.name = paste(pathways.show, names(cellchat_list[i])),
                        vertex.label.cex = 1, font.main = 2)
    dev.off()
  }
  
  # Combine chord diagrams
  chord_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_chord.png"))
  chord_imgs <- lapply(chord_filenames, function(x) ggdraw() + draw_image(x))
  combined_chord_plot <- plot_grid(plotlist = chord_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_chord.pdf"), combined_chord_plot, width = 8, height = 8)
}




