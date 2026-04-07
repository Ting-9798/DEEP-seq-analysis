####################### Seurat Analysis #####################
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

# Splice the full path
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))
#ScRNA$`orig.ident` <- factor(ScRNA$`orig.ident`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))

## Calculate red blood cell proportion
ScRNA[["percent.hb"]] <- PercentageFeatureSet(ScRNA, pattern = "^Hbb-|^Hba-")

# Generate violin plot showing QC metrics
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 5)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4, group.by = "treatment", pt.size = 0,cols = col)
dev.off()
## Remove red blood cell contamination
ScRNA <- subset(ScRNA, subset = percent.hb < 1)


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

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#6A4C93',
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA（after batch correction before clustering）.RData")


ScRNA <- ScRNA %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,spread = 2) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)

ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))

#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.0.7"
ScRNA$seurat_clusters <- ScRNA@active.ident## Select your resolution based on clustering tree
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



########### Manual cell annotation ########

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
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
  "PECAM1", "VWF", "CDH5",             # Endothelial Cells
  #"PROX1","LYVE1", "CCL21",  # Lymphatic endothelial cells 
  "TPSAB1", "TPSB2",                   # Mast cells
  "PCNA", "TOP2A", "CDK1",            # Proliferating cells
  "COL1A1", "COL1A2", "DCN",                    # Fibroblasts / Stromal
  #"ACTA2" , "MYH11", "CNN1",                           # Smooth muscle cells
  "SFRP4","FAP","MMP11","PDGFRA","PDPN","FSP1",   # CAFs
  "PDGFRB", "RGS5",   # Pericytes
  
  #"MALAT1",    # Tfh (T follicular helper)
  
  "KRT19", "MUC1", "KRT8","CDH1",         # Epithelial Cells / EOC
  "EPCAM","MUC16","WFDC2","KRT7","PAX8","CLDN4","CD24","TP53","WT1",      # Cancer Cells
  
  #"HBA1","HBA2","HBB",    # Erythroid Cells
  #"CD79A", "MS4A1","IGKC","IGHA1",                 # B cells
  "MZB1","JCHAIN",   # Plasma Cells
  #"GZMB","LILRA4","SPIB","ZFAT",         # pDCs
  "HLA-DPB1", "HLA-DRA", "HLA-DRB1", # Dendritic cells
  "CD68", "APOE", "LGALS3", "ITGAM", "PPARG" ,         # Macrophages
  "FCN1","CD300E", "NLRP3" ,"TBC1D8",               # Monocytes
  #"ITGAM","CD14","CD33","S100A8","S100A9","ARG1","NOS2","IL4R","IL10","ANXA1","LOX1",           # MDSCs
  "S100A8", "S100A9", "G0S2",     # Neutrophils
  
  "GZMA", "GZMB",  "IFNG", "CCL4","CCL5",  # Cytotoxic T cells
  #"NKG7", "CCL5", "KLRB1","GZMA","KLRF1", "PRF1",     # Natural Killer (NK) cells
  "H3F3A","ABCF1","INTS6","TRNAU1AP","ERO1B","ERCC1","CEMIP2",      # NK T cells   "BRCA1", "MCM6", "HELLS", "CDT1", "DTL",
  #"CD8A","CD8B", # CD8+ T 
  "TRBC1","TRBC2","CCL5","CD2", "CD3E",  "CD3G"    # T cells
  #"CCR7",  "CD3D","CD3E","CD4","CD8A","SELL","TCF7","LEF1",  # Naive T Cells
  
  #### Lung cancer ####
  #"PECAM1", "FLT1", "VWF", "CDH5", "CA4",              # Endothelial Cells
  #"PROX1", "FLT4", "PDPN", "LYVE1", "CCL21", "ITGA9", "NRP2", "MRC1", "SOX18", "FOXC2", # Lymphatic endothelial cells 
  #"CALB2","WT1","MSLN","UPK3B","PDPN",     # Mesothelial Cells
  #"PDGFRB", "CSPG4", "RGS5", "ANPEP", "ABCC9", "NOTCH3", "KCNJ8", # Pericytes
  #"PSD3", "FOXI1", "CFTR", "ATP6V0D2", "ATP6V1B1", "ASCL3", "TP63", "EZR",  # Pulmonary ionocytes
  #"CTSD", "CTSS", "ITGAL", "RBM47","PLXDC2", "FCN1","C5AR1","CD300E", "NLRP3" ,"TBC1D8",  # Monocytes
  #"TPSAB1", "TPSB2",                                   # Mast cells
  #"PCNA", "TOP2A", "CCNA2", "MCM5", "CDK1",            # Proliferating cells
  #"HLA-DPB1", "HLA-DRA", "NAAA", "GM2A", "HLA-DRB1", "PPT1", "CYTIP", # Dendritic cells
  #"COL1A1", "COL1A2", "DCN", "FN1",                    # Fibroblasts / Stromal
  #"MYH11", "CNN1", "SMTN",                             # Smooth muscle cells
  #"HBA1","HBA2","HBB",    # Erythroid Cells
  #"CD79A", "CD79B", "MS4A1",                           # B cells
  #"CD38","MZB1","IRF4","JCHAIN","SDC1",   # Plasma Cells
  #"GZMB","IL3RA","CLEC4C","LILRA4","TCF4","IRF7","PTCRA","SPIB","NRP1","IRF8",           # pDCs
  #"CD68", "LGALS3", "ITGAM", "APOE",                   # Macrophages
  #"MARCO","ABCA1","MRC1","CCL3", "FBP1",  "MCEMP1",       # Alveolar Macrophages
  #"ARG1","MRC1",  ## M2 Macrophages
  #"IL2RB", "NKG7", "ADAMTS14", "KLRA4", "KLRF1", "PRF1", # Natural Killer (NK) cells
  #"H4C3","SMC4","HMGN2","HMGB2","NCALD","KCNQ5","PDGFD","ZNF331","B4GALT1", "PLCB1", "BRCA1", "MCM6", "HELLS", "CDT1", "DTL", "UNG", "RMI2", # NK T cells
  #"CD2", "CD3E", "CD3G", "CCR7", "NKG7", "KLRB1", "CD3G", "TRBC2", # T cells
  #"IL2RA","FOXP3","CTLA4","IKZF2","TNFRSF18" ,     # Treg 
  #"ICOS", "GRAP2", "CSF2", "GATA3", "PDCD1",           # T helper cells (Th cells)
  #"CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "NKG7","CCL4", "CCL5", "IL2", "TBX21", # Cytotoxic T cells
  #"CCR7",  "CD3D","CD3E","CD4","CD8A","SELL","TCF7","LEF1",  # Naive T Cells
  #"SNAP25", "SYP", "SYT1", "DCX", "SYN1", "SLC17A7",   # Neurons
  #"S100A8", "S100A9", "G0S2", "LY6G6C", "MPO", "CSF3R","BCL2A1","IL1R2","FCGR3B","CSF3R",    # Neutrophils
  #"SPINK5","NXN","TIMP1","KRT5","TP63","KRT14","ITGA6","ITGB4","NGFR",   # Basal Cells
  #"FOXJ1", "TPPP3", "TUBB4B",  "TP73", "CCDC7", # Ciliated cells
  #"SCGB1A1", "SCGB3A2", "CALCB", "NOTCH2", "HES1",     # Club cells
  #"CALCA", "SYP", "RESP18", "PCSK1", "SCG5", "CHGB", "SEZ6L2",   # PNEC
  #"FMNL2","PCSK2","CACNA2D1",     ### Neuroendocrine cells
  #"EPCAM", "KRT18", "CD24", "KRT19", "FGFR2", "SPRR3","SPRR2A","SCEL",          # Epithelial Cells /EOC
  #"STEAP4","CEACAM6","SCGB1A1","MUC5B","MUC5AC","SPDEF","FOXJ1",    # Secretory Epithelial Cells
  #"FOXC1", "TP73", "DNAH5", "DNAH9", "CCDC39", "CCDC40", "TEKT1", "RSPH4A", "TUBB4B", "SPEF2", # Ciliated epithelial cells
  #'AGER', 'CAV1', 'CLIC5' ,'HOPX', 'SEMA3E', 'COL4A3',  #AT1
  #"LAMP3", "ABCA3", "SLC34A2", "LPCAT1", "SFTPC", "SFTPA1", "SFTPB", "SFTPD", "AGER", "CLDN18", "NKX2-1", "MUC1", "KRT8" # AT2 cells
  
)




cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Visualize immune cell marker gene expression using DotPlot
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

# Save DotPlot figure
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 12, height = 6)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 12, height = 6)




library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1",'#FF9999',"#66CCCC",'#4F6272',"#FF3366","#CC0066","#00CC66","#CC99CC","#FFCCCC","#9999FF","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#FF6699","#6699CC","#FFFFCC")

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

# Annotate cell types for Clusters
scedata <- RenameIdents(scedata, c(
  "0"="Macrophages",
  "1"= "DCs", 
  "2"="Fibroblasts",
  "3"= "CTL",
  "4"="Fibroblasts",
  "5"= "Fibroblasts",
  "6"=  "DCs", 
  "7"="Pericytes", 
  "8"=  "Endothelial Cells",
  "9"="Fibroblasts",
  "10"="Cancer Cells", 
  "11"="Cancer Cells",  
  "12"="Fibroblasts",
  "13"="Fibroblasts",
  "14"="DCs", 
  "15"="Monocytes",
  "16"="Cancer Cells",  
  "17"= "Fibroblasts",
  "18"= "Plasma Cells",
  "19"="Proliferating Cells", 
  "20"="Endothelial Cells",
  "21"="Cancer Cells"
  
  
)
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


file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

library(ggsci)
# Plot cell type umap plot
pdf(paste(output, "ann_umap.pdf",sep = '/'), width = 7, height = 6)
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


# Plot cell type umap plot
svg(paste(output, "ann_umap.svg",sep = '/'), width = 7, height = 6)
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

# Filter required samples
#selected_samples <- c("T1", "T2", "T3", "T4")

# Filter data with orig.ident as T1, T2, T3, T4 from scedata
#filtered_data <- subset(scedata, subset = orig.ident %in% selected_samples)

# Plot and save as PDF
#pdf(paste(output, "ann-diff-umap-sample-selected.pdf", sep = '/'), width = 22, height = 5)
#DimPlot(filtered_data, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, 
#        label = FALSE, label.size = 5, repel = TRUE, cols = col)
#dev.off()

# Plot and save as SVG
#svg(paste(output, "ann-diff-umap-sample-selected.svg", sep = '/'), width = 22, height = 5)
#DimPlot(filtered_data, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, 
#        label = FALSE, label.size = 5, repel = TRUE, cols = col)
#dev.off()





library(ggsci)
# Plot cell type umap plot
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


# Plot cell type umap plot
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



########## Add overall cell count ##########
library(ggsci)

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

# Plot cell type umap plot and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 5, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + # Center the title
  ggtitle(title_label)
dev.off()

# Plot cell type umap plot and save as SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 5, repel = TRUE, cols = col) +
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

# Set cell type colors, assuming col is a predefined color vector
celltype_colors <- col  # Assuming col is a predefined color vector

# Extract UMAP coordinates
umap_coords <- Embeddings(scedata, "umap")  # Extract UMAP coordinates
umap_data <- as.data.frame(umap_coords)
umap_data$celltype <- scedata$celltype  # Add cell type information

# Create a data frame with center point coordinates for each celltype
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# Calculate center point coordinates for each celltype
celltype_centers <- umap_df %>%
  group_by(celltype) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Calculate boundary range for each cell type to determine label placement direction
celltype_ranges <- umap_df %>%
  group_by(celltype) %>%
  summarise(
    min_x = min(UMAP1), max_x = max(UMAP1),
    min_y = min(UMAP2), max_y = max(UMAP2),
    width = max_x - min_x, height = max_y - min_y
  )

# Merge center point coordinates and boundary information
celltype_labels <- left_join(celltype_centers, celltype_ranges, by = "celltype")

# Determine label placement direction: based on cell group shape
celltype_labels <- celltype_labels %>%
  mutate(
    # Determine label placement based on cell group aspect ratio
    direction_x = ifelse(width > height, 1, 0),
    direction_y = ifelse(height > width, 1, 0),
    # Set label nudge parameters: push labels outward from cell groups
    nudge_x = ifelse(UMAP1 > mean(UMAP1), 2, -2),  # Push right-side groups left, left-side right
    nudge_y = ifelse(UMAP2 > mean(UMAP2), 2, -2)   # Push top groups down, bottom groups up
  )

# Plot cell type umap plot and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1) +
  ggrepel::geom_text_repel(
    data = celltype_labels, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 7, 
    fontface = "bold",
    box.padding = 1.5,        # Increase box padding to keep labels away from points
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

# Plot cell type umap plot and save as SVG
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




##### Add cell count for each group #####
# Count cell numbers for each treatment
cell_counts <- scedata@meta.data %>%
  group_by(treatment) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(treatment, " ( n = ", n, " cells)"))

# Build named vector for replacing facet labels
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# Plot PDF
pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),
    width = 6.5*length(unique(scedata$treatment)), height = 5.5)

DimPlot(scedata, reduction = "umap", split.by = "treatment",
        pt.size = 0.1, label = FALSE, cols = col) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # Bold black subplot title
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
    width = 6.5*length(unique(scedata$treatment)), height = 5.5)

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






####### Calculate cell proportion ###########

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1",'#FF9999',"#66CCCC",'#4F6272',"#FF3366","#CC0066","#00CC66","#CC99CC","#FFCCCC","#9999FF","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#FF6699","#6699CC","#FFFFCC")

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
#View(scedata@meta.data)

table(scedata$seurat_clusters)
sum(table(scedata$celltype))

# Set grouping order
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Ovarian cancer","Gastric cancer","Pancreatic cancer","WF-14W"))


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

#cell_counts_group$CellType <- factor(cell_counts_group$CellType, levels = c("3T3", "293T"))

p <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Counts') +
  scale_fill_manual(values = col)+
  #scale_fill_manual(values = c("293T"="#FF9999","3T3"="#A4CDE1"))+
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # Modify X-axis labels
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size
# Add cell count text labels
p <- p + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 6)

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


p <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Ratio') +
  scale_fill_manual(values = col)+
  #scale_fill_manual(values = c("293T"="#FF9999","3T3"="#A4CDE1"))+
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # Modify X-axis labels
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # Modify X-axis text size and rotate 30 degrees
        axis.text.y = element_text(size = 20),  # Modify Y-axis text size
        axis.title.y = element_text(size = 22), # Modify Y-axis title size
        legend.title = element_blank(),         # Remove legend title
        legend.text = element_text(size = 20))  # Modify legend text size
# Add cell proportion text labels
p <- p + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 6)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 6, dpi = 800)



############ Grouping ############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate proportion of each cell type within each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

# Sorting
#cell_counts_treatment$Treatment <- factor(cell_counts_treatment$Treatment, levels = c("Ovarian cancer","Gastric cancer","Pancreatic cancer","WF-14W"))

########## Draw cell count stacked bar plot ##########
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
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)

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
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)







#### Distinguish metastatic and non-metastatic cells ########
#BiocManager::install("UCell")
library(UCell)
library(Seurat)
#install.packages("viridis")
library(viridis)
library(ggplot2)
library(stringr)

# Set output directory
output <- paste(outdir, 'Ucell', sep='/')
dir.create(output, showWarnings = FALSE)

# 1. Read Seurat data
sce <- readRDS("celltype.rds")

# 2. Define gene set
geneSet <- c(
  ## Pancreatic cancer metastasis-related genes
  "SMAD4", "KRAS", "TP53", "CDKN2A", "SNAI1", "SNAI2", "TWIST1", 
  "ZEB1", "ZEB2", "FOSL1", "ITGA2", "ITGB1", "CDH1", "MMP2", "MMP9", 
  "RHOA", "RAC1", "CDC42", "CXCR4", "CXCL12","ACKR3", 
  "TGFB1", "TGFBR1", "TGFBR2", "AXL", "MET", "EGFR", "MYC", "HIF1A", "VEGFA"
)  

# 3. Calculate UCell score
sce <- AddModuleScore_UCell(sce, features = list(Meta = geneSet), name = "_UCell")
score_col  <- colnames(sce@meta.data) %>% str_subset("UCell")


# Extract score vector
scores <- sce@meta.data[[score_col]]

# 4A. Split into 16 equal quantile bins (labels B1..B16)
# Use quantile splitting to avoid empty bins due to uneven value distribution (remove duplicate quantile points)
probs <- seq(0, 1, length.out = 17)
brks_raw <- quantile(scores, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
brks <- unique(brks_raw)

# If deduplication results in fewer than 16 bins due to many tied values, fall back to pretty splitting to ensure bin count
if (length(brks) < 17) {
  brks <- pretty(range(scores, na.rm = TRUE), n = 16)
  brks <- unique(brks)
}

# Ensure at least 2 bins (fallback for extreme cases)
if (length(brks) < 2) {
  brks <- range(scores, na.rm = TRUE)
  brks[1] <- brks[1] - 1e-8
  brks[2] <- brks[2] + 1e-8
}

# Generate 16-bin (if effective breakpoints are less than 17, bin count will be slightly less than 16, due to data itself)
n_bins <- length(brks) - 1
bin_labels <- paste0("B", seq_len(n_bins))
sce$Meta_bin16 <- cut(scores, breaks = brks, include.lowest = TRUE, labels = bin_labels, right = TRUE, ordered_result = TRUE)


# 4B. 20/60/20 classification (High / Medium / Low)
q20 <- quantile(scores, 0.20, na.rm = TRUE, type = 7)
q80 <- quantile(scores, 0.80, na.rm = TRUE, type = 7)
sce$Meta_3class <- factor(
  dplyr::case_when(
    scores >= q80 ~ "Met_high",
    scores <= q20 ~ "Met_low",
    TRUE ~ "Met_medium"
  ),
  levels = c("Met_high", "Met_medium", "Met_low")
)

saveRDS(sce,  "celltype.rds")
View(sce@meta.data)

# 5. Visualization
# If UMAP/TSNE is not available, first run: sce <- RunPCA(sce); sce <- RunUMAP(sce, dims = 1:30)
# 5A. 16-bin discrete visualization (DimPlot more suitable for discrete data)
p_bin16 <- DimPlot(
  sce,
  reduction = if ("umap" %in% Reductions(sce)) "umap" else DefaultDimReduc(sce),
  group.by = "Meta_bin16",
  label = FALSE,
  cols = viridis(n = max(n_bins, 3), option = "C")
) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(title = "Metastatic cell score (16-quantile bins)", color = "Bin")

ggsave(file.path(output, "Meta_16bin_DimPlot.pdf"), plot = p_bin16,
       width = 6, height = 4, device = cairo_pdf)


# 5B. High/Medium/Low visualization
# Specify clear discrete color palette for three categories
#hml_cols <- setNames(viridis(3, option = "D"), c("Low","Medium","High"))
p_hml <- DimPlot(
  sce,
  reduction = if ("umap" %in% Reductions(sce)) "umap" else DefaultDimReduc(sce),
  group.by = "Meta_3class",
  label = FALSE,
  cols = c("#FF0066","#CCFFFF","#0066CC" )
) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.position = "right"
  ) +
  labs(title = "Metastatic cell score", color = "Class")


ggsave(filename = file.path(output, "Meta_class_DimPlot.pdf"),
       plot = p_hml, width = 5, height = 4)

# 7. Optional: keep continuous score FeaturePlot (original figure)
p_cont <- FeaturePlot(
  sce,
  features = score_col,
  order = TRUE,
  ncol = 1,
  cols = c("#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
           "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36","#F5E24B")
) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(title = "Metastatic cell score", color = "Met_score")

ggsave(filename = file.path(output, "Meta_FeaturePlot.pdf"),
       plot = p_cont, width = 4, height = 3)


### Grouping ###
# 7. Optional: keep continuous score FeaturePlot (original figure)
p_cont <- FeaturePlot(
  sce,
  features = score_col,
  split.by = 'treatment',
  order = TRUE,
  ncol = 2,
  cols = c(  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
             "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36","#F5E24B")
) +
  theme(
    strip.text   = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(color = "Met_score")

ggsave(filename = file.path(output, "Meta_FeaturePlot_treatment.pdf"),
       plot = p_cont, width = 8, height = 3)









#BiocManager::install("impute", ask = FALSE, update = FALSE)
#devtools::install_github("Japrin/sscVis") #STARTRAC tool installation
#devtools::install_github("Japrin/Startrac") #STARTRAC tool installation

library(Startrac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)


sco <- readRDS("celltype.rds")
data <- sco@meta.data
colnames(data)


#Construct input table needed for Roe calculation
#data <- data[,c(1,4,7,11)]
#colnames(data) <- c("sample","tissue","celltype")


# ---- Calculate Roe and prepare matrix ----
Roe <- calTissueDist(
  data,
  byPatient = FALSE,
  colname.cluster = "Meta_3class",
  colname.patient = "orig.ident",
  colname.tissue = "treatment",
  method = "chisq",   # "chisq", "fisher", "freq"
  min.rowSum = 0
)

mat <- as.matrix(Roe)  # Use matrix consistently to avoid pitfalls with data frame min/max
rng <- range(mat, na.rm = TRUE)

# If you want 1 as the color midpoint (typical Ro/e visual semantics), but 1 is not in the data range, use the range midpoint instead
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng)

# Color mapping: low value -> light color, midpoint (1 or range midpoint) -> orange, high value -> red
col_fun <- circlize::colorRamp2(
  c(rng[1], mid, rng[2]),
  c("#f6f8e6", "#f9a33e", "red")
)

# Legend ticks: use pretty() to generate "nice" ticks; labels correspond to at one-to-one
legend_at <- pretty(rng, n = 5)
legend_labels <- formatC(legend_at, format = "f", digits = 2)

# ---- Figure 1: Heatmap with values ----
pdf(file.path(output, "STARTRAC_Roe_value.pdf"), width = 6, height = 3)
Heatmap(
  mat,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "right",
  column_names_side = "bottom",
  show_column_names = TRUE,
  show_row_names = TRUE,
  col = col_fun,
  row_names_gp = gpar(
    fontsize = 16,
    col = col,  # Specify different color vector for row names
    fontface = "bold"      # Optional, bold font for enhanced display
  ),
  column_names_gp = gpar(fontsize = 16),
  heatmap_legend_param = list(
    title = "Ro/e value",
    at = legend_at,
    labels = legend_labels,
    title_gp = gpar(fontsize = 16),  # Increase legend title font size
    labels_gp = gpar(fontsize = 15)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 16, col = "black"))
  }
)
dev.off()




#Visualization 2, custom +++ symbol version
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0


#Visualization 3, custom +++ symbol and row cell type color version
pdf(file.path(output, "STARTRAC_Roe_colored.pdf"), width = 6, height = 3)
Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = 'right', 
        column_names_side = "bottom",
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        column_names_gp = gpar(fontsize = 16),
        row_names_gp = gpar(
          fontsize = 16,
          col = col,  # Specify different color vector for row names
          fontface = "bold"      # Optional, bold font for enhanced display
        ),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(0, max(Roe)), 
          title_gp = gpar(fontsize = 16),  # Increase legend title font size
          labels_gp = gpar(fontsize = 15),
          labels = c("0", "Max.")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          value <- Roe[i, j]
          symbol <- if(value == 0) {
            "−"
          } else if(value > 0 & value < 0.2) {
            "+/−"
          } else if(value >= 0.2 & value <= 0.8) {
            "+"
          } else if(value > 0.8 & value <= 1) {
            "++"
          } else if(value > 1) {
            "+++"
          }
          # Optimize text contrast (automatically switch based on background color)
          text_color <- ifelse(mean(col2rgb(fill)) > 127, "black", "white")
          grid.text(symbol, x, y, gp = gpar(fontsize = 16, col = text_color))
        }
)
dev.off()




#Visualization 4, bubble plot
roe_df <- as.data.frame(as.table(as.matrix(Roe))) %>% #Convert matrix to data frame
  rename(Tissue = Var2, CellType = Var1, Roe = Freq) %>% #Modify new data frame column names
  mutate(
    Enrichment = ifelse(Roe >= 1, "Enrichment", "Depletion"),  #Determine enrichment (greater than 1) or depletion (less than 1) based on Roe value
    Roe = ifelse(Roe == 0, NA, Roe)  #Handle 0 values
  ) %>% 
  filter(!is.na(Roe))  #Filter invalid values

pdf(file.path(output, "STARTRAC_Roe_bubble.pdf"), width = 5, height = 3)
ggplot(roe_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(
    aes(size = Roe, color = Enrichment)
  ) +
  scale_size_continuous(
    name = "Ro/e",
    breaks = c(0.5, 1.0, 1.5),
    range = c(1,7) #Adjust bubble size range
  ) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),  # Example plot colors
    labels = c("Enrichment", "Depletion"),
    guide = guide_legend(  
      override.aes = list(
        size = 4  #Adjust point size in legend
      )
    )
  ) +
  scale_y_discrete(limits = rev) +  # Keep y-axis order consistent with input
  labs(
    title = "STARTRAC - Cell Type Enrichment",
    x = "Tissue Group",
    y = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.box = "horizontal" #Relative position of different legends, horizontal or vertical
  )
dev.off()









#####19. Differential analysis comparing two groups #####
library(scRNAtoolVis)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggrepel)

col <- c('#437eb8','#FF6666',"#FFFFCC",'#FFCC99','#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFCCCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")

# Display expression trends of marker genes for top 6 clusters
output <- paste(outdir,'Differential_analysis(metastatic_vs_non_metastatic)', sep='/')
dir.create(output, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)
colnames(scRNAsub@meta.data)


# Find differentially expressed genes between Res and Sen groups
logFCfilter <- 0.25        # Define log2FC filter value
adjPvalFilter <- 0.05   # Define adjusted P-value filter value

# Find differentially expressed genes between Epi_cisplatin_res and Epi_other groups
scRNAsub.cluster.markers <- FindMarkers(object = scRNAsub, 
                                        ident.1 = "Met_high",
                                        ident.2 =  "Met_low",
                                        group.by = "Meta_3class", 
                                        logfc.threshold = 0, 
                                        min.pct = 0.25, 
                                        test.use = "wilcox")
scRNAsub.cluster.markers$gene <- rownames(scRNAsub.cluster.markers)

# Add significance annotation
scRNAsub.cluster.markers <- scRNAsub.cluster.markers %>%
  mutate(Significance = ifelse(p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter, 
                               ifelse(avg_log2FC > 0, "Up", "Down"), "Normal"))
write.table(scRNAsub.cluster.markers, file = file.path(output,"sig.markers_ann_Tumor_vs_Normal.txt"), sep = "\t",row.names = T, quote = FALSE)

saveRDS(scRNAsub.cluster.markers, file = file.path(output, "ScRNA.sig.markers.rds"))

# Save upregulated and downregulated genes separately
upregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")
downregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")
write.csv(upregulated_genes, file = file.path(output, "upregulated_genes_Tumor_vs_Normal.csv"), row.names = TRUE, quote = FALSE)
write.csv(downregulated_genes, file = file.path(output, "downregulated_genes_Tumor_vs_Normal.csv"), row.names = TRUE, quote = FALSE)

# Calculate number of upregulated and downregulated genes
upregulated_genes <- sum(scRNAsub.cluster.markers$Significance == "Up")
downregulated_genes <- sum(scRNAsub.cluster.markers$Significance == "Down")
total_diff_genes <- upregulated_genes + downregulated_genes

# Save data frames of upregulated and downregulated genes separately
upregulated_genes_df <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")
downregulated_genes_df <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")

# Select top 10 upregulated and downregulated genes to display labels
top_genes_upregulated <- upregulated_genes_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  head(15)
top_genes_downregulated <- downregulated_genes_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  head(15)


# Genes of interest
genes <- c(
  
  "S100A4","NEDD9","FN1","CD44","MMP19","LOXL2","LOX","SERPINE1",
  #"VIM","COL1A1","VCAN","THBS1","SPARC","RHOB","RND3",
  "TIMP2","CXCL1","CXCL2","CCL2","CCL3","PDGFRA","VEGFA","S100A10",
  "JUN","FOS","FOSL1","NFKB1","STAT3"
  
)

# Filter genes of interest
interested_genes <- scRNAsub.cluster.markers %>%
  filter(gene %in% genes)

# Draw volcano plot
p <- ggplot(scRNAsub.cluster.markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), size = 2, shape = 18) +
  scale_color_manual(values = c("#339999", "#FFCCCC", "#FF0066")) +
  geom_hline(yintercept = -log10(adjPvalFilter), linetype = "dashed") +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed") +
  #geom_text_repel(data = top_genes_upregulated, aes(label = top_genes_upregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  #geom_text_repel(data = top_genes_downregulated, aes(label = top_genes_downregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  geom_text_repel(data = interested_genes, aes(label = interested_genes$gene), size = 5, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  theme_classic() +
  labs(title = "Met_high vs Met_low", 
       x = "log2 Fold Change", y = "-log10 Adjusted P-value", color = "Significance") +
  #annotate("text", x = -1.5, y = 300, label = paste("Up-regulated:", upregulated_genes), 
  #         hjust = 0, size = 5, color = "#FF3300", fontface = "bold") +
  #annotate("text", x = -1.5, y = 285, label = paste("Down-regulated:", downregulated_genes), 
  #         hjust = 0, size = 5, color = "#6699CC", fontface = "bold") +
  #annotate("text", x = -1.5, y = 270, label = paste("DEGs:", total_diff_genes), 
  #         hjust = 0, size = 5, color = "black", fontface = "bold") +
  scale_x_continuous(
    limits = c(-2, 3),                      # Set X-axis range
    breaks = seq(-1.5, 3, by = 1),            # Set ticks
    expand = expansion(mult = c(0.05, 0.05)) # Increase middle area expansion
  ) +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0), 
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 18))

# Save image
ggsave(file.path(output, "Tumor_vs_Normal_volcano_plot.svg"), p, width = 8, height = 7, dpi = 300)
ggsave(file.path(output, "Tumor_vs_Normal_volcano_plot.pdf"), p, width = 8, height = 7, dpi = 300)



library(ggrepel)
library(ggplot2)

# Add difference calculation column
scRNAsub.cluster.markers <- scRNAsub.cluster.markers %>%
  mutate(Difference = pct.1 - pct.2)


# Draw volcano plot
volcano_plot <- ggplot(scRNAsub.cluster.markers, aes(x = Difference, y = avg_log2FC, color = Significance)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("blue", "grey", "red")) + 
  # Annotate genes of interest
  geom_label_repel(data = subset(scRNAsub.cluster.markers, gene %in% genes), 
                   aes(label = gene), 
                   color = "black", 
                   segment.colour = "black", 
                   label.padding = 0.1, 
                   segment.size = 0.3, 
                   size = 6,
                   max.overlaps = 50) + 
  geom_vline(xintercept = 0.0, linetype = 2) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),  # X-axis title font size bold
    axis.title.y = element_text(size = 20, face = "bold"),  # Y-axis title font size bold
    axis.text.x = element_text(size = 18),  # X-axis text font size
    axis.text.y = element_text(size = 18),  # Y-axis text font size
    legend.title = element_text(size = 20, face = "bold"),  # Legend title font size bold
    legend.text = element_text(size = 18),  # Legend text font size
    plot.title = element_text(size = 22, face = "bold"),  # Plot title font size bold
    plot.subtitle = element_text(size = 20),  # Subtitle font size
    legend.position = "right",  # Legend position
    legend.key.size = unit(1.5, "lines")  # Legend key size
  ) +
  xlim(-0.35, 0.6) +
  ylim(-3, 2) +
  labs(title = "Tumor vs Normal")  # Set title and subtitle

# Save as PDF and SVG
ggsave(file.path(output, "volcano_plot.pdf"), volcano_plot, width = 8, height = 6)
ggsave(file.path(output, "volcano_plot.svg"), volcano_plot, width = 8, height = 6)





############### Expression of differential genes in different groups ##################
library(Seurat)
library(tidyverse)
library(ggsci)

ScRNA <- scRNAsub
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
write.table(avg_expr_df, file = paste0(output, "/Differential_gene_expression.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)




# ✅ Draw DotPlot
# ✅ Extract genes to plot (confirm again that genes are in the object)
genes_to_plot <- genes[genes %in% rownames(ScRNA)]

# Calculate normalized expression of Tumor group relative to other groups
# First get average expression for all groups
avg_exp <- Seurat::AverageExpression(ScRNA, 
                                     features = genes_to_plot,
                                     group.by = "treatment",
                                     assays = "RNA")$RNA

# Calculate normalized expression of Tumor group relative to other groups
# Here we calculate Tumor group expression minus the average expression of other groups
other_groups <- setdiff(colnames(avg_exp), "Tumor")
tumor_normalized <- avg_exp[,"Tumor"] - rowMeans(avg_exp[,other_groups, drop = FALSE])

# Sort genes based on normalized expression
# First group by positive/negative values, then sort by absolute value
genes_positive <- names(sort(tumor_normalized[tumor_normalized > 0], decreasing = TRUE))
genes_negative <- names(sort(tumor_normalized[tumor_normalized <= 0], decreasing = FALSE))

# Merge gene order
genes_ordered <- c(genes_positive, genes_negative)

# Ensure all genes to plot are in the sorted list
genes_to_plot <- intersect(genes_ordered, genes_to_plot)



# Draw gene expression dot plot (grouped by treatment)
plot <- DotPlot(ScRNA, features = genes_to_plot, group.by = "treatment") + 
  RotatedAxis() +
  coord_flip() +
  scale_color_gradientn(colors = c('#CCCCCC', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 22), 
        axis.title.x = element_text(size = 22), 
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20))

# Save DotPlot figure
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.pdf", sep='/'), plot = plot, width = 7, height = 12)
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.svg", sep='/'), plot = plot, width = 7, height = 12)






##### Perform GSEA enrichment analysis after differential gene analysis #####
#library(org.Mm.eg.db) # Mouse database
library(org.Hs.eg.db) # Human database
library(clusterProfiler)
library(enrichplot)
library(DOSE)

scRNAsub.cluster.markers <- readRDS(file.path(output, "ScRNA.sig.markers.rds"))

# Get differential gene list from differential analysis
deg <- scRNAsub.cluster.markers[, c('avg_log2FC', 'p_val_adj')]
colnames(deg) <- c('log2FoldChange', 'pvalue')  # Change column names

# Convert SYMBOL to ENTREZID
gene <- bitr(rownames(deg), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Match logFC information
gene$logFC <- deg$log2FoldChange[match(gene$SYMBOL, rownames(deg))]

# Build genelist
geneList <- gene$logFC
names(geneList) <- gene$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)  # Sort by logFC descending

# GSEA analysis (KEGG pathways)
kk_gse <- gseKEGG(geneList = geneList,
                  organism =  "hsa",  
                  ## 'mmu',  # Mouse
                  nPerm = 1000,       ## Number of permutations: randomly shuffle gene sets 1000 times for simulation
                  minGSSize = 10,     ## Gene set for enrichment analysis must contain at least 10 genes
                  pvalueCutoff = 0.25,
                  verbose = FALSE)

# Convert ENTREZID to readable SYMBOL names
kk_gse <- DOSE::setReadable(kk_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(as.data.frame(kk_gse), file = file.path(output, "kk_gse_results.csv"))

# Filter significantly enriched pathways |NES| > 1, p-value < 0.05, FDR < 0.25
kk_gse_cut <- kk_gse[kk_gse$pvalue < 0.05 & kk_gse$p.adjust < 0.25 & abs(kk_gse$NES) > 1, ]

# Upregulated pathways
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0, ]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES, decreasing = TRUE), 10), ]

# Downregulated pathways
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0, ]
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES, decreasing = TRUE), 10), ]

# Draw GSEA plot for upregulated pathways
gseap_up <- gseaplot2(kk_gse,
                      up_gsea$ID,
                      title = up_gsea$Description[1], # Use description of first upregulated pathway as title
                      color = c("#FF4500", "#32CD32"),
                      base_size = 30, # Base font size
                      rel_heights = c(1.5, 0.5, 1), # Relative heights of subplots
                      subplots = 1:3, # Display subplots
                      ES_geom = "line", # Enrichment score line style
                      pvalue_table = F)  # Display p-value table 
gseap_up[[1]]<-gseap_up[[1]]+
  scale_color_viridis_d()+
  theme(plot.title = element_text(size = 28, face = "bold"),  # Increase title font size
        legend.title = element_text(size = 20),               # Increase legend title font size
        legend.text = element_text(size = 20))
gseap_up[[2]]<-gseap_up[[2]]+
  scale_color_viridis_d()

ggsave(file.path(output, "GSEA_up.pdf"), gseap_up, width = 15, height = 12)
ggsave(file.path(output, "GSEA_up.svg"), gseap_up, width = 15, height = 12)

# Draw GSEA plot for downregulated pathways
gseap_down <- gseaplot2(kk_gse,
                        down_gsea$ID,
                        title = "DOWN_GSEA", # Title
                        color = c("#FF4500", "#32CD32"),
                        base_size = 25, # Base font size
                        rel_heights = c(1.5, 0.5, 1), # Relative heights of subplots
                        subplots = 1:3, # Display subplots
                        ES_geom = "line", # Enrichment score line style
                        pvalue_table = F) # Display p-value table
gseap_down[[1]]<-gseap_down[[1]]+
  scale_color_viridis_d()+
  theme(plot.title = element_text(size = 28, face = "bold"),  # Increase title font size
        legend.title = element_text(size = 20),               # Increase legend title font size
        legend.text = element_text(size = 16))
gseap_down[[2]]<-gseap_down[[2]]+
  scale_color_viridis_d()
ggsave(file.path(output, "GSEA_down.pdf"), gseap_down, width = 18, height = 12)
ggsave(file.path(output, "GSEA_down.svg"), gseap_down, width = 18, height = 12)

# Ridgeplot visualization, showing top 15 pathways
ridgep <- ridgeplot(kk_gse, 
                    showCategory = 15, 
                    fill = "pvalue",  # Adjust fill color based on p-value
                    core_enrichment = TRUE,
                    label_format = 30,  # Set character length for axis labels, wrap if too long
                    orderBy = "NES", 
                    decreasing = FALSE)+
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),  
        legend.text = element_text(size = 12)) 

ggsave(file.path(output, "ridgeplot_GSEA.pdf"), ridgep, width = 10, height = 8)
ggsave(file.path(output, "ridgeplot_GSEA.svg"), ridgep, width = 10, height = 8)


####### GO and KEGG enrichment analysis
# Filter upregulated and downregulated genes
gene_up <- scRNAsub.cluster.markers$gene[scRNAsub.cluster.markers$Significance == "Up"]
gene_down <- scRNAsub.cluster.markers$gene[scRNAsub.cluster.markers$Significance == "Down"]

# Convert SYMBOL to ENTREZID
gene_up_entrez <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, 
                                                             keys = gene_up, 
                                                             columns = 'ENTREZID', 
                                                             keytype = 'SYMBOL')[,2]))
gene_down_entrez <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, 
                                                               keys = gene_down, 
                                                               columns = 'ENTREZID', 
                                                               keytype = 'SYMBOL')[,2]))

# Perform GO enrichment analysis
go_up <- enrichGO(gene = gene_up_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)
go_down <- enrichGO(gene = gene_down_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)

# Convert geneID from ENTREZID to SYMBOL
go_up@result$geneID <- sapply(strsplit(go_up@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
  paste(symbols, collapse = "/")
})

go_down@result$geneID <- sapply(strsplit(go_down@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
  paste(symbols, collapse = "/")
})

write.csv(as.data.frame(go_up), file = file.path(output, "go_up_results.csv"))
write.csv(as.data.frame(go_down), file = file.path(output, "go_down_results.csv"))

# Draw GO enrichment analysis plot
go_plot_up <- dotplot(go_up) + ggtitle("Upregulated Genes GO Enrichment")
go_plot_down <- dotplot(go_down) + ggtitle("Downregulated Genes GO Enrichment")

# Save GO analysis results
ggsave(file.path(output, 'go_enrich_up_dot.pdf'), plot = go_plot_up, width = 6, height = 6)
ggsave(file.path(output, 'go_enrich_up_dot.svg'), plot = go_plot_up, width = 6, height = 6)

ggsave(file.path(output, 'go_enrich_down_dot.pdf'), plot = go_plot_down, width = 6, height = 5)
ggsave(file.path(output, 'go_enrich_down_dot.svg'), plot = go_plot_down, width = 6, height = 5)

# Perform KEGG enrichment analysis
kegg_up <- enrichKEGG(gene = gene_up_entrez, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.1)
kegg_down <- enrichKEGG(gene = gene_down_entrez, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.1)

# Convert geneID from ENTREZID to SYMBOL
kegg_up@result$geneID <- sapply(strsplit(kegg_up@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
  paste(symbols, collapse = "/")
})

kegg_down@result$geneID <- sapply(strsplit(kegg_down@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
  paste(symbols, collapse = "/")
})

write.csv(as.data.frame(kegg_up), file = file.path(output, "kegg_up_results.csv"))
write.csv(as.data.frame(kegg_down), file = file.path(output, "kegg_down_results.csv"))

# Draw KEGG enrichment analysis plot
kegg_plot_up <- dotplot(kegg_up, showCategory = 20)   + ggtitle("Upregulated Genes KEGG Enrichment")
kegg_plot_down <- dotplot(kegg_down, showCategory = 20) + ggtitle("Downregulated Genes KEGG Enrichment")

# Save KEGG analysis results
ggsave(file.path(output, 'kegg_enrich_up_dot.pdf'), plot = kegg_plot_up, width = 7, height = 8)
ggsave(file.path(output, 'kegg_enrich_up_dot.svg'), plot = kegg_plot_up, width = 7, height = 8)

ggsave(file.path(output, 'kegg_enrich_down_dot.pdf'), plot = kegg_plot_down, width = 8, height = 8)
ggsave(file.path(output, 'kegg_enrich_down_dot.svg'), plot = kegg_plot_down, width = 8, height = 8)


# Integrate and display GO and KEGG analysis plots
combined_plot_GO <- go_plot_up + go_plot_down + plot_layout(guides = 'collect')
combined_plot_KEGG <-  kegg_plot_up + kegg_plot_down + plot_layout(guides = 'collect')

# Save integrated images
ggsave(file.path(output, 'combined_GO_dot.pdf'), plot = combined_plot_GO, width = 13, height =10)
ggsave(file.path(output, 'combined_KEGG_dot.pdf'), plot = combined_plot_KEGG, width = 13, height = 12)



# Prepare visualization data
# Extract GO analysis results
go_up_dt <- as.data.frame(go_up)
go_down_dt <- as.data.frame(go_down)

# Extract GO pathways of interest
interested_terms <- c(
  "Focal adhesion", "Ubiquitin mediated proteolysis", 
  "Regulation of actin cytoskeleton", "Hedgehog signaling pathway", 
  "Wnt signaling pathway", "Phospholipase D signaling pathway", 
  "Autophagy - animal", "ECM-receptor interaction", "Sphingolipid signaling pathway"
)

# Extract GO pathways of interest
interested_terms <- c(
  "small GTPase-mediated signal transduction",
  "protein localization to plasma membrane", "Wnt signaling pathway", "Rho protein signal transduction",
  "cellular response to epidermal growth factor stimulus", "regulation of autophagy", 
  "response to epidermal growth factor", 
  "ERK1 and ERK2 cascade", "positive regulation of MAPK cascade", 
  "positive regulation of cell projection organization", 
  "substrate-dependent cell migration", "receptor-mediated endocytosis", 
  "epidermal growth factor receptor signaling pathway", "positive regulation of autophagy", "Notch signaling pathway", 
  "intracellular receptor signaling pathway", "ERBB signaling pathway", "cellular response to fibroblast growth factor stimulus"
)



# Filter pathways of interest from GO enrichment results
go_up_dt <- go_up_dt[go_up_dt$Description %in% interested_terms, ]



# Extract KEGG analysis results
kegg_up_dt <- as.data.frame(kegg_up)
kegg_down_dt <- as.data.frame(kegg_down)

# Extract GO pathways of interest
interested_terms <- c(
  "Focal adhesion","ECM-receptor interaction","Regulation of actin cytoskeleton",
  "Rap1 signaling pathway","TGF-β signaling pathway","PI3K-Akt signaling pathway",
  "MAPK signaling pathway","NF-κB signaling pathway","HIF-1 signaling pathway",
  "Chemokine signaling pathway","IL-17 signaling pathway","Ras signaling pathway",
  "ErbB signaling pathway","EGFR-TKI resistance","Hippo signaling pathway",
  "PD-L1 expression and PD-1 checkpoint pathway in cancer","Platelet activation",
  "Complement and coagulation cascades","Leukocyte transendothelial migration",
  "Tight junction","Adherens junction","Endocytosis","Lysosome","Autophagy-animal",
  "Cytokine–cytokine receptor interaction","FoxO signaling pathway"
)

# Filter pathways of interest from GO enrichment results
kegg_up_dt <- kegg_up_dt[kegg_up_dt$Description %in% interested_terms, ]

# Set color classification
classification_colors <- c('#437eb8','#FF6666','#FFCC99','#FF9999', '#80c5d8',"#9999FF",
                           "#FFCCCC","#99CCFF","#FF3366","#CCCCFF","#CC0066","#FFFFCC",
                           "#66CCCC","#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                           "#6699CC","#CC99CC","#FF6699","#FF0000","#6666CC","#FF9966",
                           "#669999","#CC99FF","#FFCCFF")

classification_colors<- c(
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


# Text wrapping function, limit character length
wrap_text <- function(text, width = 40) {
  sapply(text, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

# Bar plot for GO enrichment analysis
plot_GO_bar <- function(dt, title) {
  dt <- dt[order(-dt$pvalue, decreasing = TRUE),]  # Sort by p-value first
  #dt <- dt[order(dt$Count, decreasing = TRUE),]  # Sort by gene count first
  
  dt <- head(dt,20)  # Select top 20 pathways
  dt$Description <- factor(wrap_text(dt$Description), levels = wrap_text(dt$Description))
  
  # Left plot: enrichment p-value
  p1 <- ggplot(dt, aes(x = Description, y = log10(p.adjust), fill = Description)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = classification_colors) +
    coord_flip() +
    ylab('-log10(P-value)') +
    xlab('') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 16, face = "bold"),
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
          axis.text.x = element_text(size = 16),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(size = 16, face = "bold"),
          plot.margin = margin(10, 10, 10, 10),
          panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add XY axis borders
  
  # Combine two plots
  p_combined <- p1 + p2 + plot_layout(widths = c(2, 1.5))
  return(p_combined)
}

# Draw GO and KEGG enrichment analysis bar plots
go_up_plot <- plot_GO_bar(go_up_dt, "Upregulated Genes GO Enrichment")
go_down_plot <- plot_GO_bar(go_down_dt, "Downregulated Genes GO Enrichment")

kegg_up_plot <- plot_GO_bar(kegg_up_dt, "Upregulated Genes KEGG Enrichment")
kegg_down_plot <- plot_GO_bar(kegg_down_dt, "Downregulated Genes KEGG Enrichment")

# Save images
ggsave(file.path(output, 'go_enrich_up_bar.pdf'), plot = go_up_plot, width = 8, height = 6)
#ggsave(file.path(output, 'go_enrich_down_bar.pdf'), plot = go_down_plot, width = 8, height = 6)
ggsave(file.path(output, 'kegg_enrich_up_bar.pdf'), plot = kegg_up_plot, width = 10, height = 8)
#ggsave(file.path(output, 'kegg_enrich_down_bar.pdf'), plot = kegg_down_plot, width = 8, height = 5)

# Combine and save
#combined_GO_plot <- go_up_plot + go_down_plot + plot_layout(guides = 'collect')
#combined_KEGG_plot <- kegg_up_plot + kegg_down_plot + plot_layout(guides = 'collect')

#ggsave(file.path(output, 'combined_GO_bar.pdf'), plot = combined_GO_plot, width = 13, height = 10)
#ggsave(file.path(output, 'combined_KEGG_bar.pdf'), plot = combined_KEGG_plot, width = 13, height = 10)












#### Batch calculate UCell for multiple gene sets and visualize ####
#BiocManager::install("UCell")
library(UCell)
library(Seurat)
library(viridis)
library(ggplot2)
library(stringr)
library(dplyr)

# =============== Basic settings ===============
if (!exists("outdir")) outdir <- "."  # If outdir is not defined, use current directory
output_root <- file.path(outdir, "Ucell")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

# 1. Read Seurat object
sce <- readRDS("celltype.rds")

# 2. Define multiple gene sets of interest (named list)
gene_sets <- list(
  
  Cancer = c("EPCAM","WFDC2","KRT8","CLDN4","CD24"),
  
  Ovarian = c("WT1","MUC16","PAX8"),
  
  Digestive_tract = c("KRT7", "KRT20","CDX2"),
  
  Metastatic = c(  "SMAD4", "KRAS", "TP53", "CDKN2A", "SNAI1", "SNAI2", "TWIST1", 
                   "ZEB1", "ZEB2", "FOSL1", "ITGA2", "ITGB1", "CDH1", "MMP2", "MMP9", 
                   "RHOA", "RAC1", "CDC42", "CXCR4", "CXCL12","ACKR3", 
                   "TGFB1", "TGFBR1", "TGFBR2", "AXL", "MET", "EGFR", "MYC", "HIF1A", "VEGFA")
  
  
)

# =============== Utility functions ===============
# Generate bin breaks (prefer quantiles, fallback to pretty if too many duplicates)
.make_breaks <- function(x, n_bins_target = 16) {
  probs <- seq(0, 1, length.out = n_bins_target + 1)
  brks_raw <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
  brks <- unique(brks_raw)
  if (length(brks) < (n_bins_target + 1)) {
    brks <- pretty(range(x, na.rm = TRUE), n = n_bins_target)
    brks <- unique(brks)
  }
  if (length(brks) < 2) {
    brks <- range(x, na.rm = TRUE)
    brks[1] <- brks[1] - 1e-8
    brks[2] <- brks[2] + 1e-8
  }
  return(brks)
}

# Select dimensionality reduction (prefer UMAP, otherwise use DefaultDimReduc)
.pick_reduction <- function(obj) {
  if ("umap" %in% Reductions(obj)) return("umap")
  DefaultDimReduc(obj)
}

# =============== Main loop: process each gene set ===============
for (set_name in names(gene_sets)) {
  message("Processing gene set: ", set_name)
  
  # Create output directory for this gene set
  out_dir <- file.path(output_root, set_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Gene availability check (only use genes present in expression matrix)
  genes_in <- intersect(gene_sets[[set_name]], rownames(sce))
  genes_miss <- setdiff(gene_sets[[set_name]], genes_in)
  if (length(genes_miss) > 0) {
    message(sprintf("[%s] %d gene(s) not found and will be skipped: %s",
                    set_name, length(genes_miss), paste(genes_miss, collapse = ", ")))
  }
  if (length(genes_in) < 2) {
    warning(sprintf("[%s] Fewer than 2 genes present in the object; skipping this set.", set_name))
    next
  }
  
  # 3. Calculate UCell score (column name will be '<set_name>_UCell')
  sce <- AddModuleScore_UCell(
    sce,
    features = setNames(list(genes_in), set_name),  # Name elements to ensure column name prefix
    name = "_UCell"
  )
  score_col <- paste0(set_name, "_UCell")
  if (!score_col %in% colnames(sce@meta.data)) {
    stop(sprintf("[%s] UCell score column not found: %s", set_name, score_col))
  }
  scores <- sce@meta.data[[score_col]]
  
  # 4A. 16 equal quantile bins
  brks <- .make_breaks(scores, n_bins_target = 16)
  n_bins <- length(brks) - 1
  bin_labels <- paste0("B", seq_len(n_bins))
  bin_col <- paste0(set_name, "_bin16")
  sce[[bin_col]] <- cut(
    scores, breaks = brks, include.lowest = TRUE, labels = bin_labels,
    right = TRUE, ordered_result = TRUE
  )
  
  # 4B. 20/60/20 classification (high/medium/low)
  q20 <- quantile(scores, 0.20, na.rm = TRUE, type = 7)
  q80 <- quantile(scores, 0.80, na.rm = TRUE, type = 7)
  class_col <- paste0(set_name, "_3class")
  sce[[class_col]] <- factor(
    dplyr::case_when(
      scores >= q80 ~ paste0(set_name, "_high"),
      scores <= q20 ~ paste0(set_name, "_low"),
      TRUE          ~ paste0(set_name, "_medium")
    ),
    levels = c(paste0(set_name, "_high"),
               paste0(set_name, "_medium"),
               paste0(set_name, "_low"))
  )
  
  # 5. Visualization
  red_use <- .pick_reduction(sce)
  
  ## 5A. 16-bin discrete visualization (DimPlot)
  p_bin16 <- DimPlot(
    sce,
    reduction = red_use,
    group.by = bin_col,
    label = FALSE,
    cols = viridis(n = max(n_bins, 3), option = "C")
  ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(title = sprintf("%s score", set_name), color = "Bin")
  
  ggsave(file.path(out_dir, paste0(set_name, "_16bin_DimPlot.pdf")),
         plot = p_bin16, width = 6, height = 4, device = cairo_pdf)
  
  ## 5B. High/Medium/Low visualization
  p_hml <- DimPlot(
    sce,
    reduction = red_use,
    group.by = class_col,
    label = FALSE,
    cols = c("#FF0066", "#CCFFFF", "#0066CC")
  ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.position = "right"
    ) +
    labs(title = sprintf("%s score", set_name), color = "Class")
  
  ggsave(file.path(out_dir, paste0(set_name, "_class_DimPlot.pdf")),
         plot = p_hml, width = 7, height = 4)
  
  ## 5C. Continuous score FeaturePlot (overall)
  p_cont <- FeaturePlot(
    sce,
    features = score_col,
    order = TRUE,
    ncol = 1,
    cols = c("#CCCCCC", "#0066CC")
  ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12),
      legend.position = "right"
    ) +
    labs(title = sprintf("%s score", set_name), color = "Score")
  
  ggsave(file.path(out_dir, paste0(set_name, "_FeaturePlot.pdf")),
         plot = p_cont, width = 4, height = 3)
  
  ## 5D. Continuous score FeaturePlot (grouped by treatment)
  # If 'treatment' column is not present, skip grouping plot
  if ("treatment" %in% colnames(sce@meta.data)) {
    p_cont_split <- FeaturePlot(
      sce,
      features = score_col,
      split.by = "treatment",
      order = TRUE,
      ncol = 2,
      cols = c("lightgrey", "#0066CC")
    ) +
      theme(
        strip.text   = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        legend.position = "right"
      ) +
      labs(color = "Score")
    
    ggsave(file.path(out_dir, paste0(set_name, "_FeaturePlot_treatment.pdf")),
           plot = p_cont_split, width = 8, height = 3)
  } else {
    message(sprintf("[%s] 'treatment' column not found in meta.data; skip split plot.", set_name))
  }
  
  
  
  ## 5E. New: Draw violin plot of UCell score grouped by celltype
  # One plot per gene set, x = celltype, y = UCell score for that gene set
  message(sprintf("[%s] Drawing violin plot grouped by celltype.", set_name))
  
  # If your celltype column name is not "celltype", modify here
  celltype_col <- "celltype"
  
  # Extract data for plotting from meta.data
  df_plot <- data.frame(
    celltype = sce@meta.data[[celltype_col]],
    score    = sce@meta.data[[score_col]],
    stringsAsFactors = FALSE
  )
  
  # Fix celltype order (alphabetical or custom order)
  df_plot$celltype <- factor(df_plot$celltype,
                             levels = sort(unique(df_plot$celltype)))
  
  p_box <- ggplot(df_plot, aes(x = celltype, y = score, fill = celltype)) +
    geom_boxplot(outlier.size = 0.3, width = 0.8) +
    scale_fill_manual(values = col) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    ) +
    ylab(paste0(set_name, " UCell score")) +
    ggtitle(sprintf("%s score", set_name))
  
  
  ggsave(file.path(out_dir, paste0(set_name, "_VlnPlot_celltype.pdf")),
         plot = p_box, width = 4, height = 3.5)
  
  
  
  ## 5F. New: Draw boxplot of UCell score grouped by treatment with P-values
  treatment_col <- "treatment"
  
  df_treat <- data.frame(
    treatment = sce@meta.data[[treatment_col]],
    score     = sce@meta.data[[score_col]],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(score))
  
  # Fix treatment order
  df_treat$treatment <- factor(
    df_treat$treatment,
    levels = sort(unique(df_treat$treatment))
  )
  
  # Pairwise comparisons
  comps_treat <- combn(levels(df_treat$treatment), 2, simplify = FALSE)
  
  # Y-axis maximum (for P-value height)
  y_max_treat <- max(df_treat$score, na.rm = TRUE)
  
  p_box_treat <- ggplot(df_treat, aes(x = treatment, y = score, fill = treatment)) +
    geom_boxplot(
      outlier.size = 0.3,
      width = 0.8,
      alpha = 0.85
    ) +
    scale_fill_manual(values = c('#437eb8', "#CC0066", "#FF6666")) +
    
    ## ⭐ Significance test (Wilcoxon)
    stat_compare_means(
      comparisons = comps_treat,
      method = "wilcox.test",
      label = "p.signif",      # **** / ** / *
      size = 5,
      step.increase = 0.13     # Automatically offset for multiple group comparisons
    ) +
    
    ## ⭐ Reserve space for P-values
    ylim(NA, y_max_treat * 1.4) +
    
    theme_bw(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
      axis.text.y  = element_text(size = 14),
      panel.grid   = element_blank()
    ) +
    ylab("UCell score") +
    ggtitle(sprintf("%s score", set_name))
  
  # Save
  ggsave(
    file.path(out_dir, paste0(set_name, "_BoxPlot_treatment_pvalue.pdf")),
    plot  = p_box_treat,
    width = 3,
    height = 4
  )
  
}

# 6. Save updated object after loop completion (includes all gene set scores and classification columns)
saveRDS(sce, "celltype.rds")

# To view meta.data
View(sce@meta.data)








# ---- marker - Tissue preference analysis - Startrac---

# ---- Install/Load ----
# BiocManager::install("impute", ask = FALSE, update = FALSE)
# devtools::install_github("Japrin/sscVis")
# devtools::install_github("Japrin/Startrac")

library(Startrac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(grid)   # for grid.text

# =============== Basic settings ===============
if (!exists("outdir")) outdir <- "."     # If outdir is not defined, default to current directory
output_root <- file.path(outdir, "STARTRAC_marker")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

# Load object and meta data
sco  <- readRDS("celltype.rds")
data <- sco@meta.data

# Gene set names to loop through (corresponding to <set>_3class columns generated in previous UCell script)
gene_set_names <- c("Cancer","Ovarian","Digestive_tract","Metastatic")

# =============== Utility functions ===============
# Given a Ro/e matrix, return color mapping function and legend ticks
.make_color_meta <- function(mat) {
  rng <- range(mat, na.rm = TRUE)
  # Use 1 as visual midpoint (if not in range, fall back to range mean)
  mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng)
  col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]),
                                  c("#f6f8e6", "#f9a33e", "red"))
  legend_at <- pretty(rng, n = 5)
  legend_labels <- formatC(legend_at, format = "f", digits = 2)
  list(col_fun = col_fun, legend_at = legend_at, legend_labels = legend_labels)
}

# Determine if input is suitable for statistics (≥2 clusters & ≥2 tissues and not all NA)
.check_inputs <- function(df, clm_cluster, clm_tissue) {
  ok <- all(c(clm_cluster, "orig.ident", clm_tissue) %in% colnames(df))
  if (!ok) return(FALSE)
  n_cluster <- df[[clm_cluster]] %>% as.character() %>% unique() %>% length()
  n_tissue  <- df[[clm_tissue]]  %>% as.character() %>% unique() %>% length()
  n_cluster >= 2 && n_tissue >= 1
}

# =============== Main loop: calculate Roe for each gene set and plot ===============
for (set_name in gene_set_names) {
  message("Processing Roe for gene set: ", set_name)
  
  # Classification column name for this set, e.g., "Immunosuppression_3class"
  cluster_col <- paste0(set_name, "_3class")
  
  # Subdirectory for output
  output <- file.path(output_root, set_name)
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  # --------- Input check ---------
  if (!.check_inputs(data, cluster_col, "treatment")) {
    warning(sprintf("[%s] Required columns missing or insufficient categories (need >=2 clusters). Skipping.",
                    set_name))
    next
  }
  # If some samples have NA in this column, filter them out
  df_use <- data %>% filter(!is.na(.data[[cluster_col]]), !is.na(treatment))
  
  # If only 1 category left after filtering NA, skip
  if (length(unique(df_use[[cluster_col]])) < 2) {
    warning(sprintf("[%s] Only 1 cluster left after filtering NA, cannot calculate Roe. Skipping.", set_name))
    next
  }
  
  # --------- Calculate Ro/e ---------
  Roe_df <- tryCatch(
    calTissueDist(
      df_use,
      byPatient       = FALSE,
      colname.cluster = cluster_col,   # Key: group by this set's 3-class
      colname.patient = "orig.ident",
      colname.tissue  = "treatment",
      method          = "chisq",       # Options "chisq", "fisher", "freq"
      min.rowSum      = 0
    ),
    error = function(e) {
      warning(sprintf("[%s] calTissueDist error: %s", set_name, e$message))
      return(NULL)
    }
  )
  if (is.null(Roe_df) || nrow(Roe_df) == 0) {
    warning(sprintf("[%s] Roe result is empty, skipping plotting.", set_name))
    next
  }
  
  # Convert to matrix to avoid pitfalls with data frame min/max
  mat <- as.matrix(Roe_df)
  
  # If all zeros, skip plotting (or you could still plot an all-zero figure)
  if (all(is.na(mat)) || max(mat, na.rm = TRUE) == 0) {
    warning(sprintf("[%s] Roe is all 0 or NA, skipping plotting.", set_name))
    next
  }
  
  # Color mapping & legend
  col_meta <- .make_color_meta(mat)
  col_fun <- col_meta$col_fun
  legend_at <- col_meta$legend_at
  legend_labels <- col_meta$legend_labels
  
  # --------- Figure 1: Heatmap with values ---------
  pdf(file.path(output, paste0(set_name, "_STARTRAC_Roe_value.pdf")),
      width = 8, height = 3)
  p1<-Heatmap(
    mat,
    show_heatmap_legend = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "right",
    column_names_side = "bottom",
    show_column_names = TRUE,
    show_row_names = TRUE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    column_names_gp = gpar(fontsize = 16),
    heatmap_legend_param = list(
      title = "Ro/e value",
      at = legend_at,
      labels = legend_labels,
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 15)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 16, col = "black"))
    }
  )
  print(p1)
  dev.off()
  
  # --------- Figure 2: +++ symbol heatmap ---------
  pdf(file.path(output, paste0(set_name, "_STARTRAC_Roe_symbols.pdf")),
      width = 8, height = 3)
  p2<-Heatmap(
    mat,
    show_heatmap_legend = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "right",
    column_names_side = "bottom",
    show_column_names = TRUE,
    show_row_names = TRUE,
    col = col_fun,
    column_names_gp = gpar(fontsize = 16),
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Ro/e",
      at = c(0, max(mat, na.rm = TRUE)),
      title_gp = gpar(fontsize = 16),
      labels_gp = gpar(fontsize = 15),
      labels = c("0", "Max.")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- mat[i, j]
      symbol <- if (is.na(value) || value == 0) {
        "−"
      } else if (value > 0 & value < 0.2) {
        "+/−"
      } else if (value >= 0.2 & value <= 0.8) {
        "+"
      } else if (value > 0.8 & value <= 1) {
        "++"
      } else if (value > 1) {
        "+++"
      }
      # Background contrast: black text on light background, white on dark
      rgb_mean <- tryCatch(mean(col2rgb(fill)), error = function(e) 255)
      text_color <- ifelse(rgb_mean > 127, "black", "white")
      grid.text(symbol, x, y, gp = gpar(fontsize = 16, col = text_color))
    }
  )
  print(p2)
  dev.off()
  
  
}



