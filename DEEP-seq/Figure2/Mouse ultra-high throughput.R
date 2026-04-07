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

folders <- c("0822-WF-1/","0822-WF-2/","0822-WF-3/","0822-WF-4/")

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



save(ScRNA, file = "ScRNA（after batch correction before clustering）.RData")



#### 7. Cell clustering and annotation ####

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3','#6A4C93',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#6A4C93',"#99FFFF",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3','#6A4C93',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066")

load("ScRNA（after batch correction before clustering）.RData")

# Cell clustering
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:60,
          n.neighbors = 15,     # Default 30, decreasing makes clusters more separated
          min.dist = 0.1,       # Default 0.3, decreasing makes clusters more compact and separated
          spread = 2,
          seed.use = 42) %>%       # Increasing spread expands distance between clusters
  RunTSNE(dims = 1:60,          # Use more principal components
          perplexity = 50,      # Default is 30, can be increased appropriately
          theta = 0.3,          # More accurate calculation
          seed.use = 42) %>%
  FindNeighbors(dims = 1:50)


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


# Display clustering, ordered by Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_tsne.pdf"), width = 6*length(unique(ScRNA$treatment)), height = 5)
DimPlot(ScRNA, reduction = "tsne", pt.size=0.1, label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# Display clustering, ordered by Non-infected and Infected
pdf(paste(OUTPUT, "split.by_cluster_tsne_sample.pdf"), width = 6*length(unique(ScRNA$orig.ident)), height = 5)
DimPlot(ScRNA, reduction = "tsne", pt.size=0.1, label = TRUE, repel = TRUE, split.by = "orig.ident", cols = col)
dev.off()

# Generate separate umap plot
pdf(paste(OUTPUT, "cluster_tsne.pdf"), width = 7, height = 5)
DimPlot(ScRNA, reduction = "tsne",pt.size=0.1, label = TRUE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

DimPlot(ScRNA, reduction = "tsne", pt.size=0.1, label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
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

table(ScRNA$seurat_clusters)
sum(table(ScRNA$seurat_clusters))


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

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 180,height = 10)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # Increase X-axis title text size
        axis.title.y = element_text(size = 18))  # Increase Y-axis title text size
dev.off()

dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=150*dpi,h=7*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis(): tilt X-axis text
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 20),  # Increase X-axis title text size
        axis.title.y = element_text(size = 20))  # Increase Y-axis title text size
dev.off()



###########Manual cell annotation########

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("/out/")
outdir <- "/out/"

output <- paste(outdir,"celltype11", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
scedata <- readRDS(file_path)


# Define marker genes for different cell types
cellmarker <- c(
  
  ### Brain ###
  # Neuron (general)
  "Rbfox3","Tubb3","Map2",
  # Excitatory
  "Slc17a7","Slc17a6",
  # Inhibitory
  "Gad1","Gad2","Slc6a1",
  # Dopaminergic
  #"Th","Slc6a3","Ddc",
  # Cholinergic
  #"Chat","Slc5a7",
  # Astrocytes
  "Gfap","Aldh1l1","Aqp4","Slc1a3","Slc1a2",
  # Oligodendrocytes
  "Mbp","Plp1","Mog","Sox10","Matn4",
  # OPC
  "Brinp3","Pdgfra","Cspg4","Vcan","Megf11",
  # Microglia
  "P2ry12","Tmem119","Cx3cr1","Sall1","Aif1",
  # Brain endothelium   Endothelial/BBB
  "Slc2a1","Pecam1","Kdr",
  # Pericyte / VSMC
  "Pdgfrb","Rgs5","Acta2","Myh11",
  # Ependymal
  "Foxj1","S100b",
  # Choroid plexus
  "Ttr","Kcne2","Folr1","Clic6",
  "Hba-a1", "Hbb-bs", "Hba-a2", "Hbb-bt",    #Erythrocytes
  "Trbc1","Trbc2","Ccl5","Cd2", "Cd3e",  "Cd3g",    # T cells
  
  
  # Smooth muscle cell 
  
  ### Kidney ###
  # Podocyte
  "Nphs1","Nphs2","Podxl","Mafb","Wt1","Magi2","Thsd7a",
  # Glomerular Endothelium
  #"Pecam1","Kdr","Emcn",
  # Mesangial
  #"Pdgfrb","Itga8","Acta2",
  #"Slc8a1","Rhcg","Slc2a9",   # CNT
  #"Slc12a1","Ptger3","Cryab",    # ATL
  "Tshz2","Bst1","Slc4a11",      # DTL
  # Proximal Tubule
  "Lrp2","Cubn","Slc34a1","Ggt1",
  # Proximal segment
  "Cdh6",
  # Medullary cells
  "Col4a1","Akap12",
  # TAL
  "Umod","Slc12a1",
  # DCT
  "Slc12a3","Pvalb","Trpm6",
  # CDP Collecting duct Principal cells
  "Aqp2","Aqp3","Aqp4","Scnn1g","Avpr2","Apela","Atp1b1","Cdh16","Fxyd2",
  # TECs Tubular epithelial cells
  "Aqp1",
  # IC A Intercalated Cell Type A
  "Atp6v1b1","Slc4a1","Pam","Adgrf5","Kit", 
  # IC B
  "Atp6v1b1","Slc26a4","Slc4a9","Atp6v1g3","Hepacam2",
  # JG / Macula densa
  #"Ren1","Nos1","Ptgs2",
  # Interstitium / Pericyte
  #"Pdgfra","Pdgfrb","Rgs5","Col1a1",
  #"Pdgfrb","Itga8","Des","Igfbp7"   # MC Mesangium cells
  # Immune in kidney
  #"Adgre1","Cd68","Lyz2","Cd3e","Cd4","Cd8a","Ms4a1","Cd79a",
  
  
  
  
  ### Liver ###
  # LSEC  Liver sinusoid endothelial cell
  "Stab2","Mrc1","F8","Clec4g",  
  # Hepatocytes
  "Alb","Apoa1","Apoa2","Tdo2",
  # Kupffer
  "Clec4f","Timd4","Vsig4","Adgre1",
  # Neutrophils / Monocytes
  "S100a8","S100a9",
  # Endothelial
  "Klf2","Car4",
  # Dendritic cells
  "Zbtb46","Xcr1","Itgam",
  "Cd79a","Cd79b", "Ms4a1",       # B cells
  # Macrophages
  "Cd68","Cd163","Adgre1","Ly6C","Apoe",
  "Trbc2","Cd3e",  "Cd3g",    # T cells
  "Ly6c2","Ccr2",  # Monocytes
  # Cholangiocytes
  "Krt19","Krt7","Epcam","Sox9",
  "Pkhd1","Epcam",   # Cholangiocyte
  "Gzma", "Ifng",   # Cytotoxic T cells
  "Mzb1","Jchain",                     # Plasma Cells
  # pDCs
  "Bst2", "Siglech", "Tcf4","Cox6a2",
  "Gzmb",   # Natural Killer (NK) cells
  
  
  
  ### Lung ###
  # AT1
  "Ager","Pdpn","Hopx",
  # AT2
  "Sftpc","Sftpa1","Sftpb","Abca3",
  "Sulf1",                 # Tuft cells
  # Club cells
  "Scgb1a1","Scgb3a2",
  "Sulf1",                 # Tuft cells
  "Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells
  # Ciliated cells
  "Foxj1","Tekt2","Dynlrb2",
  "Foxj1","Tppp3","Tubb4b","Tubb1",   #Cilliated
  # Endothelial
  "Pecam1","Kdr","Klf2","Car4",
  # Lymphatic Endothelium
  "Prox1","Lyve1","Pdpn",
  # Matrix fibroblast
  "Gpx3","Col25a1",
  # Fibroblast
  "Col1a1","Col1a2","Pdgfra","Prg4","Dcn","Acta2",
  # Smooth muscle / Pericytes
  "Acta2","Myh11","Tagln","Pdgfrb","Rgs5",
  # Alveolar Macrophage
  "Lyz2","Siglecf","Itgax","Pparg","Marco",
  # Interstitial Macrophage
  "Mertk","Lyve1","Cd68","Adgre1","Spp1",
  # Dendritic cells
  "Flt3","Zbtb46","Xcr1","Itgam",
  # Neutrophils / Monocytes
  "S100a8","S100a9"
  
  
  
  
  
  
  
  
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
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 25, height = 6)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 25, height = 6)







# Annotate cell types for clusters
scedata <- RenameIdents(scedata, c(
  "0"="Lung",
  "1"= "Kidney",
  "2"="Brain", 
  "3"="Liver",
  "4"="Liver",
  "5"= "Kidney",
  "6"="Brain",  
  "7"="Brain", 
  "8"="Brain", 
  "9"="Kidney",
  "10"="Liver", 
  "11"="Liver",
  "12"="Lung",
  "13"="Liver",
  "14"="Liver",
  "15"="Liver",
  "16"="Lung",
  "17"="Brain", 
  "18"="Liver",
  "19"="Brain", 
  "20"="Brain", 
  "21"="Liver",
  "22"="Lung",
  "23"="Liver",
  "24"="Brain", 
  "25"="Kidney",
  "26"= "Kidney",
  "27"="Lung",
  "28"="Lung",
  "29"="Liver",
  "30"= "Liver",
  "31"="Liver",
  "32"="Lung",
  "33"="Kidney",
  "34"="Lung", 
  "35"="Lung",
  "36"="Liver",
  "37"="Brain", 
  "38"="Brain", 
  "39"="Lung",
  "40"="Liver",
  "41"="Lung"
  
)
)





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
  "#2e87b6","#ac352f","#501a92","#fee089","#FF3366","#CC0066"
)

file_path <- file.path(outdir, "ScRNA（after clustering）.rds")
scedata <- readRDS(file_path)

# Annotate cell types for clusters
scedata <- RenameIdents(scedata, c(
  "0"="Fibroblasts",
  "1"= "Proximal Tubule",
  "2"="Oligodendrocytes", 
  "3"="Hepatocytes",
  "4"="LSECs",
  "5"= "Proximal Tubule",
  "6"= "Astrocytes", 
  "7"="Microglia", 
  "8"="Erythrocytes", 
  "9"="TAL",
  "10"="Neutrophils", 
  "11"="Kupffer",
  "12"="Fibroblasts",
  "13"="T Cells",
  "14"="Endothelial Cells",
  "15"="DCs",
  "16"="Alveolar macrophages",
  "17"="Brain endothelium", 
  "18"= "B cells",
  "19"="Astrocytes", 
  "20"="Choroid plexus cells",
  "21"="Endothelial Cells",
  "22"= "Tuft cells",
  "23"="Macrophages",
  "24"="Neurons", 
  "25"="IC",
  "26"= "CDP",
  "27"="AT2",
  "28"="Smooth muscle",
  "29"="Monocytes", 
  "30"= "Endothelial Cells",
  "31"="Cholangiocytes",
  "32"="Club cells",
  "33"="DCT",
  "34"="Matrix fibroblast", 
  "35"="Cilliated",
  "36"= "Plasma Cells",
  "37"="Pericytes",
  "38"="OPC", 
  "39"="AT1",
  "40"= "LSECs",
  "41"="Fibroblasts"
  
)
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
pdf(paste(output, "ann_umap.pdf",sep = '/'), width = 6, height = 4)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=0.1,label=TRUE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        #legend.position = "none",
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18),
        plot.title = element_blank())

dev.off() 

#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom"))


# Plot umap with cell types
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


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=7*length(unique(scedata$treatment)),height=7)
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

svg(paste(output, "ann-diff-umap.svg",sep = '/'),width=7*length(unique(scedata$treatment)),height=7)
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

# Filter data from scedata where orig.ident is T1, T2, T3, T4
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







########## Add overall total cell count ##########
library(ggsci)

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

# Plot umap with cell types and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 6, height = 6)
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
svg(paste(output, "ann_umap.svg", sep = '/'), width = 6, height = 6)
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




########## Add overall total cell count ##########

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

library(ggsci)
# Plot umap with cell types
pdf(paste(output, "ann_tsne.pdf",sep = '/'), width = 11, height = 6)
DimPlot(object=scedata,group.by = "celltype",reduction='tsne',pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title = element_text(face = 2,hjust = 0.03),
        #legend.position = "none",
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18))+
  ggtitle(title_label)

dev.off() 

#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom"))


# Plot umap with cell types
svg(paste(output, "ann_tsne.svg",sep = '/'), width = 11, height = 6)
DimPlot(object=scedata,group.by = "celltype",reduction='tsne',pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title = element_text(face = 2,hjust = 0.03),
        #legend.position = "none"
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18),
  )+
  ggtitle(title_label)
#        legend.position = c(0.99, 0.12),  # Move legend to bottom right
#        legend.justification = c("right", "bottom")) +
dev.off()





########## Add overall total cell count ##########
library(ggsci)

# Calculate total cell count
total_cells <- ncol(scedata)

# Construct title
title_label <- paste("( n =", total_cells, "cells )")

# Plot umap with cell types and save as PDF
pdf(paste(output, "ann_tsne.pdf", sep = '/'), width = 6, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'tsne', pt.size = 0.1, label = TRUE, 
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
svg(paste(output, "ann_tsne.svg", sep = '/'), width = 6, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'tsne', pt.size = 0.1, label = TRUE, 
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
col <- c("#41B6C4","#FECACA","#DADAEB","#FDBE85","#A4CDE1","#FCA5A5","#CCCCFF","#FF9966","#66CCCC","#FFCCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066","#9966CC",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)
#View(scedata@meta.data)

table(scedata$celltype)
sum(table(scedata$celltype))


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
ggsave(file_path, plot = p, width = 6*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 6*length(unique(scedata$orig.ident)), height = 6, dpi = 800)



###############Plot each sample individually########################
library(ggplot2)
library(dplyr)

# Define samples to plot
samples_to_plot <- c("WF-14W")

# Create output directory (if it doesn't exist)
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

for (sample_id in samples_to_plot) {
  
  # Filter current sample data
  cell_counts_group_filtered <- subset(cell_counts_group, Sample == sample_id)
  
  ### —— Cell count plot —— ###
  # Aggregate and generate legend labels (cell counts)
  cell_counts_group_agg <- aggregate(Counts ~ CellType, cell_counts_group_filtered, sum)
  cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", cell_counts_group_agg$Counts, ")")
  legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)
  
  # Create count bar plot
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
  ggsave(paste0(output, "/genecount_", sample_id, ".pdf"), plot = p_counts, width = 7, height = 8, dpi = 800)
  ggsave(paste0(output, "/genecount_", sample_id, ".svg"), plot = p_counts, width = 7, height = 8, dpi = 800)
  
  ### —— Cell ratio plot —— ###
  # Aggregate and generate legend labels (proportions)
  cell_counts_group_agg <- aggregate(Ratio ~ CellType, cell_counts_group_filtered, sum)
  cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", round(cell_counts_group_agg$Ratio * 100, 2), "%)")
  legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)
  
  # Create ratio bar plot
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
  
  # Save cell ratio plot
  ggsave(paste0(output, "/geneRatio_", sample_id, ".pdf"), plot = p_ratio, width = 7, height = 8, dpi = 800)
  ggsave(paste0(output, "/geneRatio_", sample_id, ".svg"), plot = p_ratio, width = 7, height = 8, dpi = 800)
}




############Group by treatment############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate proportion of each cell type within each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

# Sort
#cell_counts_treatment$Treatment <- factor(cell_counts_treatment$Treatment, levels = c("0d", "3d","7d", "14d"))

########## Create stacked bar plot for cell counts ##########
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
ggsave(file_path, plot = p1, width = 6.5*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 6.5*length(unique(scedata$treatment)), height = 7, dpi = 800)

########## Create stacked bar plot for cell proportions ##########
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
ggsave(file_path, plot = p2, width = 6*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 6*length(unique(scedata$treatment)), height = 7, dpi = 800)






## Merge into one group ###
library(ggplot2)
library(dplyr)

# Merge all treatment samples into one group
cell_counts_group$treatment_group <- ifelse(cell_counts_group$Sample %in% c("WF-14W", "WF-15W", "WF-16W"), "Treatment", "Control")

# Calculate merged cell counts
cell_counts_group_agg <- cell_counts_group %>%
  group_by(treatment_group, CellType) %>%
  summarise(Counts = sum(Counts), .groups = 'drop')

# Calculate cell proportions
total_counts <- sum(cell_counts_group_agg$Counts)
cell_counts_group_agg$Ratio <- cell_counts_group_agg$Counts / total_counts

# Create stacked cell count plot
p_counts_stacked <- ggplot(cell_counts_group_agg, aes(x = treatment_group, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))

# Save stacked cell count plot
ggsave(paste0(output, "/cell_counts_stacked.pdf"), plot = p_counts_stacked, width = 5, height = 8, dpi = 800)
ggsave(paste0(output, "/cell_counts_stacked.svg"), plot = p_counts_stacked, width = 5, height = 8, dpi = 800)

# Create stacked cell ratio plot
p_ratio_stacked <- ggplot(cell_counts_group_agg, aes(x = treatment_group, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))

# Save stacked cell ratio plot
ggsave(paste0(output, "/cell_ratio_stacked.pdf"), plot = p_ratio_stacked, width = 5, height = 8, dpi = 800)
ggsave(paste0(output, "/cell_ratio_stacked.svg"), plot = p_ratio_stacked, width = 5, height = 8, dpi = 800)



