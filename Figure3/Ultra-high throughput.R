
# 清空环境
rm(list = ls())

# 加载所需的包
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
#install.packages("tidydr")


setwd("D:/R/GS/WH/20250716-超高通/data/")
# 定义文件夹路径
data_dir <- "D:/R/GS/WH/20250716-超高通/data/"

folders <- c("0627-WF-RT-14W/","1128-LC-3/")

# 初始化 Seurat 对象列表
seurat.list <- list()

# 循环读取每个子文件夹中的表达矩阵文件并进行 Seurat 分析
for (folder in folders) {
  folder_path <- file.path(data_dir, folder)
  
  # 读取10X数据格式
  sample_data <- Read10X(data.dir = folder_path)
  sample_name <- basename(folder)
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 3, min.features = 200)
  
  # 添加线粒体基因比例
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-|^GRCh38_MT-|^GRCm39_mt-")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hbb-|^Hba-")
  
  # 质控过滤
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
  
  # 数据标准化与高变基因寻找
  #seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 根据文件夹名称分配 orig.ident 分组
  if (str_detect(sample_name, "WF-RT-14W")) {
    seurat_obj$orig.ident <- "WF-14W"
    seurat_obj$treatment <- "WF-14W"
  } else if (str_detect(sample_name, "LC-3")) {
    seurat_obj$orig.ident <- "Ovarian cancer"
    seurat_obj$treatment <- "Ovarian cancer"
  }else if (str_detect(sample_name, "PC_")) {
    seurat_obj$orig.ident <- "Pancreatic cancer"
    seurat_obj$treatment <- "Pancreatic cancer"
  }else if (str_detect(sample_name, "GC_")) {
    seurat_obj$orig.ident <- "Gastric cancer"
    seurat_obj$treatment <- "Gastric cancer"
  }
  
  # 添加 treatment 列
  #seurat_obj$treatment <- sample_name 
  
  
  # 将处理后的 Seurat 对象加入列表
  seurat.list[[sample_name]] <- seurat_obj
}

# 整合数据
#anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)
#combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:20)

combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)
combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
View(combined_seurat@meta.data)

# 保存合并后的 Seurat 对象
saveRDS(combined_seurat, file = "D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/combined_seurat.rds")


#######################Seurat分析#####################
# 设置输出目录
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/")
outdir <- "D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# 拼接完整路径
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))
#ScRNA$`orig.ident` <- factor(ScRNA$`orig.ident`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))

## 计算红细胞比例
ScRNA[["percent.hb"]] <- PercentageFeatureSet(ScRNA, pattern = "^Hbb-|^Hba-")

# 生成小提琴图，显示质控指标
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 5)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4, group.by = "treatment", pt.size = 0,cols = col)
dev.off()
## 去除红细胞污染
ScRNA <- subset(ScRNA, subset = percent.hb < 1)

# 生成小提琴图，显示质控指标
svg(paste(OUTPUT, "QC-BoxPlot.svg"), width = 8, height = 6)
p1 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nFeature_RNA, color = treatment)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nFeature_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p2 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nCount_RNA, color = treatment)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nCount_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

CombinePlots(plots = list(p1, p2))
dev.off()


svg(paste(OUTPUT, "QC-ViolinPlot-sample.svg"), width = 9, height = 6)
# 小提琴图1：nFeature_RNA
p1 <- ggplot(data = ScRNA@meta.data, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values =col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 小提琴图2：nCount_RNA
p2 <- ggplot(data = ScRNA@meta.data, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p3 <- ggplot(data = ScRNA@meta.data, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "percent.mt", x = "", y = "") +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 合并图像
library(patchwork)
(p1 | p2 |p3) + plot_layout(ncol = 3)
dev.off()



#### 5.表达量标准化 ####
ScRNA <- NormalizeData(ScRNA, normalization.method = "LogNormalize",
                       scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
ScRNA <- FindVariableFeatures(ScRNA, selection.method = "vst",
                              nfeatures = 2000)

#展示高变基因
pdf(paste(OUTPUT,"variable gene.pdf"),width = 9,height = 6)
top10 <- head(VariableFeatures(ScRNA), 10) 
plot1 <- VariableFeaturePlot(ScRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=3)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()

#### 6.归一化与PCA降维 ####
#归一化
ScRNA<-ScaleData(ScRNA)

#运行PCA
ScRNA<-RunPCA(ScRNA,npcs = 30)

pdf(paste(OUTPUT,"Dimplot.pdf"),width = 9,height = 6)
DimPlot(object = ScRNA, reduction = "pca", pt.size = .1, group.by = "treatment",cols = col)
dev.off()

pdf(paste(OUTPUT,"vlnplot.pdf"),width = 9,height = 6)
VlnPlot(object = ScRNA, features = "PC_1", group.by = "treatment", pt.size = 0,cols = col)
dev.off()

#PCA展示
pdf(paste(OUTPUT, "DimHeatmap.pdf"),width = 9,height = 6)
DimHeatmap(ScRNA, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# 评估PC维度
pdf(paste0(OUTPUT,"PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(ScRNA)
dev.off()


#批次矫正
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


save(ScRNA, file = "ScRNA（批次矫正后分群前）.RData")



#### 7.细胞分群与注释 ####

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#6A4C93',
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA（批次矫正后分群前）.RData")
#细胞分群
#ScRNA <- ScRNA %>% 
#  RunUMAP(dims = 1:20) %>% 
#  RunTSNE(dims = 1:20) %>%
#  FindNeighbors(dims = 1:20)

ScRNA <- ScRNA %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,spread = 2) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)

ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))
#pdf(paste(OUTPUT, "clustree.pdf"),width=10,height=9)
#library(clustree)
#clustree(ScRNA)
#dev.off()

#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.0.7"
ScRNA$seurat_clusters <- ScRNA@active.ident##根据聚类树选择你要的resolution
table(Idents(ScRNA))

#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))

# 确保 "treatment" 因子水平按照 Non-infected 和 Infected 顺序排列
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 6*length(unique(ScRNA$treatment)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=0.1, label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 6*length(unique(ScRNA$orig.ident)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=0.1, label = TRUE, repel = TRUE, split.by = "orig.ident", cols = col)
dev.off()

# 单独生成umap图
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

#将肿瘤与正常展示在一起
pdf(paste(OUTPUT, "cluster-diff_umap.pdf"),width=6,height=6)
DimPlot(ScRNA, repel = TRUE,pt.size=0.1, 
        reduction = "umap",
        group.by ="treatment")+
  scale_color_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")
dev.off()


saveRDS(ScRNA, "ScRNA（分群后）.rds")



# 选取前6细胞群展示标记基因的表达趋势
output <- paste(outdir,'细胞定位', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)

#寻找Marker基因
ScRNA.markers <- FindAllMarkers(ScRNA, only.pos = TRUE,   ### only.pos = TRUE：只寻找上调基因
                                min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ScRNA.markers,paste0(output,"./ScRNA.all.markers.csv"))

dim.use <-1:30
top5 <- ScRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top5,file=paste0(output,"/top20_marker_genes_tsne_",max(dim.use),"PC.csv"))

pdf(paste0(output,"/Heatmap_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 25,height = 20)
DoHeatmap(ScRNA, features = top5$gene,size = 2)+
  scale_fill_gradientn(colors = c("#437eb8", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # 增大X轴标题文本大小
        axis.title.y = element_text(size = 18),  # 增大Y轴标题文本大小
        legend.text = element_text(size = 14),  # 修改图例文本大小
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 30))  # 修改分组标签的颜色和大小
dev.off()

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 80,height = 7)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis()：倾斜X轴文本
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # 增大X轴标题文本大小
        axis.title.y = element_text(size = 18))  # 增大Y轴标题文本大小
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=80*dpi,h=7*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis()：倾斜X轴文本
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 20),  # 增大X轴标题文本大小
        axis.title.y = element_text(size = 20))  # 增大Y轴标题文本大小
dev.off()



###########细胞手动注释########

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/")
outdir <- "D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)


# 定义不同细胞类型的marker基因
cellmarker <- c(
  "PECAM1", "VWF", "CDH5",             # Endothelial Cells
  #"PROX1","LYVE1", "CCL21",  # Lymphatic endothelial cells 淋巴内皮细胞 
  "TPSAB1", "TPSB2",                   # Mast cells 肥大细胞
  "PCNA", "TOP2A", "CDK1",            # Proliferating cells 增殖细胞
  "COL1A1", "COL1A2", "DCN",                    # Fibroblasts / Stromal
  #"ACTA2" , "MYH11", "CNN1",                           # Smooth muscle cells
  "SFRP4","FAP","MMP11","PDGFRA","PDPN","FSP1",   # CAFs 肿瘤相关成纤维细胞
  "PDGFRB", "RGS5",   # Pericytes  周细胞
  
  #"MALAT1",    # Tfh (T follicular helper)
  
  "KRT19", "MUC1", "KRT8","CDH1",         # Epithelial Cells /EOC上皮性卵巢癌
  "EPCAM","MUC16","WFDC2","KRT7","PAX8","CLDN4","CD24","TP53","WT1",      # Cancer Cells
  
  #"HBA1","HBA2","HBB",    # Erythroid Cells
  #"CD79A", "MS4A1","IGKC","IGHA1",                 # B cells
  "MZB1","JCHAIN",   # Plasma Cells 浆细胞
  #"GZMB","LILRA4","SPIB","ZFAT",         # pDCs 浆细胞样树突状细胞
  "HLA-DPB1", "HLA-DRA", "HLA-DRB1", # Dendritic cells
  "CD68", "APOE", "LGALS3", "ITGAM", "PPARG" ,         # Macrophages
  "FCN1","CD300E", "NLRP3" ,"TBC1D8",               # Monocytes 单核细胞
  #"ITGAM","CD14","CD33","S100A8","S100A9","ARG1","NOS2","IL4R","IL10","ANXA1","LOX1",           # 髓源性抑制细胞（MDSCs）
  "S100A8", "S100A9", "G0S2",     # Neutrophils
  
  "GZMA", "GZMB",  "IFNG", "CCL4","CCL5",  # Cytotoxic T cells 细胞毒性T细胞
  #"NKG7", "CCL5", "KLRB1","GZMA","KLRF1", "PRF1",     # Natural Killer (NK) cells
  "H3F3A","ABCF1","INTS6","TRNAU1AP","ERO1B","ERCC1","CEMIP2",      # NK T cells   "BRCA1", "MCM6", "HELLS", "CDT1", "DTL",
  #"CD8A","CD8B", # CD8+ T 
  "TRBC1","TRBC2","CCL5","CD2", "CD3E",  "CD3G"    # T cells
  #"CCR7",  "CD3D","CD3E","CD4","CD8A","SELL","TCF7","LEF1",  # Naive T Cells 初始T细胞
  
  #### 肺癌 ####
  #"PECAM1", "FLT1", "VWF", "CDH5", "CA4",              # Endothelial Cells
  #"PROX1", "FLT4", "PDPN", "LYVE1", "CCL21", "ITGA9", "NRP2", "MRC1", "SOX18", "FOXC2", # Lymphatic endothelial cells 淋巴内皮细胞 
  #"CALB2","WT1","MSLN","UPK3B","PDPN",     # Mesothelial Cells
  #"PDGFRB", "CSPG4", "RGS5", "ANPEP", "ABCC9", "NOTCH3", "KCNJ8", # Pericytes
  #"PSD3", "FOXI1", "CFTR", "ATP6V0D2", "ATP6V1B1", "ASCL3", "TP63", "EZR",  # Pulmonary ionocytes 非离子细胞
  #"CTSD", "CTSS", "ITGAL", "RBM47","PLXDC2", "FCN1","C5AR1","CD300E", "NLRP3" ,"TBC1D8",  Monocytes 单核细胞
  #"TPSAB1", "TPSB2",                                   # Mast cells 肥大细胞
  #"PCNA", "TOP2A", "CCNA2", "MCM5", "CDK1",            # Proliferating cells 增殖细胞
  #"HLA-DPB1", "HLA-DRA", "NAAA", "GM2A", "HLA-DRB1", "PPT1", "CYTIP", # Dendritic cells
  #"COL1A1", "COL1A2", "DCN", "FN1",                    # Fibroblasts / Stromal
  #"MYH11", "CNN1", "SMTN",                             # Smooth muscle cells
  #"HBA1","HBA2","HBB",    # Erythroid Cells
  #"CD79A", "CD79B", "MS4A1",                           # B cells
  #"CD38","MZB1","IRF4","JCHAIN","SDC1",   # Plasma Cells 浆细胞
  #"GZMB","IL3RA","CLEC4C","LILRA4","TCF4","IRF7","PTCRA","SPIB","NRP1","IRF8",           # pDCs 浆细胞样树突状细胞
  #"CD68", "LGALS3", "ITGAM", "APOE",                   # Macrophages
  #"MARCO","ABCA1","MRC1","CCL3", "FBP1",  "MCEMP1",       # Alveolar Macrophages
  #"ARG1","MRC1",  ## M2 Macrophages
  #"IL2RB", "NKG7", "ADAMTS14", "KLRA4", "KLRF1", "PRF1", # Natural Killer (NK) cells
  #"H4C3","SMC4","HMGN2","HMGB2","NCALD","KCNQ5","PDGFD","ZNF331","B4GALT1", "PLCB1", "BRCA1", "MCM6", "HELLS", "CDT1", "DTL", "UNG", "RMI2", # NK T cells
  #"CD2", "CD3E", "CD3G", "CCR7", "NKG7", "KLRB1", "CD3G", "TRBC2", # T cells
  #"IL2RA","FOXP3","CTLA4","IKZF2","TNFRSF18" ,     # Treg 
  #"ICOS", "GRAP2", "CSF2", "GATA3", "PDCD1",           # T helper cells (Th cells)
  #"CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "NKG7","CCL4", "CCL5", "IL2", "TBX21", # Cytotoxic T cells 细胞毒性T细胞
  #"CCR7",  "CD3D","CD3E","CD4","CD8A","SELL","TCF7","LEF1",  # Naive T Cells 初始T细胞
  #"SNAP25", "SYP", "SYT1", "DCX", "SYN1", "SLC17A7",   # Neurons
  #"S100A8", "S100A9", "G0S2", "LY6G6C", "MPO", "CSF3R","BCL2A1","IL1R2","FCGR3B","CSF3R",    # Neutrophils
  #"SPINK5","NXN","TIMP1","KRT5","TP63","KRT14","ITGA6","ITGB4","NGFR",   # Basal Cells
  #"FOXJ1", "TPPP3", "TUBB4B",  "TP73", "CCDC7", # Ciliated cells 纤毛细胞
  #"SCGB1A1", "SCGB3A2", "CALCB", "NOTCH2", "HES1",     # Club cells
  #"CALCA", "SYP", "RESP18", "PCSK1", "SCG5", "CHGB", "SEZ6L2",   # PNEC (肺神经内分泌细胞)
  #"FMNL2","PCSK2","CACNA2D1",     ### 神经元内分泌细胞 Neuroendocrine cells
  #"EPCAM", "KRT18", "CD24", "KRT19", "FGFR2", "SPRR3","SPRR2A","SCEL",          # Epithelial Cells /EOC上皮性卵巢癌
  #"STEAP4","CEACAM6","SCGB1A1","MUC5B","MUC5AC","SPDEF","FOXJ1",    # Secretory Epithelial Cells 分泌性上皮细胞
  #"FOXC1", "TP73", "DNAH5", "DNAH9", "CCDC39", "CCDC40", "TEKT1", "RSPH4A", "TUBB4B", "SPEF2", # Ciliated epithelial cells 纤毛上皮细胞
  #'AGER', 'CAV1', 'CLIC5' ,'HOPX', 'SEMA3E', 'COL4A3',  #AT1
  #"LAMP3", "ABCA3", "SLC34A2", "LPCAT1", "SFTPC", "SFTPA1", "SFTPB", "SFTPD", "AGER", "CLDN18", "NKX2-1", "MUC1", "KRT8" # AT2 cells
  
)




cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# 使用DotPlot可视化免疫细胞marker基因的表达
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

# 保存DotPlot图
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
  # 深蓝→绿→浅绿 梯度
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq 浅色
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)

# 对Cluster进行细胞类型注释
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

# 将细胞类型添加到meta数据中
scedata$celltype <- scedata@active.ident
head(scedata@meta.data)

# 提取UMAP坐标和细胞类型
umap_coords <- as.data.frame(Embeddings(scedata, reduction = "umap"))
umap_coords$celltype <- scedata$celltype
#umap_coords$CB <- scedata$CB

# 保存为CSV文件
write.csv(umap_coords, "celltype_umap.csv", row.names = TRUE)

saveRDS(scedata,  "celltype.rds")


file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

library(ggsci)
# 绘制细胞类型的umap图
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

#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom"))


# 绘制细胞类型的umap图
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
#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=6*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # 增大子图标题字体
    axis.text.x = element_text(size = 16),  # X轴标签大小
    axis.text.y = element_text(size = 16),  # Y轴标签大小
    axis.title.x = element_text(size = 18, face = "bold"),  # 增大X轴标题大小
    axis.title.y = element_text(size = 18, face = "bold"),  # 增大Y轴标题大小
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # 增大标题大小
    legend.title = element_text(size = 18),  # 增大图例标题大小
    legend.text = element_text(size = 18)    # 增大图例文本大小
  )
dev.off()

svg(paste(output, "ann-diff-umap.svg",sep = '/'),width=6*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # 增大子图标题字体
    axis.text.x = element_text(size = 16),  # X轴标签大小
    axis.text.y = element_text(size = 16),  # Y轴标签大小
    axis.title.x = element_text(size = 18, face = "bold"),  # 增大X轴标题大小
    axis.title.y = element_text(size = 18, face = "bold"),  # 增大Y轴标题大小
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # 增大标题大小
    legend.title = element_text(size = 18),  # 增大图例标题大小
    legend.text = element_text(size = 18)    # 增大图例文本大小
  )
dev.off()

# 筛选需要的样本
#selected_samples <- c("T1", "T2", "T3", "T4")

# 从scedata中筛选orig.ident为T1、T2、T3、T4的数据
#filtered_data <- subset(scedata, subset = orig.ident %in% selected_samples)

# 绘制并保存为PDF
#pdf(paste(output, "ann-diff-umap-sample-selected.pdf", sep = '/'), width = 22, height = 5)
#DimPlot(filtered_data, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, 
#        label = FALSE, label.size = 5, repel = TRUE, cols = col)
#dev.off()

# 绘制并保存为SVG
#svg(paste(output, "ann-diff-umap-sample-selected.svg", sep = '/'), width = 22, height = 5)
#DimPlot(filtered_data, reduction = "umap", split.by = "orig.ident", pt.size = 0.1, 
#        label = FALSE, label.size = 5, repel = TRUE, cols = col)
#dev.off()





library(ggsci)
# 绘制细胞类型的umap图
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

#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom"))


# 绘制细胞类型的umap图
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
#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=12,height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)
dev.off()

svg(paste(output, "ann-diff-umap.svg",sep = '/'),width=12,height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=0.1,label=FALSE,label.size = 5,repel = TRUE,cols=col)
dev.off()



########## 添加整体的细胞总数 #########
library(ggsci)

# 计算细胞总数
total_cells <- ncol(scedata)

# 构建标题
title_label <- paste("( n =", total_cells, "cells )")

# 绘制细胞类型的umap图并保存为PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 5, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + # 居中显示标题
  ggtitle(title_label)
dev.off()

# 绘制细胞类型的umap图并保存为SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 7, height = 6)
DimPlot(object = scedata, group.by = "celltype", reduction = 'umap', pt.size = 0.1, label = TRUE, 
        label.size = 5, repel = TRUE, cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + # 居中显示标题
  ggtitle(title_label)
dev.off()



###### 添加细胞标签 ##########
library(ggsci)
library(Seurat)
library(dplyr)
library(ggrepel)

# 计算细胞总数
total_cells <- ncol(scedata)

# 构建标题
title_label <- paste("( n =", total_cells, "cells )")

# 设置细胞类型的颜色，假设 col 是一个包含颜色的向量，已经与细胞类型对应
celltype_colors <- col  # 假设 col 是预先定义的颜色向量

# 提取UMAP坐标
umap_coords <- Embeddings(scedata, "umap")  # 提取UMAP坐标
umap_data <- as.data.frame(umap_coords)
umap_data$celltype <- scedata$celltype  # 添加细胞类型信息

# 创建一个含有每个celltype的中心点坐标的数据框
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# 计算每个celltype的中心点坐标
celltype_centers <- umap_df %>%
  group_by(celltype) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# 计算每个细胞类型的边界范围，用于确定标签放置方向
celltype_ranges <- umap_df %>%
  group_by(celltype) %>%
  summarise(
    min_x = min(UMAP1), max_x = max(UMAP1),
    min_y = min(UMAP2), max_y = max(UMAP2),
    width = max_x - min_x, height = max_y - min_y
  )

# 合并中心点坐标和边界信息
celltype_labels <- left_join(celltype_centers, celltype_ranges, by = "celltype")

# 确定标签放置方向：基于细胞群的形状
celltype_labels <- celltype_labels %>%
  mutate(
    # 根据细胞群的宽高比决定标签放置方向
    direction_x = ifelse(width > height, 1, 0),
    direction_y = ifelse(height > width, 1, 0),
    # 设置标签的nudge参数：将标签推向细胞群的外侧
    nudge_x = ifelse(UMAP1 > mean(UMAP1), 2, -2),  # 右侧细胞群向左推，左侧向右推
    nudge_y = ifelse(UMAP2 > mean(UMAP2), 2, -2)   # 上方细胞群向下推，下方向上推
  )

# 绘制细胞类型的umap图并保存为PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1) +
  ggrepel::geom_text_repel(
    data = celltype_labels, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 7, 
    fontface = "bold",
    box.padding = 1.5,        # 增加边框填充，让标签远离点
    point.padding = 0.8,      # 增加点填充
    nudge_x = celltype_labels$nudge_x,  # 水平方向偏移
    nudge_y = celltype_labels$nudge_y,  # 垂直方向偏移
    min.segment.length = 1,   # 增加最小线段长度
    segment.size = 0.5,       # 线段粗细
    segment.color = "grey40", # 线段颜色
    segment.alpha = 0.7,      # 线段透明度
    force = 2,                # 增加排斥力
    max.iter = 10000,         # 增加最大迭代次数
    direction = "both",       # 允许双向调整
    seed = 123,               # 设置随机种子保证可重复性
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

# 绘制细胞类型的umap图并保存为SVG
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




##### 添加各组的细胞总数 #####
# 统计每个 treatment 的细胞数
cell_counts <- scedata@meta.data %>%
  group_by(treatment) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(treatment, " ( n = ", n, " cells)"))

# 构建命名向量用于替换 facet 标签
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# 绘制 PDF
pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),
    width = 6.5*length(unique(scedata$treatment)), height = 5.5)

DimPlot(scedata, reduction = "umap", split.by = "treatment",
        pt.size = 0.1, label = FALSE, cols = col) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),  # 子图标题加粗黑体
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
dev.off()

# 同理，SVG 输出
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






## 将rds转变为h5ad
library(Seurat)
library(SeuratDisk)

# 1. 读取 Seurat 对象
sc_data <- readRDS("celltype.rds")

# 2. 【关键】把 factor 转为 character
sc_data$celltype  <- as.character(sc_data$celltype)
sc_data$treatment <- as.character(sc_data$treatment)

# 可选：确认
str(sc_data@meta.data$celltype)

# 3. 设置 assay
DefaultAssay(sc_data) <- "RNA"

# 5. 精简对象 —— 显式保留 umap & pca
sc_data <- DietSeurat(
  sc_data,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = c("pca", "umap")  # ⭐ 核心修改
)

# 5. 保存并转换
SaveH5Seurat(
  sc_data,
  filename = "ScRNA.h5seurat",
  overwrite = TRUE
)

Convert(
  "ScRNA.h5seurat",
  dest = "h5ad",
  overwrite = TRUE
)






############## 载入CNV信息 #############

library(Seurat)
library(Matrix)
library(data.table)

# 读取矩阵
expr <- readMM("cnv_matrix.mtx")

# 读取基因名
genes <- fread("cnv_genes.tsv", header = FALSE)
rownames(expr) <- genes$V1

# 读取元数据
meta <- fread("cnv_meta.tsv", header = FALSE)
colnames(meta) <- c("CB","orig.ident",'treatment', 'celltype', 'cnv_leiden', 'cnv_score',"cnv_status")

# 对齐 cell 名
colnames(expr) <- meta$CB
rownames(meta) <- meta$CB
meta$CB <- NULL

# 创建 Seurat 对象
seu <- CreateSeuratObject(
  counts = expr,
  meta.data = meta,
  assay = "RNA"
)

# 1️⃣ 设置默认身份
Idents(seu) <- "celltype"

# 2️⃣ 检查 assay
DefaultAssay(seu)

# 3️⃣ 看一下 metadata
head(seu@meta.data)


# 指定 treatment 顺序
treatment_order <- c("Adjacent normal", "0d-Ctrl", "14d-Si")

# 转换为有序因子
seu$treatment <- factor(
  seu$treatment,
  levels = treatment_order
)

# 同步到 meta.data（保险起见）
seu@meta.data$treatment <- seu$treatment


# 保存
saveRDS(seu, "cnv.rds")









###########绘制分群注释点图
library(ggh4x)

# 定义不同细胞类型的marker基因
cellmarker <- c(
  "CD68", "APOE",  "ITGAM",        # Macrophages
  "HLA-DPB1", "HLA-DRA", "HLA-DRB1", # Dendritic cells
  "COL1A1", "COL1A2", "DCN",                    # Fibroblasts / Stromal
  "GZMA", "GZMB",  "IFNG","CCL5",  # Cytotoxic T cells 细胞毒性T细胞
  "PDGFRB", "RGS5",   # Pericytes  周细胞
  "PECAM1", "VWF", "CDH5",             # Endothelial Cell
  #"KRT19", "MUC1", "KRT8","CDH1",         # Epithelial Cells /EOC上皮性卵巢癌
  "EPCAM","MUC16","WFDC2","KRT7","PAX8","CLDN4","CD24","WT1",      # Cancer Cells
  "FCN1","CD300E", "NLRP3" ,"TBC1D8",               # Monocytes 单核细胞
  "MZB1","JCHAIN",   # Plasma Cells 浆细胞
  "PCNA", "TOP2A", "CDK1"            # Proliferating cells 增殖细胞
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# 生成DotPlot
p <- DotPlot(scedata, features = unique(cellmarker))

# 获取数据
dat <- p$data

# 获取cluster注释并合并
anno <- distinct(data.frame(id = scedata$celltype, celltype = scedata$seurat_clusters))
colnames(anno) <- c("id", "celltype")
df <- left_join(dat, anno, by = "id")

# 定义cluster顺序
cluster.order <- c(0,2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)
df$celltype <- factor(df$celltype, levels = cluster.order)

# 绘制点图
p <- ggplot(df, aes(features.plot, interaction(celltype, id), size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank(),
    axis.text.y = element_text(color = 'black', size = 14),
    axis.text.x = element_text(color = 'black', size = 12, angle = 90, hjust = 1, vjust = 0.5)
  ) +
  scale_fill_gradientn(colours = c('#5749a0',  '#00bbb1', '#bef0b0', '#fdf4af', '#f9b64b', '#ec840e', '#ca443d', '#FF6666')) +
  guides(y = "axis_nested") +
  theme(
    ggh4x.axis.nesttext.y = element_text(colour = c('#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0', '#24796C', '#DAA51B', '#2F8AC4', '#764E9F', '#ED645A', '#CC3A8E')),
    ggh4x.axis.nestline.y = element_line(size = 2)
  )

# 保存图像
ggsave(plot = p, filename = paste0(output, "/ann_DotPlot.pdf"), width = 10, height = 8)
ggsave(plot = p, filename = paste0(output, "/ann_DotPlot.svg"), width = 10, height = 8)


#########绘制细胞注释热图############
# change annotation color
library("scales")
library(ggsci)
library(scRNAtoolVis)

# 设定注释颜色
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

pdf(file = paste(output, "ann_Heatmap.pdf",sep = '/'), width = 10, height = 8)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # 自定义高值颜色
dev.off()

svg(file = paste(output, "ann_Heatmap.svg",sep = '/'), width = 10, height = 8)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # 自定义高值颜色
dev.off()






#########marker基因在细胞中的表达趋势#######
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 挑选差异细胞展示
output <- paste(outdir,'cluster', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
#file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)


# 定义不同细胞类型的marker基因
cellmarker <- c(
  "CD68",       # Macrophages
  "HLA-DPB1",  # Dendritic cells
  "COL1A1",                  # Fibroblasts / Stromal
  "GZMA", # Cytotoxic T cells 细胞毒性T细胞
  "PDGFRB",  # Pericytes  周细胞
  "PECAM1",             # Endothelial Cell
  #"KRT19",# Epithelial Cells /EOC上皮性卵巢癌
  "MUC16","WFDC2",      # Cancer Cells
  "FCN1",            # Monocytes 单核细胞
  "MZB1"   # Plasma Cells 浆细胞
  
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# 创建存储所有 RidgePlot、VlnPlot 和 FeaturePlot 的列表
ridge_plots <- list()
vln_plots <- list()
feature_plots <- list()

# 遍历免疫细胞标志基因列表
for (gene in cellmarker) {
  # RidgePlot
  ridge_plots[[gene]] <- RidgePlot(ScRNA, features = gene, ncol = 1, cols = col) +
    theme(legend.position = "none")  
  
  # VlnPlot
  vln_plots[[gene]] <- VlnPlot(ScRNA, features = gene, ncol = 1, pt.size = 0, cols = col) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") 
  
  # 获取基因表达数据，确保返回为数值向量
  gene_expr <- as.numeric(FetchData(ScRNA, vars = gene, assay = "RNA")[[gene]])  # 转换为数值向量
  
  # 检查基因表达是否成功提取
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # 将基因表达数据加入元数据
  ScRNA@meta.data[[paste0(gene, "_expr")]] <- gene_expr  # 添加表达量信息到元数据
  
  # 使用数值向量进行排序，避免直接索引 ScRNA 对象
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # 提取排序后的细胞名称
  
  # FeaturePlot 按排序后的顺序绘制，并使用连续型颜色梯度
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", 
                                       cells = cell_names_ordered,  # 指定细胞顺序
                                       ncol = 1) + 
    scale_color_gradientn(colors = c("#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366")) +  # 设置连续型颜色梯度
    theme(legend.position = "right", 
          plot.title = element_text(size = 22, face = "bold"),  # 增大标题文字大小并加粗
          legend.text = element_text(size = 18)) +  # 增大图例文字大小
    NoAxes()  # 删除坐标轴
  
}


# 保存 RidgePlot 图
#pdf(paste0(out, "cellmarker_RidgePlot.pdf"), width = 25, height = 12)
#print(cowplot::plot_grid(plotlist = ridge_plots, ncol = 4))
#dev.off()

# 保存 FeaturePlot 图
pdf(paste0(output, "/marker_FeaturePlot_umap.pdf"), width = 25, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 5))
dev.off()

svg(paste0(output, "/marker_FeaturePlot_umap.svg"), width = 25, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 5))
dev.off()

pdf(paste0(output, "/marker_VlnPlot_umap.pdf"), width = 28, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =5))
dev.off()
svg(paste0(output, "/marker_VlnPlot_umap.svg"), width = 28, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =5))
dev.off()





#######################肿瘤标记物##########################
#########marker基因在细胞中的表达趋势#######
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 挑选差异细胞展示
output <- paste(outdir,'marker', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
#file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)


# 设置T细胞激活相关基因
cellmarker <- c(
  "EPCAM","WFDC2","KRT8","CLDN4","CD24",    # Cancer Cells  "TP53",
  "WT1","MUC16","PAX8",  #,"ESR1","PGR"   ## 卵巢
  "KRT7", "KRT20","CDX2"  ## 消化道  "MUC6"
  
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# 创建存储所有 RidgePlot、VlnPlot 和 FeaturePlot 的列表
ridge_plots <- list()
vln_plots <- list()
feature_plots <- list()

# 遍历免疫细胞标志基因列表
for (gene in cellmarker) {
  # RidgePlot
  ridge_plots[[gene]] <- RidgePlot(ScRNA, features = gene, ncol = 1, cols = col) +
    theme(legend.position = "none")  
  
  # VlnPlot
  vln_plots[[gene]] <- VlnPlot(ScRNA, features = gene, ncol = 1, pt.size = 0, cols = col) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") 
  
  # 获取基因表达数据，确保返回为数值向量
  gene_expr <- as.numeric(FetchData(ScRNA, vars = gene, assay = "RNA")[[gene]])  # 转换为数值向量
  
  # 检查基因表达是否成功提取
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # 将基因表达数据加入元数据
  ScRNA@meta.data[[paste0(gene, "_expr")]] <- gene_expr  # 添加表达量信息到元数据
  
  # 使用数值向量进行排序，避免直接索引 ScRNA 对象
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # 提取排序后的细胞名称
  
  # FeaturePlot 按排序后的顺序绘制，并使用连续型颜色梯度
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", pt.size=0.01,
                                       cells = cell_names_ordered,  # 指定细胞顺序
                                       ncol = 1) + 
    scale_color_gradientn(colors = c("#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366")) +  # 设置连续型颜色梯度
    #  "#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366"
    theme(legend.position = "right", 
          plot.title = element_text(size = 22, face = "bold"),  # 增大标题文字大小并加粗
          legend.text = element_text(size = 18)) +  # 增大图例文字大小
    NoAxes()  # 删除坐标轴
  
}


# 保存 RidgePlot 图
#pdf(paste0(out, "cellmarker_RidgePlot.pdf"), width = 25, height = 12)
#print(cowplot::plot_grid(plotlist = ridge_plots, ncol = 4))
#dev.off()

# 保存 FeaturePlot 图
pdf(paste0(output, "/spacial_FeaturePlot_umap.pdf"), width = 20, height = 12)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 4))
dev.off()


pdf(paste0(output, "/spacial_VlnPlot_umap.pdf"), width = 22, height = 10)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =4))
dev.off()





# 设置颜色（如果样本较多请酌情调整颜色列表）
col_sample <- c('#FF9999',"#A4CDE1","#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
                '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
                "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
                "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("0d", "3d","7d", "14d"))
# 创建存储新的 violin plots（按样本分组）
vln_plots_sample <- list()

for (gene in cellmarker) {
  # 绘制每个基因在不同样本中的小提琴图
  vln_plots_sample[[gene]] <- VlnPlot(
    ScRNA,
    features = gene,
    group.by = "treatment",     # <- 按样本分组
    pt.size = 0,
    cols = col_sample
  ) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle(paste(gene))
}

# 保存每个样本中基因表达的小提琴图
pdf(paste0(output, "/cellmarker_VlnPlot_bySample.pdf"), width = 12, height = 10)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 4))
dev.off()

svg(paste0(output, "/cellmarker_VlnPlot_bySample.svg"), width = 12, height = 10)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 4))
dev.off()



############绘制基因表达热图###########
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)

col <- c("#99CCFF","#FF3366","#66CCCC","#FF9933","#CC0066")

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
#ScRNA <- readRDS("ScRNA（分群后）.rds")

# 设置分组顺序
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("0d", "3d", "7d", "14d"))

# 设置待分析的T细胞激活相关基因
genes <- c(
  "EPCAM","WFDC2","KRT8","CLDN4","CD24",    # Cancer Cells  "TP53",
  "WT1","MUC16","PAX8",  #,"ESR1","PGR"   ## 卵巢
  "KRT7", "KRT20","CDX2"  ## 消化道  "MUC6"
  
)



# 提取表达矩阵并计算每组 treatment 的平均表达
avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment", assays = "RNA")$RNA

# 转换为矩阵
avg_expr_mat <- as.matrix(avg_expr)

# 确保行名和 gene 对应
rownames(avg_expr_mat) <- rownames(avg_expr)

# 按指定基因顺序排序
gene_order <- genes
#avg_expr_mat <- avg_expr_mat[gene_order, , drop = FALSE]

# Z-score 转换
mat_scaled <- t(scale(t(avg_expr_mat)))

# 设置分组信息
#group <- colnames(avg_expr_mat)  # 即 treatment 分组
#group_anno <- HeatmapAnnotation(
#  Treatment = factor(group),
#  #Treatment = factor(group, levels = c("0d", "3d", "7d", "14d")),
#  col = list(Treatment = c("WF-14W" = "#99CCFF","Tumor" =  "#FF3366","3d" = "#66CCCC", "7d" = "#FF9933", "14d" = "#CC0066")))

# 保存热图（PDF格式）
pdf(file.path(output, "cellmarker_Heatmap_byTreatment.pdf"), width = 5, height = 4)
Heatmap(mat_scaled,
        name = "Z-score",
        #top_annotation = group_anno,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 12),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = colorRampPalette(c("#47CFD1", "white", "#F04D28"))(100),
        heatmap_legend_param = list(
          title = "Z-score",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 10))
) %>% draw()
dev.off()




###############差异基因表达-点图##################
library(Seurat)
library(tidyverse)
library(ggsci)

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
#ScRNA <- readRDS("celltype.rds")


# 过滤在表达矩阵中实际存在的基因
genes <- genes[genes %in% rownames(ScRNA)]

# 提取treatment信息
ScRNA$treatment <- as.factor(ScRNA@meta.data$treatment)
treatment_groups <- levels(ScRNA$treatment)

# 获取表达矩阵
expr_matrix <- GetAssayData(ScRNA, slot = "data")[genes, ]


avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment")$RNA
#avg_expr_selected <- avg_expr[, treatment_groups]

# 保存平均表达表
avg_expr_df <- avg_expr %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Gene")
write.table(avg_expr_df, file = paste0(output, "/差异基因表达量.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)


# 绘制基因表达量点图（根据treatment分组）
plot <- DotPlot(ScRNA, features = unique(genes), group.by = "treatment") + 
  RotatedAxis() +
  coord_flip() +
  scale_color_gradientn(colors = c('#CCCCCC', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18))

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.pdf", sep='/'), plot = plot, width = 5, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.svg", sep='/'), plot = plot, width = 5, height = 5)




############绘制基因表达uamp###########
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)


output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
ScRNA <- readRDS("celltype.rds")

# 设置分组顺序
ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Ovarian cancer","Gastric cancer","Pancreatic cancer","WF-14W"))


#####12.展示已知的细胞marker基因的表达情况####
# 定义要绘制的基因列表
genes <- c(
  #"EPCAM","WFDC2","KRT8","CLDN4","CD24","MKI67","TOP2A",      # Cancer Cells  "TP53",
  "WT1","MUC16","PAX8",  #,"ESR1","PGR"   ## 卵巢
  "KRT7", "KRT20","CDX2"  ## 消化道  "MUC6"
)

#"MYO1E","MYC","KLF6","SAT1","JUN",

# 循环绘制并保存每个基因的特征图
for (gene in genes) {
  # 提取表达数据
  feature_data <- FetchData(ScRNA, vars = c(gene, "UMAP_1", "UMAP_2", "treatment"))
  
  # 按表达量排序
  feature_data <- feature_data %>% arrange(!!sym(gene))
  
  # 分组绘图
  p1 <- ggplot(feature_data, aes(x = UMAP_1, y = UMAP_2, color = !!sym(gene))) +
    geom_point(size = 0.2, alpha = 0.8) +
    scale_color_gradientn(colors = c('#E5D2DD',"#FF0066")) +
    theme_minimal() +
    facet_wrap(~treatment) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "right",
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    labs(title = paste0(gene), x = "UMAP_1", y = "UMAP_2", color = gene)
  
  ggsave(filename = file.path(output, paste0("FeaturePlot_", gene, "(不同组).pdf")), plot = p1, device = "pdf", width = 6, height = 3)
}






######### 消化道 #########
#########marker基因在细胞中的表达趋势#######
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 挑选差异细胞展示
output <- paste(outdir,'marker', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
#file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)


# 设置T细胞激活相关基因
cellmarker <- c(
  "VIL1","GCC1",  #"SPOCK1","APOC1","HER2",    ## 胃癌
  "CEA","CEACAM5","CEACAM6","VIL1","SPOCK1","APOC1","CXCR4","HER2",     ## 胃癌
  "ERCC1","KIAA1199","CEMIP","HERG1","KRT20",     ## 结直肠癌
  "FUT3","FUT6","KRT19","EGFR",    ## 胆囊癌 / 胆道癌
  "FUT3","CEACAM5","CEACAM6","SPAN1","TIMP1",        ## 胰腺癌
  "APOA4","TIMP1",            ## 胰腺癌
  "AFP","F2","CTNNB1"    # 肝癌
  
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# 创建存储所有 RidgePlot、VlnPlot 和 FeaturePlot 的列表
ridge_plots <- list()
vln_plots <- list()
feature_plots <- list()

# 遍历免疫细胞标志基因列表
for (gene in cellmarker) {
  # RidgePlot
  ridge_plots[[gene]] <- RidgePlot(ScRNA, features = gene, ncol = 1, cols = col) +
    theme(legend.position = "none")  
  
  # VlnPlot
  vln_plots[[gene]] <- VlnPlot(ScRNA, features = gene, ncol = 1, pt.size = 0, cols = col) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") 
  
  # 获取基因表达数据，确保返回为数值向量
  gene_expr <- as.numeric(FetchData(ScRNA, vars = gene, assay = "RNA")[[gene]])  # 转换为数值向量
  
  # 检查基因表达是否成功提取
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # 将基因表达数据加入元数据
  ScRNA@meta.data[[paste0(gene, "_expr")]] <- gene_expr  # 添加表达量信息到元数据
  
  # 使用数值向量进行排序，避免直接索引 ScRNA 对象
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # 提取排序后的细胞名称
  
  # FeaturePlot 按排序后的顺序绘制，并使用连续型颜色梯度
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", pt.size=0.01,
                                       cells = cell_names_ordered,  # 指定细胞顺序
                                       ncol = 1) + 
    scale_color_gradientn(colors = c("#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366")) +  # 设置连续型颜色梯度
    #  "#663399", "#3366CC", "#66CCCC", "#FFCC66", "#FF3366"
    theme(legend.position = "right", 
          plot.title = element_text(size = 22, face = "bold"),  # 增大标题文字大小并加粗
          legend.text = element_text(size = 18)) +  # 增大图例文字大小
    NoAxes()  # 删除坐标轴
  
}


# 保存 RidgePlot 图
#pdf(paste0(out, "cellmarker_RidgePlot.pdf"), width = 25, height = 12)
#print(cowplot::plot_grid(plotlist = ridge_plots, ncol = 4))
#dev.off()

# 保存 FeaturePlot 图
pdf(paste0(output, "/spacial_FeaturePlot_umap11.pdf"), width = 20, height = 20)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 4))
dev.off()



pdf(paste0(output, "/spacial_VlnPlot_umap11.pdf"), width = 22, height = 15)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =4))
dev.off()





# 设置颜色（如果样本较多请酌情调整颜色列表）
col_sample <- c('#FF9999',"#A4CDE1","#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
                '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
                "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
                "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("0d", "3d","7d", "14d"))
# 创建存储新的 violin plots（按样本分组）
vln_plots_sample <- list()

for (gene in cellmarker) {
  # 绘制每个基因在不同样本中的小提琴图
  vln_plots_sample[[gene]] <- VlnPlot(
    ScRNA,
    features = gene,
    group.by = "treatment",     # <- 按样本分组
    pt.size = 0,
    cols = col_sample
  ) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle(paste(gene))
}

# 保存每个样本中基因表达的小提琴图
pdf(paste0(output, "/cellmarker_VlnPlot_bySample11.pdf"), width = 12, height = 10)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 4))
dev.off()

svg(paste0(output, "/cellmarker_VlnPlot_bySample11.svg"), width = 12, height = 10)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 4))
dev.off()




############绘制基因表达热图###########
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)

col <- c("#99CCFF","#FF3366","#66CCCC","#FF9933","#CC0066")

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
#ScRNA <- readRDS("ScRNA（分群后）.rds")

# 设置分组顺序
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("0d", "3d", "7d", "14d"))

# 设置待分析的T细胞激活相关基因
genes <- c(
  
  "VIL1","GCC1",  #"SPOCK1","APOC1","HER2",    ## 胃癌
  "CEA","CEACAM5","CEACAM6","VIL1","SPOCK1","APOC1","CXCR4","HER2",     ## 胃癌
  "ERCC1","KIAA1199","CEMIP","HERG1","KRT20",     ## 结直肠癌
  "FUT3","FUT6","KRT19","EGFR",    ## 胆囊癌 / 胆道癌
  "FUT3","CEACAM5","CEACAM6","SPAN1","TIMP1",        ## 胰腺癌
  "APOA4","TIMP1",            ## 胰腺癌
  "AFP","F2","CTNNB1"    # 肝癌
)



# 提取表达矩阵并计算每组 treatment 的平均表达
avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment", assays = "RNA")$RNA

# 转换为矩阵
avg_expr_mat <- as.matrix(avg_expr)

# 确保行名和 gene 对应
rownames(avg_expr_mat) <- rownames(avg_expr)

# 按指定基因顺序排序
gene_order <- genes
#avg_expr_mat <- avg_expr_mat[gene_order, , drop = FALSE]

# Z-score 转换
mat_scaled <- t(scale(t(avg_expr_mat)))

# 设置分组信息
group <- colnames(avg_expr_mat)  # 即 treatment 分组
group_anno <- HeatmapAnnotation(
  Treatment = factor(group),
  #Treatment = factor(group, levels = c("0d", "3d", "7d", "14d")),
  col = list(Treatment = c("WF-14W" = "#99CCFF","Tumor" =  "#FF3366","3d" = "#66CCCC", "7d" = "#FF9933", "14d" = "#CC0066"))
)

# 保存热图（PDF格式）
pdf(file.path(output, "cellmarker_Heatmap_byTreatment11.pdf"), width = 6, height = 5)
Heatmap(mat_scaled,
        name = "Z-score",
        top_annotation = group_anno,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 12),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = colorRampPalette(c("#47CFD1", "white", "#F04D28"))(100),
        heatmap_legend_param = list(
          title = "Z-score",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 10))
) %>% draw()
dev.off()




###############差异基因表达-点图##################
library(Seurat)
library(tidyverse)
library(ggsci)

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
#ScRNA <- readRDS("celltype.rds")


# 过滤在表达矩阵中实际存在的基因
genes <- genes[genes %in% rownames(ScRNA)]

# 提取treatment信息
ScRNA$treatment <- as.factor(ScRNA@meta.data$treatment)
treatment_groups <- levels(ScRNA$treatment)

# 获取表达矩阵
expr_matrix <- GetAssayData(ScRNA, slot = "data")[genes, ]


avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment")$RNA
avg_expr_selected <- avg_expr[, treatment_groups]

# 保存平均表达表
avg_expr_df <- avg_expr_selected %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Gene")
write.table(avg_expr_df, file = paste0(output, "/差异基因表达量.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)


# 绘制基因表达量点图（根据treatment分组）
plot <- DotPlot(ScRNA, features = unique(genes), group.by = "treatment") + 
  RotatedAxis() +
  coord_flip() +
  scale_color_gradientn(colors = c('#CCCCCC', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18))

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_by_treatment11.pdf", sep='/'), plot = plot, width = 5, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_by_treatment11.svg", sep='/'), plot = plot, width = 5, height = 5)




#################绘制小提琴图################


#View(ScRNA@meta.data)

# 筛选存在于数据集中的marker基因
existing_markers <- cellmarker[cellmarker %in% rownames(ScRNA[["RNA"]]@data)]
existing_markers <- unique(existing_markers)

# 提取表达数据并转换为适合绘图的数据格式
vln.df <- as.data.frame(ScRNA[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

## 计算该基因全局 99th 分位数（在已去 0 的前提下），并去除 >99th 的“彪高”值
#vln.df <- vln.df %>% filter(!is.na(exp), exp > 0)
q05 <- quantile(vln.df$exp, 0.0001, na.rm = TRUE, names = FALSE)
q99 <- quantile(vln.df$exp, 0.99, na.rm = TRUE, names = FALSE)
vln.df <- vln.df %>% filter(exp > q05 & exp <= q99)

# 继续原有步骤
anno <- ScRNA@meta.data[, c("CB", "treatment")]
vln.df <- inner_join(vln.df, anno, by = "CB")

# 绘制Violin Plot，将X轴和Y轴调换
plot <- vln.df %>%
  ggplot(aes(exp, treatment)) +
  geom_violin(aes(fill = treatment), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # 调整facet_grid，以基因作为列
  scale_fill_manual(values = col_sample) +
  scale_x_continuous("") + scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 25),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    strip.text = element_text(size = 18, face = "bold"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# 保存Violin Plot
ggsave(filename = paste(output, "spacial_ViolinPlot.pdf", sep = '/'), plot = plot, width = 20,height=4,limitsize = FALSE)
ggsave(filename = paste(output, "spacial_ViolinPlot.svg", sep = '/'), plot = plot, width = 20,height=4,limitsize = FALSE)








####### 计算细胞比例 ###########

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
  # 深蓝→绿→浅绿 梯度
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq 浅色
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

# 设置分组顺序
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Ovarian cancer","Gastric cancer","Pancreatic cancer","WF-14W"))


######## 计算所有样本不同细胞群的细胞数
cell_counts <- as.data.frame(table(Idents(scedata)))
colnames(cell_counts) <- c("CellType", "Counts")

# 按细胞数从大到小排序
cell_counts <- cell_counts[order(-cell_counts$Counts), ]
# 保存所有细胞群的细胞数到指定目录
write.csv(cell_counts, paste(output, "cell_counts.csv", sep='/'), row.names = FALSE)

# 挑选前11种细胞保存到指定目录
cell_counts_top9 <- head(cell_counts, 11)
write.csv(cell_counts_top9, paste(output, "cell_counts_top9.csv", sep='/'), row.names = FALSE)

# 加载所需包
library(ggplot2)
# 绘制柱状图
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

# 保存图片为png格式
ggsave(paste(output, "cell_type_distribution.pdf", sep='/'), plot = p, width = 7, height = 6, dpi = 800)
ggsave(paste(output, "cell_type_distribution.svg", sep='/'), plot = p, width = 7, height = 6, dpi = 800)


# 计算各组不同细胞群的细胞数
# 按样本分组计算每种细胞类型的数量
cell_counts_group <- as.data.frame(table(scedata$orig.ident, Idents(scedata)))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# 添加分组信息 (假设分组变量为 `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # 确保分组信息的唯一性
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# 计算每个样本中每种细胞类型的比例
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
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # 修改X轴标签
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 20),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小
# 添加细胞数文本标签
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
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # 修改X轴标签
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 20),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小
# 添加细胞比例文本标签
p <- p + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 6)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 6, dpi = 800)



############分组############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# 计算每个处理组中每种细胞类型的比例
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

# 排序
#cell_counts_treatment$Treatment <- factor(cell_counts_treatment$Treatment, levels = c("Ovarian cancer","Gastric cancer","Pancreatic cancer","WF-14W"))

########## 绘制细胞数堆叠图 ##########
p1 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 22),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)

########## 绘制细胞比例堆叠图 ##########
p2 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 20),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小

file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 6, dpi = 800)





#################绘制细胞比例折线面积图###############
# 替换为你提供的颜色方案
cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
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
  # 深蓝→绿→浅绿 梯度
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq 浅色
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)

# 绘制数量折线面积图
p_counts_area <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, group = CellType)) +
  stat_summary(geom = 'line', fun = 'mean', color = 'white', linewidth = 1) +
  geom_area(aes(fill = CellType)) +
  scale_fill_manual(values = col) +
  labs(x = NULL, y = "Counts") +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.01, 0.01), name = "Counts") +  # ✅ 使用正确的 scale_y 并添加 Y 轴标题
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 16, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 16),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  geom_vline(aes(xintercept = "3d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "7d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "14d"), linetype = "dashed", size = 1, colour = "white")

ggsave(paste0(output, "/genecount_treatment_area.pdf"), plot = p_counts_area, width = 6, height = 5, dpi = 800)


# cell_counts_treatment$Ratio 应该已存在
p_ratio_area <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, group = CellType)) +
  stat_summary(geom = 'line', fun = 'mean', color = 'white', linewidth = 1) +
  geom_area(aes(fill = CellType)) +
  scale_fill_manual(values = col) +
  labs(x = NULL, y = "Ratio") +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.01, 0.01), name = "Ratio") +  # ✅ 使用正确的 scale_y 并添加 Y 轴标题
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 16, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 16),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  geom_vline(aes(xintercept = "3d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "7d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "14d"), linetype = "dashed", size = 1, colour = "white")

ggsave(paste0(output, "/geneRatio_treatment_area.pdf"), plot = p_ratio_area, width = 6, height = 5, dpi = 800)











####### 计算细胞比例 ###########
col <- c("#ca0020", "#FFCCCC", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1","#FF3366","#66CCCC","#9999FF",'#4F6272','#58A4C3',
         '#57C3F3', '#E59CC4','#437eb8',  "#FFCCCC","#CCFFCC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "cnv.rds")
scedata <- readRDS(file_path)

table(scedata$orig.ident)
view(scedata@meta.data)


############################################################
## 一、按 treatment 统计 copykat.pred：数量 & 比例并作图
############################################################

########## 1.1 统计 treatment × copykat.pred ##########
copy_treatment <- as.data.frame(table(scedata$orig.ident, scedata$cnv_status))
colnames(copy_treatment) <- c("Treatment", "cnv_status", "Counts")

# 计算每个 treatment 组内不同 copykat.pred 的比例
copy_treatment <- copy_treatment %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(Ratio = Counts / sum(Counts))

# 如需指定 treatment 的顺序可在此设置，例如：
# copy_treatment$Treatment <- factor(copy_treatment$Treatment,
#                                    levels = c("0d", "3d","7d", "14d"))


########## 1.2 按 treatment 绘制 copykat.pred 细胞数堆叠图 ##########
p_treat_counts <- ggplot(copy_treatment, aes(x = Treatment, y = Counts, fill = cnv_status)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))
p_treat_counts <- p_treat_counts + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- file.path(output, "copykatCount_treatment.pdf")
ggsave(file_path, plot = p_treat_counts,
       width = 3 * length(unique(copy_treatment$Treatment)), height = 6, dpi = 800)

file_path <- file.path(output, "copykatCount_treatment.svg")
ggsave(file_path, plot = p_treat_counts,
       width = 3 * length(unique(copy_treatment$Treatment)), height = 6, dpi = 800)


########## 1.3 按 treatment 绘制 copykat.pred 细胞比例堆叠图 ##########
p_treat_ratio <- ggplot(copy_treatment, aes(x = Treatment, y = Ratio, fill = cnv_status)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))
p_treat_ratio <- p_treat_ratio + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- file.path(output, "copykatRatio_treatment.pdf")
ggsave(file_path, plot = p_treat_ratio,
       width = 3 * length(unique(copy_treatment$Treatment)), height = 6, dpi = 800)

file_path <- file.path(output, "copykatRatio_treatment.svg")
ggsave(file_path, plot = p_treat_ratio,
       width = 3 * length(unique(copy_treatment$Treatment)), height = 6, dpi = 800)






############################################################
## 一、按 treatment 统计 copykat.pred：数量 & 比例并作图
############################################################

########## 1.1 统计 treatment × copykat.pred ##########
copy_treatment <- as.data.frame(table(scedata$celltype, scedata$cnv_status))
colnames(copy_treatment) <- c("celltype", "cnv_status", "Counts")

# 计算每个 treatment 组内不同 copykat.pred 的比例
copy_treatment <- copy_treatment %>%
  dplyr::group_by(celltype) %>%
  dplyr::mutate(Ratio = Counts / sum(Counts))

# 如需指定 treatment 的顺序可在此设置，例如：
# copy_treatment$Treatment <- factor(copy_treatment$Treatment,
#                                    levels = c("0d", "3d","7d", "14d"))


########## 1.2 按 treatment 绘制 copykat.pred 细胞数堆叠图 ##########
p_treat_counts <- ggplot(copy_treatment, aes(x = celltype, y = Counts, fill = cnv_status)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))
#p_treat_counts <- p_treat_counts + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- file.path(output, "copykatCount_celltype.pdf")
ggsave(file_path, plot = p_treat_counts,
       width = 1 * length(unique(copy_treatment$celltype)), height = 6, dpi = 800)

file_path <- file.path(output, "copykatCount_celltype.svg")
ggsave(file_path, plot = p_treat_counts,
       width = 1 * length(unique(copy_treatment$celltype)), height = 6, dpi = 800)


########## 1.3 按 treatment 绘制 copykat.pred 细胞比例堆叠图 ##########
p_treat_ratio <- ggplot(copy_treatment, aes(x = celltype, y = Ratio, fill = cnv_status)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))
#p_treat_ratio <- p_treat_ratio + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- file.path(output, "copykatRatio_celltype.pdf")
ggsave(file_path, plot = p_treat_ratio,
       width = 1 * length(unique(copy_treatment$celltype)), height = 6, dpi = 800)

file_path <- file.path(output, "copykatRatio_celltype.svg")
ggsave(file_path, plot = p_treat_ratio,
       width = 1 * length(unique(copy_treatment$celltype)), height = 6, dpi = 800)





#### 细胞类型的组织偏好性 ####

#BiocManager::install("impute", ask = FALSE, update = FALSE)
#devtools::install_github("Japrin/sscVis") #STARTRAC工具安装
#devtools::install_github("Japrin/Startrac") #STARTRAC工具安装

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1","#FF3366","#66CCCC","#9999FF",'#4F6272','#58A4C3',
         '#57C3F3', '#E59CC4','#437eb8',  "#FFCCCC","#CCFFCC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


library(Startrac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)


# 定义输出目录
output <- paste(outdir, 'STARTRAC_celltype_cnv', sep='/')
dir.create(output)

sco <- readRDS("cnv.rds")
data <- sco@meta.data
colnames(data)


#构建Roe计算需要的输入表格
#data <- data[,c(1,4,7,11)]
#colnames(data) <- c("sample","tissue","celltype")


# ---- 计算 Roe 并准备矩阵 ----
Roe <- calTissueDist(
  data,
  byPatient = FALSE,
  colname.cluster = "cnv_status",
  colname.patient = "orig.ident",
  colname.tissue = "orig.ident",
  method = "chisq",   # "chisq", "fisher", "freq"
  min.rowSum = 0
)

mat <- as.matrix(Roe)  # 后续统一用矩阵，避免数据框在 min/max 上的坑
rng <- range(mat, na.rm = TRUE)

# 若你希望以 1 为颜色中点（典型的 Ro/e 视觉语义），但 1 不在数据范围内，则改用范围中点
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng)

# 颜色映射：低值 -> 浅色，中点(1 或范围中点) -> 橙色，高值 -> 红色
col_fun <- circlize::colorRamp2(
  c(rng[1], mid, rng[2]),
  c("#f6f8e6", "#f9a33e", "red")
)

# 图例刻度：用 pretty() 生成“好看”的刻度；labels 与 at 一一对应
legend_at <- pretty(rng, n = 5)
legend_labels <- formatC(legend_at, format = "f", digits = 2)

# ---- 图1：带数值热图 ----
pdf(file.path(output, "STARTRAC_Roe_value.pdf"), width = 7, height = 3)
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
    col = col,  #指定行名的不同颜色向量
    fontface = "bold"      #可选，加粗字体增强显示
  ),
  column_names_gp = gpar(fontsize = 16),
  heatmap_legend_param = list(
    title = "Ro/e value",
    at = legend_at,
    labels = legend_labels,
    title_gp = gpar(fontsize = 16),  # 增大图例标题的字体大小
    labels_gp = gpar(fontsize = 15)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 16, col = "black"))
  }
)
dev.off()




#可视化2，自定义+++符号版本
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0


#可视化3，自定义+++符号及行细胞类型颜色版本
pdf(file.path(output, "STARTRAC_Roe_colored.pdf"), width = 7, height = 3)
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
          col = col,  #指定行名的不同颜色向量
          fontface = "bold"      #可选，加粗字体增强显示
        ),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(0, max(Roe)), 
          title_gp = gpar(fontsize = 16),  # 增大图例标题的字体大小
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
          # 优化文本对比度（根据背景色自动切换）
          text_color <- ifelse(mean(col2rgb(fill)) > 127, "black", "white")
          grid.text(symbol, x, y, gp = gpar(fontsize = 16, col = text_color))
        }
)
dev.off()




#可视化4，气泡图
roe_df <- as.data.frame(as.table(as.matrix(Roe))) %>% #将矩阵转换为数据框
  rename(Tissue = Var2, CellType = Var1, Roe = Freq) %>% #修改新建数据框列名
  mutate(
    Enrichment = ifelse(Roe >= 1, "Enrichment", "Depletion"),  #根据Roe值判断富集（大于1）或耗竭（小于1）
    Roe = ifelse(Roe == 0, NA, Roe)  #处理0值
  ) %>% 
  filter(!is.na(Roe))  #过滤无效值

pdf(file.path(output, "STARTRAC_Roe_bubble.pdf"), width = 5, height = 3)
ggplot(roe_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(
    aes(size = Roe, color = Enrichment)
  ) +
  scale_size_continuous(
    name = "Ro/e",
    breaks = c(0.5, 1.0, 1.5),
    range = c(1,7) #调整气泡大小范围
  ) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),  # 示例图配色
    labels = c("Enrichment", "Depletion"),
    guide = guide_legend(  
      override.aes = list(
        size = 4  #调整图例中点的大小
      )
    )
  ) +
  scale_y_discrete(limits = rev) +  # 保持y轴顺序与输入一致
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
    legend.box = "horizontal" #不同的图例相对位置,水平或垂直
  )
dev.off()










######### 通路分析--PROGENy #########
## ================== 环境准备 ==================
# 如未安装，请先安装 progeny 等包
# BiocManager::install("progeny")
# devtools::install_github("saezlab/progeny")

library(Seurat)
library(dplyr)
library(tidyr)
library(progeny)
library(pheatmap)
library(ggpubr)
library(ggsci)
library(tibble)

## 你自己的输出目录 outdir 需在脚本前面定义
## 假设 outdir 已有，这里只创建 progeny 子目录
output <- file.path(outdir, "progeny")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

## ================== 读取 Seurat 对象 ==================
file_path <- file.path(outdir, "celltype.rds")
sce <- readRDS(file_path)   # 用 sce 这个名字，与参考代码一致

## 设置分群为 celltype
Idents(sce) <- "treatment"

## 指定 CellType 顺序
#celltype_order <- c("Non-malignant cells", "Malignant cells")


## 做一个 cell 与 celltype 对照表（后面合并用）
CellsClusters <- data.frame(
  Cell     = colnames(sce),
  CellType = factor(as.character(Idents(sce))))



## ================== PROGENy通路活性分析 ==================
## 计算PROGENy通路活性分数并添加到Seurat对象
sce <- progeny(sce, scale = FALSE, organism = "Human", top = 500, 
               perm = 1, return_assay = TRUE)

## 对通路活性分数进行标准化
sce <- Seurat::ScaleData(sce, assay = "progeny")

## ================== 数据处理 ==================
## 将PROGENy分数转换为数据框便于处理
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(sce, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

## 将PROGENy分数与细胞类型信息合并
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## 按细胞类型和通路汇总活性分数
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(
    avg = mean(Activity), 
    std = sd(Activity),
    .groups = 'drop'  # 避免分组警告
  )

## ================== 热图数据准备 ==================
## 将数据转换为热图所需的矩阵格式
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  column_to_rownames("CellType") %>%
  as.matrix()

## ================== 颜色配置 ==================
paletteLength <- 100
myColor <- colorRampPalette(c('#57C3F3', "white", "#FF3366"))(paletteLength)

## 设置颜色断点，确保白色对应0值
progenyBreaks <- c(
  seq(min(summarized_progeny_scores_df), 0, 
      length.out = ceiling(paletteLength/2) + 1),
  seq(max(summarized_progeny_scores_df)/paletteLength, 
      max(summarized_progeny_scores_df), 
      length.out = floor(paletteLength/2))
)

## ================== 绘制热图 ==================
progeny_hmap <- pheatmap(
  t(summarized_progeny_scores_df),
  fontsize = 14,
  fontsize_row = 14,
  fontsize_col = 14,
  color = myColor, 
  breaks = progenyBreaks,
  main = "PROGENy",
  angle_col = "45",
  treeheight_col = 20,
  treeheight_row = 20,
  border_color = NA,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE
)

## ================== 保存热图 ==================
# 保存为PDF格式
pdf_file <- file.path(output, "PROGENy_heatmap.pdf")
pdf(pdf_file, width = 5, height = 5)
print(progeny_hmap)
dev.off()


## ================== 保存结果数据 ==================
# 保存汇总的通路活性分数
write.csv(summarized_progeny_scores_df, 
          file.path(output, "PROGENy_pathway_scores_by_celltype.csv"))

# 保存详细的通路活性数据
write.csv(progeny_scores_df, 
          file.path(output, "PROGENy_detailed_scores.csv"))




## ================== 绘制UMAP图 ==================
library(viridis)

## 设置默认assay为progeny
DefaultAssay(sce) <- 'progeny'

## 获取所有PROGENy通路名称
progeny_pathways <- rownames(GetAssayData(sce, assay = "progeny"))


# 绘制NFkB通路
p1 <- FeaturePlot(sce, features = "Androgen", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("Androgen") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

# 绘制MAPK通路
p2 <- FeaturePlot(sce, features = "TGFb", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis::turbo(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("TGFb") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

# 绘制MAPK通路
p3 <- FeaturePlot(sce, features = "EGFR", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis::turbo(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("EGFR") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

# 绘制MAPK通路
p4 <- FeaturePlot(sce, features = "EGFR", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis::turbo(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("EGFR") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

# 绘制MAPK通路
p5 <- FeaturePlot(sce, features = "EGFR", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis::turbo(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("EGFR") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))



# 保存单个通路的UMAP图
pdf(file.path(output, "PROGENy_UMAP_individual.pdf"), width = 5, height = 4)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()



# 单独保存Hypoxia通路UMAP图
pdf(file.path(output, "Androgen_pathway_UMAP.pdf"), width = 5, height = 4)
print(p1)
dev.off()

# 单独保存NFkB通路UMAP图
pdf(file.path(output, "TGFb_pathway_UMAP.pdf"), width = 5, height = 4)
print(p2)
dev.off()

# 单独保存TNFa通路UMAP图
pdf(file.path(output, "EGFR_pathway_UMAP.pdf"), width = 5, height = 4)
print(p3)
dev.off()

# 单独保存TNFa通路UMAP图
pdf(file.path(output, "Hypoxia_pathway_UMAP.pdf"), width = 5, height = 4)
print(p4)
dev.off()

# 单独保存TNFa通路UMAP图
pdf(file.path(output, "MAPK_pathway_UMAP.pdf"), width = 5, height = 4)
print(p4)
dev.off()





## ================== （可选）小提琴图：指定通路在不同细胞类型中的活性 ==================
# 设置需要比较的细胞类型及成对比较
# 这里以你给的例子为模板，可根据自己的细胞类型名称修改

col <- c("#DC050C",  "#1965B0","#FB8072", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A","#9999FF","#00CC66",
         "#A4CDE1",'#FF9999',"#66CCCC",'#4F6272',"#FF3366","#CC0066","#CC99CC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#FF6699","#6699CC","#FFFFCC")


library(dplyr)
library(ggpubr)
library(RColorBrewer)

# 选择要画的小提琴图通路
pathways_to_plot <- c("Androgen", "TGFb", "EGFR", "Hypoxia", "MAPK")

# 指定 CellType 顺序
#celltype_order <- c("Non-malignant cells", "Malignant cells")

# 生成颜色（与顺序一致）
my_colors <- c(
  "Non-malignant cells" = "#FFCCCC", 
  "Malignant cells"     = "#ca0020"
)

for (pw in pathways_to_plot) {
  
  pathway_df <- progeny_scores_df %>%
    filter(Pathway == pw) %>%
    mutate(
      Activity = as.numeric(Activity),
      CellType = factor(CellType)
    )
  
  p <- ggviolin(
    data   = pathway_df,
    x      = "CellType",
    y      = "Activity",
    fill   = "CellType",
    add    = "boxplot",
    add.params = list(fill = "white", width = 0.1)
  ) +
    scale_fill_manual(values = my_colors) +
    
    # 添加显著性 P 值
    stat_compare_means(
      comparisons = list(c("Non-malignant cells", "Malignant cells")),
      method = "wilcox.test",
      label = "p.format",
      size = 5
    ) +
    
    ggtitle(paste0(pw, " pathway activity")) +
    ylab("PROGENy activity score") +
    xlab("") +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 16,angle = 30, hjust = 1, vjust = 1),
      axis.text.y  = element_text(size = 16),
      axis.title.y = element_text(size = 18, face = "bold"),
      plot.title   = element_text(size = 20, hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(output, paste0("progeny_", pw, "_violin.pdf")),
    plot     = p,
    width    = 4,
    height   = 5.5
  )
}








############### AUC分析#############
# 安装并加载必要的包
#BiocManager::install("AUCell", force = TRUE)
#BiocManager::install("msigdbr", force = TRUE)
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(viridis)
library(msigdbr)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(Seurat)
col <- c('#FF6666','#E5D2DD',"#BC8F8F",'#FFCC99','#4F6272','#58A4C3',"#CC0066",
         '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  "#FF3366",'#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300", '#F9BB72', 
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 设置输出目录
output <- paste(outdir,"AUC分析", sep='/')
dir.create(output)

# 加载单细胞数据
file_path <- file.path(outdir, "celltype.rds")
ScRNA <- readRDS(file_path)

# 获取人类的KEGG基因集，并转为list格式
human_KEGG <- msigdbr(species = "Homo sapiens",
                      category = "C2",
                      subcategory = "KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)
human_KEGG_Set <- human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
print(unique(human_KEGG$gs_name))


# 获取小鼠的KEGG基因集，并转为list格式
#mouse_KEGG <- msigdbr(species = "Mus musculus",
#                      category = "C2",
#                      subcategory = "KEGG") %>% 
#  dplyr::select(gs_name, gene_symbol)
#mouse_KEGG_Set <- mouse_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
#print(unique(mouse_KEGG$gs_name))

# 选择多个通路进行分析
selected_pathways <- c("KEGG_IL_17_SIGNALING_PATHWAY","KEGG_TNF_SIGNALING_PATHWAY",
                       "KEGG_NF_KAPPA_B_SIGNALING_PATHWAY","KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                       "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY","KEGG_CHEMOKINE_SIGNALING_PATHWAY",
                       "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","KEGG_PROTEOGLYCANS_IN_CANCER",
                       "KEGG_PI3K_AKT_SIGNALING_PATHWAY","KEGG_FOCAL_ADHESION","KEGG_ECM_RECEPTOR_INTERACTION",
                       "KEGG_SMALL_CELL_LUNG_CANCER","KEGG_NON_SMALL_CELL_LUNG_CANCER",
                       "KEGG_MAPK_SIGNALING_PATHWAY","KEGG_TRANSCRIPTIONAL_MISREGULATION_IN_CANCER",
                       "KEGG_HIF_1_SIGNALING_PATHWAY","KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                       "KEGG_PD_L1_EXPRESSION_AND_PD_1_CHECKPOINT_PATHWAY_IN_CANCER",
                       "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                       "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY")

# 选择多个通路进行分析
selected_pathways <- c(
  ##炎症相关通路
  "KEGG_TNF_SIGNALING_PATHWAY","KEGG_NF_KAPPA_B_SIGNALING_PATHWAY","KEGG_IL_17_SIGNALING_PATHWAY","KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES","KEGG_CHEMOKINE_SIGNALING_PATHWAY","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_ARACHIDONIC_ACID_METABOLISM","KEGG_JAK_STAT_SIGNALING_PATHWAY",
  ##肺癌相关通路
  "KEGG_PATHWAYS_IN_CANCER","KEGG_SMALL_CELL_LUNG_CANCER","KEGG_NON_SMALL_CELL_LUNG_CANCER","KEGG_P53_SIGNALING_PATHWAY",
  "KEGG_PI3K_AKT_SIGNALING_PATHWAY","KEGG_HIPPO_SIGNALING_PATHWAY","KEGG_WNT_SIGNALING_PATHWAY","KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "KEGG_ERBB_SIGNALING_PATHWAY","KEGG_ECM_RECEPTOR_INTERACTION","KEGG_FOCAL_ADHESION",
  ##耐药性相关通路
  "KEGG_PI3K_AKT_SIGNALING_PATHWAY","KEGG_NF_KAPPA_B_SIGNALING_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY",
  "KEGG_MTOR_SIGNALING_PATHWAY","KEGG_WNT_SIGNALING_PATHWAY","KEGG_NOTCH_SIGNALING_PATHWAY",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY","KEGG_APOPTOSIS","KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
  ##糖酵解代谢相关通路
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS","KEGG_CITRATE_CYCLE_TCA_CYCLE","KEGG_PENTOSE_PHOSPHATE_PATHWAY",
  "KEGG_PYRUVATE_METABOLISM","KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM","KEGG_GALACTOSE_METABOLISM",
  "KEGG_STARCH_AND_SUCROSE_METABOLISM","KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
  ##免疫逃逸相关通路
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_PD_L1_EXPRESSION_AND_PD_1_CHECKPOINT_PATHWAY_IN_CANCER",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY","KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY","KEGG_TH1_AND_TH2_CELL_DIFFERENTIATION",
  "KEGG_TH17_CELL_DIFFERENTIATION","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
)

selected_features <- lapply(selected_pathways, function(pathway) {
  list(human_KEGG_Set[[pathway]])
})
names(selected_features) <- selected_pathways


# 创建空的metadata列以存储每个通路的AUC值
for (pathway in selected_pathways) {
  ScRNA[[pathway]] <- NA
}

# 对每个通路进行AUC分析
cells_rankings <- AUCell_buildRankings(ScRNA@assays$RNA@data, splitByBlocks = TRUE)

# 检查每个通路是否有对应的基因集
valid_pathways <- selected_pathways[sapply(selected_pathways, function(pathway) {
  !is.null(human_KEGG_Set[[pathway]]) && length(human_KEGG_Set[[pathway]]) > 0
})]
if (length(valid_pathways) == 0) {
  stop("No valid pathways found. Please check the selected_pathways or the KEGG data.")
}
print(paste("Valid pathways:", valid_pathways))

for (pathway in valid_pathways) {
  message("Processing pathway: ", pathway)
  
  gene_set <- list(human_KEGG_Set[[pathway]])
  names(gene_set) <- pathway
  
  cells_AUC <- AUCell_calcAUC(gene_set, cells_rankings, 
                              aucMaxRank = nrow(cells_rankings) * 0.1)
  
  ScRNA[[pathway]] <- as.numeric(getAUC(cells_AUC)[1, ])
}

################绘制小提琴图###########
library(ggpubr)
for (pathway in valid_pathways) {
  # 数据格式调整为适合ggviolin
  df <- ScRNA@meta.data
  df$AUC <- df[[paste0(pathway)]]
  p <- ggviolin(df, 
                x = "celltype", 
                y = "AUC", 
                width = 0.8, 
                color = "black", # 外框颜色
                fill = "celltype", # 填充颜色
                #palette = "npg", # 使用更美观的调色板
                add = "mean_sd", # 显示均值和标准差
                bxp.errorbar = TRUE, # 显示误差条
                bxp.errorbar.width = 0.5, # 调整误差条宽度
                size = 1, # 外框线条粗细
                outlier.shape = NA, # 不显示异常值点
                xlab = FALSE) +
    scale_fill_manual(values = col) + 
    labs(title = pathway, y = "AUC", x = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), # 删除网格线
          panel.border = element_rect(color = "black", fill = NA, size = 1), # 添加边框
          axis.text.x = element_text(angle = 45, hjust = 1, size = 20), # X轴文本大小
          axis.text.y = element_text(size = 20), # Y轴文本大小
          axis.title.y = element_text(size = 24), # Y轴标题大小
          legend.title = element_text(size = 22), # 图例标题大小
          legend.text = element_text(size = 20), # 图例文本大小
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold"),
          legend.position = "none")
  
  ggsave(filename = paste0(output, "/violin_", pathway, ".pdf"),plot = p, width = 10,height = 6)
  ggsave(filename = paste0(output, "/violin_", pathway, ".svg"),plot = p, width = 10,height = 6)
}


# 提取UMAP坐标数据
umap <- data.frame(ScRNA@meta.data, ScRNA@reductions$umap@cell.embeddings)

# 绘制UMAP图
for (pathway in valid_pathways) {
  geneSet <- paste0(pathway)
  # 将 pathway 对应的 AUCell 数据添加到 umap 数据框中
  umap$AUCell <- umap[[geneSet]]
  # 按表达量从低到高排序
  umap <- umap[order(umap$AUCell, na.last = TRUE), ]
  cell_type_med <- umap %>%
    group_by(celltype) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2))
  p <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = AUCell)) +
    scale_color_gradientn(
      colors = c("#333399", "#66CCCC", "#FFCC66", "#FF3366"), # 使用白色到粉红的渐变
      limits = c(min(umap$AUCell, na.rm = TRUE), max(umap$AUCell, na.rm = TRUE))) +
    ggrepel::geom_label_repel(
      aes(label = celltype),
      fontface = "bold",
      data = cell_type_med,
      point.padding = unit(0.5, "lines"),
      size = 6) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),       # 删除网格线
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_blank(), # 删除背景
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), # 增大标题
      axis.title.x = element_text(size = 22, face = "bold"),  # 增大X轴标题文本
      axis.title.y = element_text(size = 22, face = "bold"),  # 增大Y轴标题文本
      axis.text.x = element_text(size = 20),   # 增大X轴刻度文本
      axis.text.y = element_text(size = 20),   # 增大Y轴刻度文本
      legend.title = element_text(size = 22, face = "bold"), # 增大图例标题文本
      legend.text = element_text(size = 20),   # 增大图例文本
      legend.key.size = unit(1, "lines")     # 调整图例键的大小
    ) +
    labs(title = pathway, color = "AUC")
  
  ggsave(filename = paste0(output, "/umap_", pathway, ".pdf"), plot = p, width = 8, height = 6)
  ggsave(filename = paste0(output, "/umap_", pathway, ".svg"), plot = p, width = 8, height = 6)
}








########## 代谢分析 ###############

#install.packages(c("devtools", "data.table", "wesanderson", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
#如果有些包安装不上，尝试使用BiocManager::install进行安装，比如GSVA包
#BiocManager::install("GSVA")
#devtools::install_github("cran/loe")
#devtools::install_github("YosefLab/VISION@v2.1.0") #Please note that the version would be v2.1.0
#安装scMetabolism
#devtools::install_github("wu-yc/scMetabolism")

library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(gridExtra)
# 安装wesanderson包
#install.packages("wesanderson")
library(wesanderson)

col <- c('#FF9999','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'代谢分析', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)
summary(scRNAsub$seurat_clusters)


# 设置细胞类型为身份信息
Idents(scRNAsub) <- "treatment"

# 筛选"Fib_cisplatin_res"和"Fib_other"细胞类型
#fib_subset <- subset(scRNAsub, idents = c("Fib_cisplatin_res", "Fib_other"))

# 运行scMetabolism分析
countexp.Seurat <- sc.metabolism.Seurat(obj = scRNAsub,  # Seurat对象
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 2, 
                                        metabolism.type = "KEGG")

# 提取score结果
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4, 1:4]

# 将score中barcode的点转为下划线
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))

# 验证barcode一致性
identical(colnames(score_change), rownames(countexp.Seurat@meta.data))

# 合并元数据和代谢评分
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data, t(score_change))
#View(countexp.Seurat@meta.data)
print(rownames(score))

# 展示前10个代谢通路
input.pathway <- rownames(countexp.Seurat@assays$METABOLISM$score)[1:30]
# 指定分析的代谢通路
input.pathway <- c(
  #铂耐药核心通路：
  "Glutathione metabolism", "Purine metabolism", "Pyrimidine metabolism", "Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Fatty acid biosynthesis", "Fatty acid degradation",
  
  #早期到晚期进展核心通路：
  "Glycolysis / Gluconeogenesis", "Mucin type O-glycan biosynthesis", "Glycosaminoglycan degradation", "Arachidonic acid metabolism",
  
  "Glycolysis / Gluconeogenesis",
  "Citrate cycle (TCA cycle)",
  "Pentose phosphate pathway",
  "Pyruvate metabolism",
  "Oxidative phosphorylation",
  "Alanine, aspartate and glutamate metabolism",
  "Cysteine and methionine metabolism",
  "Glutathione metabolism",
  "Fatty acid biosynthesis",
  "Fatty acid elongation",
  "Fatty acid degradation",
  "Glycerolipid metabolism",
  "Glycerophospholipid metabolism",
  "Sphingolipid metabolism",
  "Arachidonic acid metabolism",
  "Steroid biosynthesis",
  "Steroid hormone biosynthesis",
  "Purine metabolism",
  "Pyrimidine metabolism",
  "Metabolism of xenobiotics by cytochrome P450",
  "Drug metabolism - cytochrome P450",
  "Drug metabolism - other enzymes"
  
)

#"Glycolysis / Gluconeogenesis","Drug metabolism - cytochrome P450","Metabolism of xenobiotics by cytochrome P450","Fatty acid biosynthesis", 
#"Fructose and mannose metabolism","Arachidonic acid metabolism",

view(countexp.Seurat@meta.data)




# 绘制DotPlot
pdf(file.path(output, "Metabolism_DotPlot.pdf"), width = 6, height = 7)
DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "treatment",  
                   norm = "y") +
  scale_color_gradient(low = "#66CCCC", high = "#FF3366") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14))
dev.off()

svg(file.path(output,"Metabolism_DotPlot.svg"), width = 6, height = 7)
DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "treatment",  
                   norm = "y") +
  scale_color_gradient(low = "#66CCCC", high = "#FF3366") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14))
dev.off()










input.pathway <- c(   
  
  #铂耐药核心通路：
  "Glutathione metabolism", "Purine metabolism", "Pyrimidine metabolism", "Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Fatty acid biosynthesis", "Fatty acid degradation",
  
  #早期到晚期进展核心通路：
  "Glycolysis / Gluconeogenesis", "Mucin type O-glycan biosynthesis", "Glycosaminoglycan degradation", "Arachidonic acid metabolism"
)

# 循环遍历每个通路并生成FeaturePlot和VlnPlot组合图
for (pathway in input.pathway) {
  # 生成FeaturePlot
  p1 <- FeaturePlot(countexp.Seurat, reduction = "umap", features = pathway,split.by = 'treatment', 
                    cols = c("#660066","#00CC99","#FFFF99"))
  
  # 生成VlnPlot
  p2 <- VlnPlot(countexp.Seurat, features = pathway, cols = col, pt.size = 0)
  
  # 组合并显示
  combined_plot <- p1 + p2
  print(combined_plot)
  
  # 保存图像到文件 (可选)
  #ggsave(filename = paste0(pathway, "_Feature_VlnPlot.pdf"), plot = combined_plot, width = 5, height = 6)
  #ggsave(filename = paste0(pathway, "_Feature_VlnPlot.svg"), plot = combined_plot, width = 5, height = 6)
  
}




input.pathway <- c( 
  #铂耐药核心通路：
  "Glutathione metabolism","Pyrimidine metabolism", "Oxidative phosphorylation", 
  
  #早期到晚期进展核心通路：
  "Arachidonic acid metabolism","Glycolysis / Gluconeogenesis"
  
)


# 绘制代谢通路的tSNE图
plot_list <- list()
for (pathway in input.pathway) {
  p <- DimPlot.metabolism(obj = countexp.Seurat, 
                          pathway = pathway, 
                          dimention.reduction.type = "umap",  
                          dimention.reduction.run = F, size = 1) +
    ggtitle(pathway)+
    theme(
      legend.title = element_blank(),          # 删除图例标题
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # 设置主标题大小、加粗、居中
      axis.title = element_text(size = 14),    # 设置轴标题大小
      axis.text = element_text(size = 14),     # 设置坐标轴刻度大小
      legend.text = element_text(size = 14)    # 设置图例文本大小
    )
  plot_list[[pathway]] <- p
}
pdf(file.path(output, "Metabolism_combined.pdf"), width = 25, height = 4)
grid.arrange(grobs = plot_list, ncol = 5)
dev.off()
svg(file.path(output, "Metabolism_combined.svg"), width = 25, height = 4)
grid.arrange(grobs = plot_list, ncol = 5)
dev.off()








# 绘制代谢通路的箱线图
boxplot_list <- list()
for (pathway in input.pathway) {
  p <- BoxPlot.metabolism(obj = countexp.Seurat, 
                          pathway = pathway, 
                          phenotype = "stage", 
                          ncol = 1) +                          
    scale_fill_manual(values =col)+                                        
    theme(
      axis.title = element_text(size = 14),                    
      axis.text = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 13, face = "bold"))+
    stat_compare_means(method = "wilcox.test",label = "p.format", label.x.npc = 'middle',label.x = 1.5, size = 5) # 添加p值
  boxplot_list[[pathway]] <- p
}

pdf(file.path(output,"Metabolism_BoxPlot_combined.pdf"), width = 20, height = 5)
grid.arrange(grobs = boxplot_list, ncol = 5)
dev.off()
svg(file.path(output,"Metabolism_BoxPlot_combined.svg"), width = 20, height = 5)
grid.arrange(grobs = boxplot_list, ncol = 5)
dev.off()








# 方法1：使用FeaturePlot按treatment分组显示
plot_list_feature <- list()
for (pathway in input.pathway) {
  p <- FeaturePlot(countexp.Seurat, 
                   features = pathway,
                   split.by = "treatment",  # 按treatment分组
                   reduction = "umap",
                   cols = c("lightgrey", "red"),  # 设置颜色
                   pt.size = 1) +
    ggtitle(pathway) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 12)  # 设置分组标签文字大小
    )
  plot_list_feature[[pathway]] <- p
}

# 保存FeaturePlot结果
pdf(file.path(output, "Metabolism_by_treatment_FeaturePlot.pdf"), width = 30, height = 8)
for (i in seq_along(plot_list_feature)) {
  print(plot_list_feature[[i]])
}
dev.off()

svg(file.path(output, "Metabolism_by_treatment_FeaturePlot.svg"), width = 30, height = 8)
for (i in seq_along(plot_list_feature)) {
  print(plot_list_feature[[i]])
}
dev.off()

# 方法2：使用VlnPlot按treatment分组显示代谢通路活性
plot_list_vln <- list()
for (pathway in input.pathway) {
  p <- VlnPlot(countexp.Seurat, 
               features = pathway,
               group.by = "treatment",  # 按treatment分组
               pt.size = 0) +  # pt.size = 0 不显示点
    ggtitle(pathway) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10)
    )
  plot_list_vln[[pathway]] <- p
}

# 保存VlnPlot结果
pdf(file.path(output, "Metabolism_by_treatment_Violin.pdf"), width = 30, height = 8)
grid.arrange(grobs = plot_list_vln, ncol = 4)
dev.off()

svg(file.path(output, "Metabolism_by_treatment_Violin.svg"), width = 30, height = 8)
grid.arrange(grobs = plot_list_vln, ncol = 4)
dev.off()

# 方法3：分别绘制每个treatment组的UMAP图
treatment_groups <- unique(countexp.Seurat@meta.data$treatment)

# 为每个代谢通路创建按treatment分组的图
for (pathway in input.pathway) {
  plot_list_treatment <- list()
  
  for (treatment_group in treatment_groups) {
    # 创建当前treatment组的子集
    cells_to_keep <- WhichCells(countexp.Seurat, expression = treatment == treatment_group)
    subset_obj <- subset(countexp.Seurat, cells = cells_to_keep)
    
    p <- FeaturePlot(subset_obj, 
                     features = pathway,
                     reduction = "umap",
                     cols = c("lightgrey", "red"),
                     pt.size = 1) +
      ggtitle(paste(pathway, "-", treatment_group)) +
      theme(
        legend.title = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)
      )
    plot_list_treatment[[treatment_group]] <- p
  }
  
  # 保存当前代谢通路的所有treatment组图
  pdf(file.path(output, paste0("Metabolism_", gsub("/", "_", gsub(" ", "_", pathway)), "_by_treatment.pdf")), 
      width = 8 * length(treatment_groups), height = 8)
  grid.arrange(grobs = plot_list_treatment, ncol = length(treatment_groups))
  dev.off()
}






#BiocManager::install("impute", ask = FALSE, update = FALSE)
#devtools::install_github("Japrin/sscVis") #STARTRAC工具安装
#devtools::install_github("Japrin/Startrac") #STARTRAC工具安装

library(Startrac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)

output <- paste(outdir,'组织偏好性', sep='/')
dir.create(output)


sco <- readRDS("celltype.rds")
data <- sco@meta.data
colnames(data)


#构建Roe计算需要的输入表格
#data <- data[,c(1,4,7,11)]
#colnames(data) <- c("sample","tissue","celltype")


# ---- 计算 Roe 并准备矩阵 ----
Roe <- calTissueDist(
  data,
  byPatient = FALSE,
  colname.cluster = "celltype",
  colname.patient = "orig.ident",
  colname.tissue = "treatment",
  method = "chisq",   # "chisq", "fisher", "freq"
  min.rowSum = 0
)

mat <- as.matrix(Roe)  # 后续统一用矩阵，避免数据框在 min/max 上的坑
rng <- range(mat, na.rm = TRUE)

# 若你希望以 1 为颜色中点（典型的 Ro/e 视觉语义），但 1 不在数据范围内，则改用范围中点
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng)

# 颜色映射：低值 -> 浅色，中点(1 或范围中点) -> 橙色，高值 -> 红色
col_fun <- circlize::colorRamp2(
  c(rng[1], mid, rng[2]),
  c("#f6f8e6", "#f9a33e", "red")
)

# 图例刻度：用 pretty() 生成“好看”的刻度；labels 与 at 一一对应
legend_at <- pretty(rng, n = 5)
legend_labels <- formatC(legend_at, format = "f", digits = 2)

# ---- 图1：带数值热图 ----
pdf(file.path(output, "STARTRAC_Roe_value.pdf"), width = 6, height = 6)
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
    col = col,  #指定行名的不同颜色向量
    fontface = "bold"      #可选，加粗字体增强显示
  ),
  column_names_gp = gpar(fontsize = 16),
  heatmap_legend_param = list(
    title = "Ro/e value",
    at = legend_at,
    labels = legend_labels,
    title_gp = gpar(fontsize = 16),  # 增大图例标题的字体大小
    labels_gp = gpar(fontsize = 15)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 16, col = "black"))
  }
)
dev.off()




#可视化2，自定义+++符号版本
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0

pdf(file.path(output, "STARTRAC_Roe.pdf"), width = 6, height = 6)
Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = 'right',
        column_names_side = "bottom",
        show_column_names = TRUE,
        show_row_names = TRUE,
        #row_split = 3,
        col = col_fun,
        row_names_gp = gpar(
          fontsize = 16,
          col = col,  #指定行名的不同颜色向量
          fontface = "bold"      #可选，加粗字体增强显示
        ),
        column_names_gp = gpar(fontsize = 16),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(0, max(Roe)),
          title_gp = gpar(fontsize = 16),  # 增大图例标题的字体大小
          labels_gp = gpar(fontsize = 15),
          labels = c("0", "Max.")  # 对应标签
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
          # 显示符号
          grid.text(symbol, x, y, gp = gpar(
            fontsize = 16, 
            col = "black")
          )
        }
)
dev.off()



#可视化3，自定义+++符号及行细胞类型颜色版本
pdf(file.path(output, "STARTRAC_Roe_colored.pdf"), width = 6, height = 6)
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
          col = col,  #指定行名的不同颜色向量
          fontface = "bold"      #可选，加粗字体增强显示
        ),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(0, max(Roe)), 
          title_gp = gpar(fontsize = 16),  # 增大图例标题的字体大小
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
          # 优化文本对比度（根据背景色自动切换）
          text_color <- ifelse(mean(col2rgb(fill)) > 127, "black", "white")
          grid.text(symbol, x, y, gp = gpar(fontsize = 16, col = text_color))
        }
)
dev.off()




#可视化4，气泡图
roe_df <- as.data.frame(as.table(as.matrix(Roe))) %>% #将矩阵转换为数据框
  rename(Tissue = Var2, CellType = Var1, Roe = Freq) %>% #修改新建数据框列名
  mutate(
    Enrichment = ifelse(Roe >= 1, "Enrichment", "Depletion"),  #根据Roe值判断富集（大于1）或耗竭（小于1）
    Roe = ifelse(Roe == 0, NA, Roe)  #处理0值
  ) %>% 
  filter(!is.na(Roe))  #过滤无效值

pdf(file.path(output, "STARTRAC_Roe_bubble.pdf"), width = 9, height = 6)
ggplot(roe_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(
    aes(size = Roe, color = Enrichment)
  ) +
  scale_size_continuous(
    name = "Ro/e",
    breaks = c(0.5, 1.0, 1.5),
    range = c(1,7) #调整气泡大小范围
  ) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),  # 示例图配色
    labels = c("Enrichment", "Depletion"),
    guide = guide_legend(  
      override.aes = list(
        size = 4  #调整图例中点的大小
      )
    )
  ) +
  scale_y_discrete(limits = rev) +  # 保持y轴顺序与输入一致
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
    legend.box = "horizontal" #不同的图例相对位置,水平或垂直
  )
dev.off()







#### 划分转移和非转移细胞 ########
#BiocManager::install("UCell")
library(UCell)
library(Seurat)
#install.packages("viridis")
library(viridis)
library(ggplot2)
library(stringr)

# 设置输出目录
output <- paste(outdir, 'Ucell', sep='/')
dir.create(output, showWarnings = FALSE)

# 1. 读取 Seurat 数据
sce <- readRDS("celltype.rds")

# 2. 定义基因集
geneSet <- c(
  ## 胰腺癌转移相关基因
  "SMAD4", "KRAS", "TP53", "CDKN2A", "SNAI1", "SNAI2", "TWIST1", 
  "ZEB1", "ZEB2", "FOSL1", "ITGA2", "ITGB1", "CDH1", "MMP2", "MMP9", 
  "RHOA", "RAC1", "CDC42", "CXCR4", "CXCL12","ACKR3", 
  "TGFB1", "TGFBR1", "TGFBR2", "AXL", "MET", "EGFR", "MYC", "HIF1A", "VEGFA"
)  

# 3. 计算 UCell 评分
sce <- AddModuleScore_UCell(sce, features = list(Meta = geneSet), name = "_UCell")
score_col  <- colnames(sce@meta.data) %>% str_subset("UCell")


# 取出分数向量
scores <- sce@meta.data[[score_col]]

# 4A. 分成 16 个等分位箱（labels 为 B1..B16）
# 使用分位数切分，避免数值分布不均带来的空箱（去重以防重复分位点）
probs <- seq(0, 1, length.out = 17)
brks_raw <- quantile(scores, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
brks <- unique(brks_raw)

# 如果因为大量并列值导致分位点去重后箱数 < 16，则退而使用 pretty 切分保证箱数
if (length(brks) < 17) {
  brks <- pretty(range(scores, na.rm = TRUE), n = 16)
  brks <- unique(brks)
}

# 再次确保箱数至少为 2（极端情况兜底）
if (length(brks) < 2) {
  brks <- range(scores, na.rm = TRUE)
  brks[1] <- brks[1] - 1e-8
  brks[2] <- brks[2] + 1e-8
}

# 生成 16-bin（若有效断点不足 17，bin 数会略少于 16，这是数据本身导致的）
n_bins <- length(brks) - 1
bin_labels <- paste0("B", seq_len(n_bins))
sce$Meta_bin16 <- cut(scores, breaks = brks, include.lowest = TRUE, labels = bin_labels, right = TRUE, ordered_result = TRUE)


# 4B. 20/60/20 分级（High / Medium / Low）
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

# 5. 可视化
# 如果没有 UMAP/TSNE，请先运行：sce <- RunPCA(sce); sce <- RunUMAP(sce, dims = 1:30)
# 5A. 16-bin 离散可视化（DimPlot 更适合离散型）
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


# 5B. High/Medium/Low 可视化
# 给三类指定清晰的离散色板
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

# 7. 如仍想保留连续分数的 FeaturePlot（原图），可选保留：
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


### 分组 ###
# 7. 如仍想保留连续分数的 FeaturePlot（原图），可选保留：
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
#devtools::install_github("Japrin/sscVis") #STARTRAC工具安装
#devtools::install_github("Japrin/Startrac") #STARTRAC工具安装

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


#构建Roe计算需要的输入表格
#data <- data[,c(1,4,7,11)]
#colnames(data) <- c("sample","tissue","celltype")


# ---- 计算 Roe 并准备矩阵 ----
Roe <- calTissueDist(
  data,
  byPatient = FALSE,
  colname.cluster = "Meta_3class",
  colname.patient = "orig.ident",
  colname.tissue = "treatment",
  method = "chisq",   # "chisq", "fisher", "freq"
  min.rowSum = 0
)

mat <- as.matrix(Roe)  # 后续统一用矩阵，避免数据框在 min/max 上的坑
rng <- range(mat, na.rm = TRUE)

# 若你希望以 1 为颜色中点（典型的 Ro/e 视觉语义），但 1 不在数据范围内，则改用范围中点
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng)

# 颜色映射：低值 -> 浅色，中点(1 或范围中点) -> 橙色，高值 -> 红色
col_fun <- circlize::colorRamp2(
  c(rng[1], mid, rng[2]),
  c("#f6f8e6", "#f9a33e", "red")
)

# 图例刻度：用 pretty() 生成“好看”的刻度；labels 与 at 一一对应
legend_at <- pretty(rng, n = 5)
legend_labels <- formatC(legend_at, format = "f", digits = 2)

# ---- 图1：带数值热图 ----
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
    col = col,  #指定行名的不同颜色向量
    fontface = "bold"      #可选，加粗字体增强显示
  ),
  column_names_gp = gpar(fontsize = 16),
  heatmap_legend_param = list(
    title = "Ro/e value",
    at = legend_at,
    labels = legend_labels,
    title_gp = gpar(fontsize = 16),  # 增大图例标题的字体大小
    labels_gp = gpar(fontsize = 15)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 16, col = "black"))
  }
)
dev.off()




#可视化2，自定义+++符号版本
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0


#可视化3，自定义+++符号及行细胞类型颜色版本
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
          col = col,  #指定行名的不同颜色向量
          fontface = "bold"      #可选，加粗字体增强显示
        ),
        heatmap_legend_param = list(
          title = "Ro/e",
          at = c(0, max(Roe)), 
          title_gp = gpar(fontsize = 16),  # 增大图例标题的字体大小
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
          # 优化文本对比度（根据背景色自动切换）
          text_color <- ifelse(mean(col2rgb(fill)) > 127, "black", "white")
          grid.text(symbol, x, y, gp = gpar(fontsize = 16, col = text_color))
        }
)
dev.off()




#可视化4，气泡图
roe_df <- as.data.frame(as.table(as.matrix(Roe))) %>% #将矩阵转换为数据框
  rename(Tissue = Var2, CellType = Var1, Roe = Freq) %>% #修改新建数据框列名
  mutate(
    Enrichment = ifelse(Roe >= 1, "Enrichment", "Depletion"),  #根据Roe值判断富集（大于1）或耗竭（小于1）
    Roe = ifelse(Roe == 0, NA, Roe)  #处理0值
  ) %>% 
  filter(!is.na(Roe))  #过滤无效值

pdf(file.path(output, "STARTRAC_Roe_bubble.pdf"), width = 5, height = 3)
ggplot(roe_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(
    aes(size = Roe, color = Enrichment)
  ) +
  scale_size_continuous(
    name = "Ro/e",
    breaks = c(0.5, 1.0, 1.5),
    range = c(1,7) #调整气泡大小范围
  ) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),  # 示例图配色
    labels = c("Enrichment", "Depletion"),
    guide = guide_legend(  
      override.aes = list(
        size = 4  #调整图例中点的大小
      )
    )
  ) +
  scale_y_discrete(limits = rev) +  # 保持y轴顺序与输入一致
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
    legend.box = "horizontal" #不同的图例相对位置,水平或垂直
  )
dev.off()









#####19.比较两组间的差异分析#####
library(scRNAtoolVis)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggrepel)

col <- c('#437eb8','#FF6666',"#FFFFCC",'#FFCC99','#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFCCCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")

# 选取前6细胞群展示标记基因的表达趋势
output <- paste(outdir,'差异分析(转移vs非转移)', sep='/')
dir.create(output, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)
colnames(scRNAsub@meta.data)


# 寻找 Res 和 Sen 组之间的差异基因
logFCfilter <- 0.25        # 定义 log2FC 过滤值
adjPvalFilter <- 0.05   # 定义矫正后 P 值过滤值

# 寻找 Epi_cisplatin_res 和 Epi_other 组之间的差异基因
scRNAsub.cluster.markers <- FindMarkers(object = scRNAsub, 
                                        ident.1 = "Met_high",
                                        ident.2 =  "Met_low",
                                        group.by = "Meta_3class", 
                                        logfc.threshold = 0, 
                                        min.pct = 0.25, 
                                        test.use = "wilcox")
scRNAsub.cluster.markers$gene <- rownames(scRNAsub.cluster.markers)

# 添加显著性标注
scRNAsub.cluster.markers <- scRNAsub.cluster.markers %>%
  mutate(Significance = ifelse(p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter, 
                               ifelse(avg_log2FC > 0, "Up", "Down"), "Normal"))
write.table(scRNAsub.cluster.markers, file = file.path(output,"sig.markers_ann_Tumor_vs_Normal.txt"), sep = "\t",row.names = T, quote = FALSE)

saveRDS(scRNAsub.cluster.markers, file = file.path(output, "ScRNA.sig.markers.rds"))

# 分别保存上调基因和下调基因
upregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")
downregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")
write.csv(upregulated_genes, file = file.path(output, "upregulated_genes_Tumor_vs_Normal.csv"), row.names = TRUE, quote = FALSE)
write.csv(downregulated_genes, file = file.path(output, "downregulated_genes_Tumor_vs_Normal.csv"), row.names = TRUE, quote = FALSE)

# 计算上调和下调基因数目
upregulated_genes <- sum(scRNAsub.cluster.markers$Significance == "Up")
downregulated_genes <- sum(scRNAsub.cluster.markers$Significance == "Down")
total_diff_genes <- upregulated_genes + downregulated_genes

# 分别保存上调基因和下调基因的数据框
upregulated_genes_df <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")
downregulated_genes_df <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")

# 筛选出要显示标签的前10个上调和下调基因
top_genes_upregulated <- upregulated_genes_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  head(15)
top_genes_downregulated <- downregulated_genes_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  head(15)


# 感兴趣的基因列表
genes <- c(
  
  "S100A4","NEDD9","FN1","CD44","MMP19","LOXL2","LOX","SERPINE1",
  #"VIM","COL1A1","VCAN","THBS1","SPARC","RHOB","RND3",
  "TIMP2","CXCL1","CXCL2","CCL2","CCL3","PDGFRA","VEGFA","S100A10",
  "JUN","FOS","FOSL1","NFKB1","STAT3"
  
)

# 筛选出感兴趣的基因
interested_genes <- scRNAsub.cluster.markers %>%
  filter(gene %in% genes)

# 绘制火山图
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
    limits = c(-2, 3),                      # 设置 X 轴范围
    breaks = seq(-1.5, 3, by = 1),            # 设置刻度
    expand = expansion(mult = c(0.05, 0.05)) # 增加中间区域扩展
  ) +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0), 
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 18))

# 保存图片
ggsave(file.path(output, "Tumor_vs_Normal_volcano_plot.svg"), p, width = 8, height = 7, dpi = 300)
ggsave(file.path(output, "Tumor_vs_Normal_volcano_plot.pdf"), p, width = 8, height = 7, dpi = 300)



library(ggrepel)
library(ggplot2)

# 添加差异计算列
scRNAsub.cluster.markers <- scRNAsub.cluster.markers %>%
  mutate(Difference = pct.1 - pct.2)


# 绘制火山图
volcano_plot <- ggplot(scRNAsub.cluster.markers, aes(x = Difference, y = avg_log2FC, color = Significance)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("blue", "grey", "red")) + 
  # 标注感兴趣的基因
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
    axis.title.x = element_text(size = 20, face = "bold"),  # X轴标题字体大小加粗
    axis.title.y = element_text(size = 20, face = "bold"),  # Y轴标题字体大小加粗
    axis.text.x = element_text(size = 18),  # X轴文本字体大小
    axis.text.y = element_text(size = 18),  # Y轴文本字体大小
    legend.title = element_text(size = 20, face = "bold"),  # 图例标题字体大小加粗
    legend.text = element_text(size = 18),  # 图例文本字体大小
    plot.title = element_text(size = 22, face = "bold"),  # 图标题字体大小加粗
    plot.subtitle = element_text(size = 20),  # 子标题字体大小
    legend.position = "right",  # 图例位置
    legend.key.size = unit(1.5, "lines")  # 图例项的大小
  ) +
  xlim(-0.35, 0.6) +
  ylim(-3, 2) +
  labs(title = "Tumor vs Normal")  # 设置标题和子标题

# 保存为PDF和SVG
ggsave(file.path(output, "volcano_plot.pdf"), volcano_plot, width = 8, height = 6)
ggsave(file.path(output, "volcano_plot.svg"), volcano_plot, width = 8, height = 6)





###############差异基因在不同组中表达##################
library(Seurat)
library(tidyverse)
library(ggsci)

ScRNA <- scRNAsub
# 过滤在表达矩阵中实际存在的基因
genes <- genes[genes %in% rownames(ScRNA)]

# 提取treatment信息
ScRNA$treatment <- as.factor(ScRNA@meta.data$treatment)
treatment_groups <- levels(ScRNA$treatment)

# 获取表达矩阵
expr_matrix <- GetAssayData(ScRNA, slot = "data")[genes, ]


avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment")$RNA
avg_expr_selected <- avg_expr[, treatment_groups]

# 保存平均表达表
avg_expr_df <- avg_expr_selected %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Gene")
write.table(avg_expr_df, file = paste0(output, "/差异基因表达量.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)




# ✅ 绘制 DotPlot
# ✅ 提取参与绘图的基因（再次确认基因是否在对象中）
genes_to_plot <- genes[genes %in% rownames(ScRNA)]

# 计算Tumor组相对于其他组的标准化表达量
# 首先获取所有组的平均表达量
avg_exp <- Seurat::AverageExpression(ScRNA, 
                                     features = genes_to_plot,
                                     group.by = "treatment",
                                     assays = "RNA")$RNA

# 计算Tumor组相对于其他组的标准化表达量
# 这里我们计算Tumor组表达量减去其他组平均表达量
other_groups <- setdiff(colnames(avg_exp), "Tumor")
tumor_normalized <- avg_exp[,"Tumor"] - rowMeans(avg_exp[,other_groups, drop = FALSE])

# 根据标准化表达量排序基因
# 先按是否为正值分组，再按绝对值大小排序
genes_positive <- names(sort(tumor_normalized[tumor_normalized > 0], decreasing = TRUE))
genes_negative <- names(sort(tumor_normalized[tumor_normalized <= 0], decreasing = FALSE))

# 合并基因顺序
genes_ordered <- c(genes_positive, genes_negative)

# 确保所有要绘制的基因都在排序列表中
genes_to_plot <- intersect(genes_ordered, genes_to_plot)



# 绘制基因表达量点图（根据treatment分组）
plot <- DotPlot(ScRNA, features = genes_to_plot, group.by = "treatment") + 
  RotatedAxis() +
  coord_flip() +
  scale_color_gradientn(colors = c('#CCCCCC', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 22), 
        axis.title.x = element_text(size = 22), 
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20))

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.pdf", sep='/'), plot = plot, width = 7, height = 12)
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.svg", sep='/'), plot = plot, width = 7, height = 12)







# 设置颜色（如果样本较多请酌情调整颜色列表）
col_sample <- c('#FF9999',"#A4CDE1","#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
                '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
                "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
                "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 感兴趣的基因列表
cellmarker <- c(
  
  "S100A4","NEDD9","FN1","CD44","MMP19","LOXL2","LOX","SPP1","SERPINE1",
  "VIM","COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL5A2","COL6A3",
  "VCAN","THBS1","SPARC","TIMP2","RHOB","RND3","CXCL8","CXCL1","CXCL2","CXCL3",
  "CCL2","CCL3","CCL4","PDGFRA","PDGFRB","VEGFA","S100A10","S100A11",
  "JUN","FOS","FOSL1","NFKB1","STAT3"
)



# 颜色映射（确保相同处理组颜色一致）
treatment_levels <- unique(scRNAsub@meta.data$treatment)
treatment_colors <- c()

for (lvl in treatment_levels) {
  if (grepl("Ctrl$", lvl)) {
    treatment_colors <- c(treatment_colors, "#A4CDE1")   # 蓝色
  } else if (grepl("RNAlater$", lvl)) {
    treatment_colors <- c(treatment_colors, "#66CCCC")   # 橙色
  } else if (grepl("Si$", lvl)) {
    treatment_colors <- c(treatment_colors, '#FF9999')   # 绿色
  } else {
    treatment_colors <- c(treatment_colors, "#AAAAAA")   # 默认灰色
  }
}


treatment_colors <- sapply(treatment_levels, function(lvl) {
  if (grepl("Tumor$", lvl)) {
    "#1E90FF"
  } else if (grepl("14d-Ctrl$", lvl)) {
    "#009999"
  } else if (grepl("RNAlater$", lvl)) {
    "#33CCCC"
  } else if (grepl("WF-14W$", lvl)) {
    "#FF3366"
  } else {
    "#9E9E9E"
  }
})


names(treatment_colors) <- treatment_levels


#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("0d", "3d","7d", "14d"))
# 创建存储新的 violin plots（按样本分组）
vln_plots_sample <- list()

for (gene in cellmarker) {
  # 绘制每个基因在不同样本中的小提琴图
  vln_plots_sample[[gene]] <- VlnPlot(
    ScRNA,
    features = gene,
    group.by = "treatment",     # <- 按样本分组
    pt.size = 0,
    cols = treatment_colors
  ) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle(paste(gene))
}

# 保存每个样本中基因表达的小提琴图
pdf(paste0(output, "/cellmarker_VlnPlot_bySample.pdf"), width = 20, height = 36)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 4))
dev.off()

svg(paste0(output, "/cellmarker_VlnPlot_bySample.svg"), width = 20, height = 36)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 4))
dev.off()





#################绘制小提琴图################


#View(ScRNA@meta.data)

# 筛选存在于数据集中的marker基因
existing_markers <- cellmarker[cellmarker %in% rownames(ScRNA[["RNA"]]@data)]
existing_markers <- unique(existing_markers)

# 提取表达数据并转换为适合绘图的数据格式
vln.df <- as.data.frame(ScRNA[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

## 计算该基因全局 99th 分位数（在已去 0 的前提下），并去除 >99th 的“彪高”值
#vln.df <- vln.df %>% filter(!is.na(exp), exp > 0)
q05 <- quantile(vln.df$exp, 0.0001, na.rm = TRUE, names = FALSE)
q99 <- quantile(vln.df$exp, 0.99, na.rm = TRUE, names = FALSE)
#vln.df <- vln.df %>% filter(exp > q05 & exp <= q99)

# 继续原有步骤
anno <- ScRNA@meta.data[, c("CB", "treatment")]
vln.df <- inner_join(vln.df, anno, by = "CB")

# 绘制Violin Plot，将X轴和Y轴调换
plot <- vln.df %>%
  ggplot(aes(exp, treatment)) +
  geom_violin(aes(fill = treatment), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # 调整facet_grid，以基因作为列
  scale_fill_manual(values = treatment_colors) +
  scale_x_continuous("") + scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5,size = 25),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    strip.text = element_text(size = 18, face = "bold"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# 保存Violin Plot
ggsave(filename = paste(output, "spacial_ViolinPlot.pdf", sep = '/'), plot = plot, width = 40,height=4,limitsize = FALSE)
ggsave(filename = paste(output, "spacial_ViolinPlot.svg", sep = '/'), plot = plot, width = 40,height=4,limitsize = FALSE)






#####12.展示已知的细胞marker基因的表达情况####
# 定义要绘制的基因列表
genes <- c("EZR","ABCC3","TIMP1","GDF15","CACNA1D","CD55","PTPRD")

# 循环绘制并保存每个基因的特征图
for (gene in genes) {
  # 提取表达数据
  feature_data <- FetchData(scRNAsub, vars = c(gene, "UMAP_1", "UMAP_2", "treatment"))
  
  # 按表达量排序
  feature_data <- feature_data %>% arrange(!!sym(gene))
  
  # 分组绘图
  p1 <- ggplot(feature_data, aes(x = UMAP_1, y = UMAP_2, color = !!sym(gene))) +
    geom_point(size = 0.8, alpha = 0.8) +
    scale_color_gradientn(colors = c("#330066","#660066","#00CC99","#FFFF99")) +
    theme_minimal() +
    facet_wrap(~treatment) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "right",
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    labs(title = paste0(gene), x = "UMAP_1", y = "UMAP_2", color = gene)
  
  ggsave(filename = file.path(output, paste0("FeaturePlot_", gene, "(不同组).pdf")), plot = p1, device = "pdf", width = 6, height = 3)
}




##### 差异基因分析完成后进行GSEA富集分析 #####
#library(org.Mm.eg.db) # 小鼠数据库
library(org.Hs.eg.db) # 小鼠数据库
library(clusterProfiler)
library(enrichplot)
library(DOSE)

scRNAsub.cluster.markers <- readRDS(file.path(output, "ScRNA.sig.markers.rds"))

# 从差异分析中获取差异基因列表
deg <- scRNAsub.cluster.markers[, c('avg_log2FC', 'p_val_adj')]
colnames(deg) <- c('log2FoldChange', 'pvalue')  # 更改列名

# SYMBOL转换为ENTREZID
gene <- bitr(rownames(deg), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 匹配logFC信息
gene$logFC <- deg$log2FoldChange[match(gene$SYMBOL, rownames(deg))]

# 构建genelist
geneList <- gene$logFC
names(geneList) <- gene$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)  # 根据logFC降序排序

# GSEA分析（KEGG通路）
kk_gse <- gseKEGG(geneList = geneList,
                  organism =  "hsa",  
                  ## 'mmu',  # 小鼠
                  nPerm = 1000,       ##置换次数：表示随机打乱基因集1000次进行模拟
                  minGSSize = 10,     ##用于富集分析的基因集至少需要包含10个基因。
                  pvalueCutoff = 0.25,
                  verbose = FALSE)

# 将ENTREZID转换为可读的SYMBOL名称
kk_gse <- DOSE::setReadable(kk_gse, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(as.data.frame(kk_gse), file = file.path(output, "kk_gse_results.csv"))

# 筛选显著富集的通路 |NES| > 1，p值 < 0.05, FDR < 0.25
kk_gse_cut <- kk_gse[kk_gse$pvalue < 0.05 & kk_gse$p.adjust < 0.25 & abs(kk_gse$NES) > 1, ]

# 上调通路
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0, ]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES, decreasing = TRUE), 10), ]

# 下调通路
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0, ]
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES, decreasing = TRUE), 10), ]

# 绘制上调通路的GSEA图
gseap_up <- gseaplot2(kk_gse,
                      up_gsea$ID,
                      title = up_gsea$Description[1], # 使用上调通路第一个的描述作为标题
                      color = c("#FF4500", "#32CD32"),
                      base_size = 30, # 基础字体大小
                      rel_heights = c(1.5, 0.5, 1), # 副图的相对高度
                      subplots = 1:3, # 显示子图
                      ES_geom = "line", # enrichment score线条样式
                      pvalue_table = F)  # 显示p值表 
gseap_up[[1]]<-gseap_up[[1]]+
  scale_color_viridis_d()+
  theme(plot.title = element_text(size = 28, face = "bold"),  # 增大标题字体大小
        legend.title = element_text(size = 20),               # 增大图例标题字体大小
        legend.text = element_text(size = 20))
gseap_up[[2]]<-gseap_up[[2]]+
  scale_color_viridis_d()

ggsave(file.path(output, "GSEA_up.pdf"), gseap_up, width = 15, height = 12)
ggsave(file.path(output, "GSEA_up.svg"), gseap_up, width = 15, height = 12)

# 绘制下调通路的GSEA图
gseap_down <- gseaplot2(kk_gse,
                        down_gsea$ID,
                        title = "DOWN_GSEA", # 标题
                        color = c("#FF4500", "#32CD32"),
                        base_size = 25, # 基础字体大小
                        rel_heights = c(1.5, 0.5, 1), # 副图的相对高度
                        subplots = 1:3, # 显示子图
                        ES_geom = "line", # enrichment score线条样式
                        pvalue_table = F) # 显示p值表
gseap_down[[1]]<-gseap_down[[1]]+
  scale_color_viridis_d()+
  theme(plot.title = element_text(size = 28, face = "bold"),  # 增大标题字体大小
        legend.title = element_text(size = 20),               # 增大图例标题字体大小
        legend.text = element_text(size = 16))
gseap_down[[2]]<-gseap_down[[2]]+
  scale_color_viridis_d()
ggsave(file.path(output, "GSEA_down.pdf"), gseap_down, width = 18, height = 12)
ggsave(file.path(output, "GSEA_down.svg"), gseap_down, width = 18, height = 12)

# ridgeplot 可视化，显示前15个通路
ridgep <- ridgeplot(kk_gse, 
                    showCategory = 15, 
                    fill = "pvalue",  # 根据p值调整填充颜色
                    core_enrichment = TRUE,
                    label_format = 30,  # 设置轴标签的字符长度，过长则换行
                    orderBy = "NES", 
                    decreasing = FALSE)+
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),  
        legend.text = element_text(size = 12)) 

ggsave(file.path(output, "ridgeplot_GSEA.pdf"), ridgep, width = 10, height = 8)
ggsave(file.path(output, "ridgeplot_GSEA.svg"), ridgep, width = 10, height = 8)


#######GO和KEGG富集分析
# 筛选上调和下调基因
gene_up <- scRNAsub.cluster.markers$gene[scRNAsub.cluster.markers$Significance == "Up"]
gene_down <- scRNAsub.cluster.markers$gene[scRNAsub.cluster.markers$Significance == "Down"]

# 将 SYMBOL 转换为 ENTREZID
gene_up_entrez <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, 
                                                             keys = gene_up, 
                                                             columns = 'ENTREZID', 
                                                             keytype = 'SYMBOL')[,2]))
gene_down_entrez <- as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db, 
                                                               keys = gene_down, 
                                                               columns = 'ENTREZID', 
                                                               keytype = 'SYMBOL')[,2]))

# 进行GO富集分析
go_up <- enrichGO(gene = gene_up_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)
go_down <- enrichGO(gene = gene_down_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)

# 将geneID从ENTREZID转为SYMBOL
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

# 绘制GO富集分析图
go_plot_up <- dotplot(go_up) + ggtitle("Upregulated Genes GO Enrichment")
go_plot_down <- dotplot(go_down) + ggtitle("Downregulated Genes GO Enrichment")

# 保存GO分析结果图
ggsave(file.path(output, 'go_enrich_up_dot.pdf'), plot = go_plot_up, width = 6, height = 6)
ggsave(file.path(output, 'go_enrich_up_dot.svg'), plot = go_plot_up, width = 6, height = 6)

ggsave(file.path(output, 'go_enrich_down_dot.pdf'), plot = go_plot_down, width = 6, height = 5)
ggsave(file.path(output, 'go_enrich_down_dot.svg'), plot = go_plot_down, width = 6, height = 5)

# 进行KEGG富集分析
kegg_up <- enrichKEGG(gene = gene_up_entrez, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.1)
kegg_down <- enrichKEGG(gene = gene_down_entrez, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.1)

# 将geneID从ENTREZID转为SYMBOL
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

# 绘制KEGG富集分析图
kegg_plot_up <- dotplot(kegg_up, showCategory = 20)   + ggtitle("Upregulated Genes KEGG Enrichment")
kegg_plot_down <- dotplot(kegg_down, showCategory = 20) + ggtitle("Downregulated Genes KEGG Enrichment")

# 保存KEGG分析结果图
ggsave(file.path(output, 'kegg_enrich_up_dot.pdf'), plot = kegg_plot_up, width = 7, height = 8)
ggsave(file.path(output, 'kegg_enrich_up_dot.svg'), plot = kegg_plot_up, width = 7, height = 8)

ggsave(file.path(output, 'kegg_enrich_down_dot.pdf'), plot = kegg_plot_down, width = 8, height = 8)
ggsave(file.path(output, 'kegg_enrich_down_dot.svg'), plot = kegg_plot_down, width = 8, height = 8)


# 整合并展示GO和KEGG分析图
combined_plot_GO <- go_plot_up + go_plot_down + plot_layout(guides = 'collect')
combined_plot_KEGG <-  kegg_plot_up + kegg_plot_down + plot_layout(guides = 'collect')

# 保存整合后的图像
ggsave(file.path(output, 'combined_GO_dot.pdf'), plot = combined_plot_GO, width = 13, height =10)
ggsave(file.path(output, 'combined_KEGG_dot.pdf'), plot = combined_plot_KEGG, width = 13, height = 12)



# 准备可视化数据
# 提取GO分析结果
go_up_dt <- as.data.frame(go_up)
go_down_dt <- as.data.frame(go_down)

# 提取感兴趣的GO通路
interested_terms <- c(
  "Focal adhesion", "Ubiquitin mediated proteolysis", 
  "Regulation of actin cytoskeleton", "Hedgehog signaling pathway", 
  "Wnt signaling pathway", "Phospholipase D signaling pathway", 
  "Autophagy - animal", "ECM-receptor interaction", "Sphingolipid signaling pathway"
)

# 提取感兴趣的GO通路
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



# 从GO富集结果中筛选感兴趣的通路
go_up_dt <- go_up_dt[go_up_dt$Description %in% interested_terms, ]



# 提取KEGG分析结果
kegg_up_dt <- as.data.frame(kegg_up)
kegg_down_dt <- as.data.frame(kegg_down)

# 提取感兴趣的GO通路
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

# 从GO富集结果中筛选感兴趣的通路
kegg_up_dt <- kegg_up_dt[kegg_up_dt$Description %in% interested_terms, ]

# 设置颜色分类
classification_colors <- c('#437eb8','#FF6666','#FFCC99','#FF9999', '#80c5d8',"#9999FF",
                           "#FFCCCC","#99CCFF","#FF3366","#CCCCFF","#CC0066","#FFFFCC",
                           "#66CCCC","#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                           "#6699CC","#CC99CC","#FF6699","#FF0000","#6666CC","#FF9966",
                           "#669999","#CC99FF","#FFCCFF")

classification_colors<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#094867",'#E59CC4',"#6666CC",
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # 深蓝→绿→浅绿 梯度
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq 浅色
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


# 文本换行功能，限制字符长度
wrap_text <- function(text, width = 40) {
  sapply(text, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

# GO富集分析的条形图
plot_GO_bar <- function(dt, title) {
  dt <- dt[order(-dt$pvalue, decreasing = TRUE),]  # 先按基因数量排序
  #dt <- dt[order(dt$Count, decreasing = TRUE),]  # 先按基因数量排序
  
  dt <- head(dt,20)  # 选取前20个通路
  dt$Description <- factor(wrap_text(dt$Description), levels = wrap_text(dt$Description))
  
  # 左图：富集p值
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
  
  # 右图：基因数量
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
          panel.border = element_rect(color = "black", fill = NA, size = 1))  # 添加XY轴边框
  
  # 组合两个图
  p_combined <- p1 + p2 + plot_layout(widths = c(2, 1.5))
  return(p_combined)
}

# 绘制GO和KEGG富集分析条形图
go_up_plot <- plot_GO_bar(go_up_dt, "Upregulated Genes GO Enrichment")
go_down_plot <- plot_GO_bar(go_down_dt, "Downregulated Genes GO Enrichment")

kegg_up_plot <- plot_GO_bar(kegg_up_dt, "Upregulated Genes KEGG Enrichment")
kegg_down_plot <- plot_GO_bar(kegg_down_dt, "Downregulated Genes KEGG Enrichment")

# 保存图像
ggsave(file.path(output, 'go_enrich_up_bar.pdf'), plot = go_up_plot, width = 8, height = 6)
#ggsave(file.path(output, 'go_enrich_down_bar.pdf'), plot = go_down_plot, width = 8, height = 6)
ggsave(file.path(output, 'kegg_enrich_up_bar.pdf'), plot = kegg_up_plot, width = 10, height = 8)
#ggsave(file.path(output, 'kegg_enrich_down_bar.pdf'), plot = kegg_down_plot, width = 8, height = 5)

# 组合并保存
#combined_GO_plot <- go_up_plot + go_down_plot + plot_layout(guides = 'collect')
#combined_KEGG_plot <- kegg_up_plot + kegg_down_plot + plot_layout(guides = 'collect')

#ggsave(file.path(output, 'combined_GO_bar.pdf'), plot = combined_GO_plot, width = 13, height = 10)
#ggsave(file.path(output, 'combined_KEGG_bar.pdf'), plot = combined_KEGG_plot, width = 13, height = 10)












#### 批量计算多个基因集的 UCell 并可视化 ####
#BiocManager::install("UCell")
library(UCell)
library(Seurat)
library(viridis)
library(ggplot2)
library(stringr)
library(dplyr)

# =============== 基础设置 ===============
if (!exists("outdir")) outdir <- "."  # 若未定义 outdir，则使用当前目录
output_root <- file.path(outdir, "Ucell")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

# 1. 读取 Seurat 对象
sce <- readRDS("celltype.rds")

# 2. 定义多个感兴趣基因集（命名 list）
gene_sets <- list(

  Cancer = c("EPCAM","WFDC2","KRT8","CLDN4","CD24"),
  
  Ovarian = c("WT1","MUC16","PAX8"),
  
  Digestive_tract = c("KRT7", "KRT20","CDX2"),
  
  Metastatic = c(  "SMAD4", "KRAS", "TP53", "CDKN2A", "SNAI1", "SNAI2", "TWIST1", 
                   "ZEB1", "ZEB2", "FOSL1", "ITGA2", "ITGB1", "CDH1", "MMP2", "MMP9", 
                   "RHOA", "RAC1", "CDC42", "CXCR4", "CXCL12","ACKR3", 
                   "TGFB1", "TGFBR1", "TGFBR2", "AXL", "MET", "EGFR", "MYC", "HIF1A", "VEGFA")


)

# =============== 工具函数 ===============
# 生成分箱断点（优先分位数，若重复过多则 pretty 兜底）
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

# 选择降维（优先 UMAP，否则使用 DefaultDimReduc）
.pick_reduction <- function(obj) {
  if ("umap" %in% Reductions(obj)) return("umap")
  DefaultDimReduc(obj)
}

# =============== 主循环：对每个基因集处理 ===============
for (set_name in names(gene_sets)) {
  message("Processing gene set: ", set_name)
  
  # 为该基因集创建输出目录
  out_dir <- file.path(output_root, set_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 基因可用性检查（仅使用落在表达矩阵中的基因）
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
  
  # 3. 计算 UCell 评分（列名将是 '<set_name>_UCell'）
  sce <- AddModuleScore_UCell(
    sce,
    features = setNames(list(genes_in), set_name),  # 命名元素，保证列名前缀
    name = "_UCell"
  )
  score_col <- paste0(set_name, "_UCell")
  if (!score_col %in% colnames(sce@meta.data)) {
    stop(sprintf("[%s] UCell score column not found: %s", set_name, score_col))
  }
  scores <- sce@meta.data[[score_col]]
  
  # 4A. 16 等分位分箱
  brks <- .make_breaks(scores, n_bins_target = 16)
  n_bins <- length(brks) - 1
  bin_labels <- paste0("B", seq_len(n_bins))
  bin_col <- paste0(set_name, "_bin16")
  sce[[bin_col]] <- cut(
    scores, breaks = brks, include.lowest = TRUE, labels = bin_labels,
    right = TRUE, ordered_result = TRUE
  )
  
  # 4B. 20/60/20 分级（高/中/低）
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
  
  # 5. 可视化
  red_use <- .pick_reduction(sce)
  
  ## 5A. 16-bin 离散可视化（DimPlot）
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
  
  ## 5B. High/Medium/Low 可视化
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
  
  ## 5C. 连续分数 FeaturePlot（整体）
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
  
  ## 5D. 连续分数 FeaturePlot（按 treatment 分组）
  # 若没有 'treatment' 列，将跳过分组图
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
  
  
  
  ## 5E. 新增：按 celltype 绘制 UCell 评分小提琴图
  # 每个基因集一个图，x 为 celltype，y 为该基因集 UCell 评分
  message(sprintf("[%s] 绘制按 celltype 分组的小提琴图。", set_name))
  
  # 如果你的 celltype 列名不是 "celltype"，在此修改即可
  celltype_col <- "celltype"
  
  # 从 meta.data 抽出要画图的数据
  df_plot <- data.frame(
    celltype = sce@meta.data[[celltype_col]],
    score    = sce@meta.data[[score_col]],
    stringsAsFactors = FALSE
  )
  
  # 固定 celltype 的顺序（按字母排序或自己指定顺序）
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
  
  
  
  ## 5F. 新增：按 treatment 绘制 UCell 评分箱线图并添加 P 值
  treatment_col <- "treatment"
  
  df_treat <- data.frame(
    treatment = sce@meta.data[[treatment_col]],
    score     = sce@meta.data[[score_col]],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(score))
  
  # 固定 treatment 顺序
  df_treat$treatment <- factor(
    df_treat$treatment,
    levels = sort(unique(df_treat$treatment))
  )
  
  # 两两比较
  comps_treat <- combn(levels(df_treat$treatment), 2, simplify = FALSE)
  
  # y 轴最大值（用于 P 值高度）
  y_max_treat <- max(df_treat$score, na.rm = TRUE)
  
  p_box_treat <- ggplot(df_treat, aes(x = treatment, y = score, fill = treatment)) +
    geom_boxplot(
      outlier.size = 0.3,
      width = 0.8,
      alpha = 0.85
    ) +
    scale_fill_manual(values = c('#437eb8', "#CC0066", "#FF6666")) +
    
    ## ⭐ 显著性检验（Wilcoxon）
    stat_compare_means(
      comparisons = comps_treat,
      method = "wilcox.test",
      label = "p.signif",      # **** / ** / *
      size = 5,
      step.increase = 0.13     # 多组比较自动向上错开
    ) +
    
    ## ⭐ 给 P 值预留空间
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
  
  # 保存
  ggsave(
    file.path(out_dir, paste0(set_name, "_BoxPlot_treatment_pvalue.pdf")),
    plot  = p_box_treat,
    width = 3,
    height = 4
  )
  
}

# 6. 循环完成后保存更新的对象（已包含所有基因集的评分与分级列）
saveRDS(sce, "celltype.rds")

# 如需查看 meta.data
View(sce@meta.data)








# ---- marker- 组织偏好性分析-Startrac---

# ---- 安装/加载 ----
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

# =============== 基础设置 ===============
if (!exists("outdir")) outdir <- "."     # 若未定义 outdir，则默认当前目录
output_root <- file.path(outdir, "STARTRAC_marker")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

# 载入对象与 meta
sco  <- readRDS("celltype.rds")
data <- sco@meta.data

# 需要循环的基因集名称（对应你前一步UCell脚本生成的 <set>_3class 列）
gene_set_names <- c("Cancer","Ovarian","Digestive_tract","Metastatic")

# =============== 工具函数 ===============
# 给定一个 Ro/e 矩阵，返回配色函数与图例刻度
.make_color_meta <- function(mat) {
  rng <- range(mat, na.rm = TRUE)
  # 以 1 作为视觉中点（若不在范围中，则退化为范围均值）
  mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng)
  col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]),
                                  c("#f6f8e6", "#f9a33e", "red"))
  legend_at <- pretty(rng, n = 5)
  legend_labels <- formatC(legend_at, format = "f", digits = 2)
  list(col_fun = col_fun, legend_at = legend_at, legend_labels = legend_labels)
}

# 判定输入是否适合做统计（≥2 个 cluster & ≥2 个 tissue 且非全 NA）
.check_inputs <- function(df, clm_cluster, clm_tissue) {
  ok <- all(c(clm_cluster, "orig.ident", clm_tissue) %in% colnames(df))
  if (!ok) return(FALSE)
  n_cluster <- df[[clm_cluster]] %>% as.character() %>% unique() %>% length()
  n_tissue  <- df[[clm_tissue]]  %>% as.character() %>% unique() %>% length()
  n_cluster >= 2 && n_tissue >= 1
}

# =============== 主循环：对每个基因集计算 Roe 并绘图 ===============
for (set_name in gene_set_names) {
  message("Processing Roe for gene set: ", set_name)
  
  # 本集合的分类列名，如 "Immunosuppression_3class"
  cluster_col <- paste0(set_name, "_3class")
  
  # 存放输出的子目录
  output <- file.path(output_root, set_name)
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  # --------- 输入检查 ---------
  if (!.check_inputs(data, cluster_col, "treatment")) {
    warning(sprintf("[%s] 需要的列缺失或类别数量不足（需要 >=2 个 cluster）。跳过。",
                    set_name))
    next
  }
  # 若某些样本该列是 NA，先过滤
  df_use <- data %>% filter(!is.na(.data[[cluster_col]]), !is.na(treatment))
  
  # 如果过滤后类别只剩 1 个，也跳过
  if (length(unique(df_use[[cluster_col]])) < 2) {
    warning(sprintf("[%s] 过滤 NA 后 cluster 只剩 1 个，无法计算 Roe。跳过。", set_name))
    next
  }
  
  # --------- 计算 Ro/e ---------
  Roe_df <- tryCatch(
    calTissueDist(
      df_use,
      byPatient       = FALSE,
      colname.cluster = cluster_col,   # 关键：按该集合的3-class分组
      colname.patient = "orig.ident",
      colname.tissue  = "treatment",
      method          = "chisq",       # 可选 "chisq", "fisher", "freq"
      min.rowSum      = 0
    ),
    error = function(e) {
      warning(sprintf("[%s] calTissueDist 出错：%s", set_name, e$message))
      return(NULL)
    }
  )
  if (is.null(Roe_df) || nrow(Roe_df) == 0) {
    warning(sprintf("[%s] Roe 结果为空，跳过绘图。", set_name))
    next
  }
  
  # 统一转矩阵，避免数据框在 min/max 上的坑
  mat <- as.matrix(Roe_df)
  
  # 若全零则跳过绘图（或你也可以改成仍然画一个全零图）
  if (all(is.na(mat)) || max(mat, na.rm = TRUE) == 0) {
    warning(sprintf("[%s] Roe 全为 0 或 NA，跳过绘图。", set_name))
    next
  }
  
  # 颜色映射 & 图例
  col_meta <- .make_color_meta(mat)
  col_fun <- col_meta$col_fun
  legend_at <- col_meta$legend_at
  legend_labels <- col_meta$legend_labels
  
  # --------- 图1：带数值热图 ---------
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
  
  # --------- 图2：+++ 符号热图 ---------
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
      # 背景对比度：浅色背景用黑字，深色背景用白字
      rgb_mean <- tryCatch(mean(col2rgb(fill)), error = function(e) 255)
      text_color <- ifelse(rgb_mean > 127, "black", "white")
      grid.text(symbol, x, y, gp = gpar(fontsize = 16, col = text_color))
    }
  )
  print(p2)
  dev.off()
  
  
}







#####################拟时序分析 monocle2 ###############################
library(Seurat)
library(monocle)
library(igraph)
#devtools::install_version("igraph", version = "2.0.8", repos = "http://cran.us.r-project.org")

#BiocManager::install("monocle")
#install.packages("igraph")
dpi=300

col <- c("#CC0066","#A4CDE1",'#FF9999',"#66CCCC","#CCFFCC","#FFFFCC",'#E5D2DD',"#FFCCCC","#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'monocle2', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data <- readRDS(file_path)
colnames(data@meta.data)
#View(data@meta.data)

# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
#data2 <- subset(data1, idents = c("cDC1", "Activated DC"))
#table(data@meta.data$celltype)

#data <- subset(data2, subset = treatment == c("AS-DC","BM-DC"))

# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
#Cells.sub <- subset(data1@meta.data, celltype == c("cDC1", "Activated DC"))
#summary(Cells.sub$celltype)
#data <- subset(data1, cells=row.names(Cells.sub))
#View(data@meta.data)


##提取原始的表达矩阵并稀疏化：UMI count
#expr_matrix<-as(as.matrix(data@assays$RNA@data), 'sparseMatrix')
expr_matrix <- data@assays$RNA@data
## （可选但强烈建议）去掉全零基因，减少体量 & 提高稳健性
nz_genes <- Matrix::rowSums(expr_matrix) > 0
expr_matrix <- expr_matrix[nz_genes, , drop = FALSE]
print(head(expr_matrix[,1:4]))

##提取表型信息，即细胞信息
p_data<-data@meta.data
rownames(p_data)<-colnames(expr_matrix)
head(data@active.ident)

# 仅保留M1和M2细胞
#p_data <- p_data[p_data$celltype %in% c("M1", "M2"), ]
#expr_matrix <- expr_matrix[, rownames(p_data)]

##提取基因信息
f_data<-data.frame(gene_id=rownames(expr_matrix),gene_short_name=rownames(expr_matrix))
rownames(f_data)<-rownames(expr_matrix)
print(head(f_data))

##构建monocle2对象
fd<-new("AnnotatedDataFrame", data = f_data)
pd<-new("AnnotatedDataFrame", data =p_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 1,
                      expressionFamily = negbinomial.size())

#估计文库大小及分散度--归一化处理
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## 排序基因选择及可视化
##计算每个基因表达的细胞数目
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

##选择在5个及以上的细胞中有表达的基因，加速计算
express_genes <- row.names(subset(fData(cds),num_cells_expressed>=5)) #subset：取子集
str(express_genes)
head(fData(cds))

##对剩余基因做差异分析，差异分析基于不同的celltype, cores=2表示调用2个CPU
diff <- differentialGeneTest(cds[express_genes,],fullModelFormulaStr="~Meta_3class",cores=2) 
head(diff)

##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01)
head(deg)

write.csv(deg,file=paste0(output,"/monocle.DEG.csv"),row.names=FALSE)


## 轨迹构建基因（排序基因）可视化
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
pdf(paste0(output,"/ordergenes.pdf"))
plot_ordering_genes(cds)
dev.off()
ggsave(paste0(output,"/ordergenes.png"),plot_ordering_genes(cds),width = dpi*6, height = dpi*6, units = "px",type='cairo')

##数据降维
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)


######## 轨迹可视化
## Pseudotime表示拟时值，State表示细胞状态，celltype表示细胞类型
## 以不同的类型进行轨迹可视化
types <- c("Pseudotime","State","celltype","treatment","Meta_3class")
for (type in types) {
  if (type == "Pseudotime") {
    # 对于连续型数据（Pseudotime），使用渐变色
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else {
    # 对于离散型数据（State 和 celltype），使用自定义颜色
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # 设置离散颜色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))+
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  ggsave(filename = paste(output, paste0("monocle_", type, ".pdf", sep=""), sep="/"), 
         plot = plot_cell_traj, width=7, height=5)
  
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep=""), sep="/"), 
         plot = plot_cell_traj, width = 7, height = 5)
}



saveRDS(cds,"monocle2.rds")


#####拟时序后计算细胞比例#####

col <- c("#CC0066","#A4CDE1",'#FF9999',"#66CCCC","#CCFFCC","#FFFFCC",'#E5D2DD',"#FFCCCC","#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'monocle2', sep='/')
dir.create(output)

cds<-readRDS("monocle2.rds")

# 提取每个细胞的状态信息
state_data <- pData(cds)$State
celltype_data <- pData(cds)$celltype

# 计算各个状态下不同细胞类型的数量
cell_counts_state <- as.data.frame(table(celltype_data, state_data))
colnames(cell_counts_state) <- c("CellType", "State", "Counts")

# 计算每个状态下不同细胞类型的比例
cell_counts_state$Ratio <- ave(cell_counts_state$Counts, cell_counts_state$State, FUN = function(x) x / sum(x))

# 绘制各状态下不同细胞类型的比例柱状图
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

# 保存图像
ggsave(paste0(output, "/state_celltype_proportion.pdf"), plot = p, width = 7)
ggsave(paste0(output, "/state_celltype_proportion.svg"), plot = p, width = 7)





#根节点的确认
cds <- orderCells(cds,root_state=2)
###重新执行可视化，发现拟时的方向发生改变
types <- c("Pseudotime","State","celltype","treatment")

# 自定义颜色
custom_colors <- c(
  "Tumor" = '#FF6666',       # 红色
  "Normal" = '#E5D2DD',     # 绿色
  "BC-cancer" = '#FF6666',  # 蓝色
  "BC-no_cancer" = '#E5D2DD' # 橙色
)


for (type in types) {
  if (type == "Pseudotime") {
    # 对于连续型数据（Pseudotime），使用渐变色
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # 对于离散型数据（State），使用颜色向量 col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # 使用 col 作为颜色
      theme(legend.text = element_text(size = 14),  # 调整图例文本大小
            legend.title = element_text(size = 16),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # 对于离散型数据（celltype 和 treatment），使用自定义颜色 custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # 使用 custom_colors 作为颜色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  # 保存为 PDF 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".pdf", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 7, height = 5)
  
  # 保存为 SVG 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 7, height = 5)
}




##以细胞状态上色（拆分）
# 生成拆分状态的轨迹图
plot_cell_traj_facet <- plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap("~State", nrow = 1) +
  scale_color_manual(values =col)+
  theme(legend.text = element_text(size = 14),  # 调整图例文本大小
        legend.title = element_text(size = 14))  # 调整图例标题大小

ggsave(filename = paste0(output, "/monocle_state_facet.pdf"), plot = plot_cell_traj_facet, width=10, height=6)
ggsave(filename = paste0(output, "/monocle_state_facet.svg"), plot = plot_cell_traj_facet, width=10, height=6)

## 基于表达趋势热图
topgene <- ordergene[1:50]

interested_genes <- c("Cd44","Itgb1","Itga4","Itgae","Ilk","Rock1","Rock2","Cdh1","Ctnnb1",
                      "Junb","Batf3","Irf8","Id2","Il12a","Il12b","Il15","Il15ra")

plot_pseu_heatmap <- plot_pseudotime_heatmap(cds[topgene,],num_clusters = 2,cores = 1,
                                             show_rownames = T,return_heatmap=T,hmcols = colorRampPalette(c("#1f77b4", "#ffffff", "#FF3366"))(100))
pdf(paste0(output,"/monocle_pheatmap1.pdf"), width = 5, height = 6)
print(plot_pseu_heatmap)
dev.off()
ggsave(paste0(output,"/monocle_pheatmap1.svg"),plot_pseu_heatmap, width = 5, height = 6)



## 关键驱动基因的表达变化图
#选择前4个基因
keygenes <- ordergene[1:8] 
#cds_subset <- cds[keygenes,]

print(ordergene)

# 挑选感兴趣的基因（用户需要提供感兴趣的基因名称或索引）
interested_genes <- c(
  #"Ly96","Slco5a1",
  "Cd44","Itgb1","Itga4","Itgae","Ilk","Rock1","Rock2","Cdh1","Ctnnb1",   
  "Junb","Batf3","Irf8","Id2","Il12a","Il12b","Il15","Il15ra"
  
)

keygenes <- intersect(interested_genes, rownames(cds))  # 确保基因在数据集中存在

# 检查选中的基因
if (length(keygenes) == 0) {
  stop("未找到感兴趣的基因，请检查基因名称是否正确。")
} else {
  message("找到以下感兴趣的基因：", paste(keygenes, collapse = ", "))
}

# 子集数据
cds_subset <- cds[keygenes,]

for (type in types) {
  if (type == "Pseudotime") {
    # 对于连续型数据使用渐变色
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 16),  
            legend.title = element_text(size = 16),
            legend.position = "top",
            axis.title.x = element_text(size = 18),  
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),  
            axis.text.y = element_text(size = 16),
            strip.text = element_text(size = 18)) +
      facet_wrap(~ gene_short_name, ncol = 4)
  } else {
    # 对于离散型数据使用自定义颜色
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_manual(values = col) +  # 应用离散颜色
      theme(legend.text = element_text(size = 16),  
            legend.title = element_text(size = 16),
            legend.position = "top",
            axis.title.x = element_text(size = 18),  
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),  
            axis.text.y = element_text(size = 16),
            strip.text = element_text(size = 18)) +
      facet_wrap(~ gene_short_name, ncol = 4)+
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  pdf(file = paste0(output, "/keygene_", type, ".pdf"),width = 10, height = 9)
  print(plot_cell_pseu)
  dev.off()
  
  ggsave(filename = paste0(output, "/keygene_", type, ".svg"), plot = plot_cell_pseu, width = 10, height = 9)
}





## 挖掘对细胞分化有重要影响的开关基因-geneswitches
## load libraries
devtools::install_github("SGDDNB/GeneSwitches")
library(GeneSwitches)
library(SingleCellExperiment)


cardiac_monocle2 <- cds

monocle:::plot_cell_trajectory(cardiac_monocle2, color_by = "State")
plot_monocle_State(cardiac_monocle2)


## Input log-normalized gene expression, Monocle2 pseudo-time and dimensionality reduction
## Path1 containing cells in states 3,2,1
sce_p1 <- convert_monocle2(monocle2_obj = cardiac_monocle2, 
                           states = c(3,2,1), expdata = logexpdata)
## Path2 containing cells in states 3,2,5
sce_p2 <- convert_monocle2(monocle2_obj = cardiac_monocle2, 
                           states = c(3,2,5), expdata = logexpdata)


##二元化基因表达
### binarize gene expression using gene-specific thresholds
sce_p1 <- binarize_exp(sce_p1, ncores = 3)

### 检查二值化的阈值
h <- hist(assays(sce_p1)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
      xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
  abline(v=0.2, col="blue")}

# 选择0.2作为阈值
bn_cutoff <- 0.2
sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = bn_cutoff)


# 对这些潜在的开关基因进行逻辑回归分析和McFadden’s Pseudo R^2拟时间相关性分析
## fit logistic regression and find the switching pseudo-time point for each gene
## with downsampling. This step takes less than 1 mins
sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE)



## 对排序后的开关基因进行可视化，首先会过滤在细胞中表达小于10%，FDR>0.05,和McFadden's Pseudo R^2（<0.03）(拟合不佳)的基因。
## 选择前15个最佳拟合的基因
sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, topnum = 15)
## 选择前20个最佳拟合的蛋白或TF相关的基因
sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, topnum = 20,
                                genelists = gs_genelists, genetype = c("Surface proteins", "TFs"))
##合并去重
sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])

# 沿着伪时间线绘制所选基因。开启的基因绘制在线条上方，而关闭的基因绘制在线条下方。
plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3)

# 使用monocle的降维聚类来可视化基因表达和逻辑回归拟合图。
plot_gene_exp(sce_p1, gene = "VIM", reduction = "monocleRD", downsample = F)


## 沿伪时间线的排序路径
## 过滤r^2小于0.1的基因
sg_pw <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.1)
## 使用超几何分布，确定时间点
switch_pw <- find_switch_pathway(rowData(sce_p1), sig_FDR = 0.05,
                                 pathways = msigdb_h_c2_c5, sg_pw)
## 删除冗余的路径
switch_pw_reduce <- reduce_pathways(switch_pw, pathways = msigdb_h_c2_c5, 
                                    redundant_pw_rate = 0.8)
##
plot_pathway_density(switch_pw_reduce[1:10,], sg_pw, orderbytime = TRUE)

# 根据自己的研究选择特定的term进行展示
sg_vis <- filter_switchgenes(sce_p1, topnum = 50, pathway_name = c("HALLMARK_MYOGENESIS",
                                                                   "GO_CARDIAC_MUSCLE_TISSUE_DEVELOPMENT"))
plot_timeline_ggplot(sg_vis, timedata=sce_p1$Pseudotime, txtsize=3)











