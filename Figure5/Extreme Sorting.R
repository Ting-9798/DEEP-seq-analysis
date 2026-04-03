
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


setwd("D:/R/GS/WH/20250901-肺癌极限分选/data/")
# 定义文件夹路径
data_dir <- "D:/R/GS/WH/20250901-肺癌极限分选/data/"

folders <- c("0829-0828-WF-1/")
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
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # 质控过滤
  #seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
  #seurat_obj <- subset(seurat_obj, subset = percent.mt < 25)
  
  # 数据标准化与高变基因寻找
  #seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 添加 treatment 列
  seurat_obj$treatment <- sample_name 
  # 根据文件夹名称分配 treatment 分组
  if (str_detect(sample_name, "WF-1")) {
    seurat_obj$treatment <- "WF-1"
  } else if (str_detect(sample_name, "WF-4")) {
    seurat_obj$treatment <- "PBMC-2"
  } else if (str_detect(sample_name, "WF-5")) {
    seurat_obj$treatment <- "PBMC-3"
  }else {
    seurat_obj$treatment <- "Unknown"
  }
  
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
saveRDS(combined_seurat, file = "D:/R/GS/WH/20250901-肺癌极限分选/out/combined_seurat.rds")


#######################Seurat分析#####################
# 设置输出目录
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("D:/R/GS/WH/20250901-肺癌极限分选/out/")
outdir <- "D:/R/GS/WH/20250901-肺癌极限分选/out/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# 拼接完整路径
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)

# 生成小提琴图，显示质控指标
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 5)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4, group.by = "treatment", pt.size = 0,cols = col)
dev.off()

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


svg(paste(OUTPUT, "QC-ViolinPlot-sample.svg"), width = 16, height = 6)
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



save(ScRNA, file = "ScRNA（批次矫正后分群前）.RData")



#### 7.细胞分群与注释 ####

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#6A4C93',
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA（批次矫正后分群前）.RData")
#细胞分群
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:20,spread = 0.5) %>% 
  #RunTSNE(dims = 1:20) %>%
  FindNeighbors(dims = 1:20)

ScRNA<-FindClusters(ScRNA,resolution =seq(from = 0.1, 
                                          to = 1.0, 
                                          by = 0.1))
#pdf(paste(OUTPUT, "clustree.pdf"),width=10,height=9)
#library(clustree)
#clustree(ScRNA)
#dev.off()

#Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident##根据聚类树选择你要的resolution
table(Idents(ScRNA))

#ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("WF-50", "WF-100","WF-200", "WF-300"))

# 确保 "treatment" 因子水平按照 Non-infected 和 Infected 顺序排列
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 6*length(unique(ScRNA$treatment)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=2,label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 6*length(unique(ScRNA$orig.ident)), height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=2, label = TRUE, repel = TRUE, split.by = "orig.ident", cols = col)
dev.off()

# 单独生成umap图
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)
DimPlot(ScRNA, reduction = "umap",pt.size=2, label = TRUE, repel = TRUE, cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

DimPlot(ScRNA, reduction = "umap", pt.size=2, label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
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
DimPlot(ScRNA, repel = TRUE,pt.size=2, 
        reduction = "umap",
        group.by ="treatment")+
  scale_color_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")
dev.off()

saveRDS(ScRNA, "ScRNA（分群后）.rds")



### 计算双胞率-DoubletFinder ####
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder',force = TRUE)
library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)


output <- paste(outdir,"双胞", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
keloid <- readRDS(file_path)


## (2)pK Identification ----------------------------------------------------------
#这是一个测试最佳参数的过程，运行速度慢
sweep.res.list_keloid <- paramSweep(keloid, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_keloid)
sweep.stats_keloid <- summarizeSweep(sweep.res.list_keloid, GT = FALSE)
bcmvn_keloid <- find.pK(sweep.stats_keloid) #可以看到最佳参数的点
## 所以最佳的参数是：
mpK<-as.numeric(as.character(bcmvn_keloid$pK[which.max(bcmvn_keloid$BCmetric)]))




## (3) Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- keloid@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
DoubletRate = ncol(keloid)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
DoubletRate = 0.042104
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 


#nExp_poi <- round(0.008 *nrow(keloid@meta.data)) 
nExp_poi <- round(DoubletRate*length(keloid$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## (4)最后，使用确定好的参数鉴定Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
keloid <- doubletFinder(keloid, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, sct = F)
keloid <- doubletFinder(keloid, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, sct = F)  # reuse.pANN = FALSE,
# 使用nExp = nExp_poi和nExp = nExp_poi.adj,分别进行doublets鉴定，以便后续确定哪些细胞是Doublet-High Confidience


## Plot results ---------------------------------------------------------------------------
View(keloid@meta.data)

keloid@meta.data[,"DF_hi.lo"] <- keloid@meta.data$DF.classifications_0.25_0.3_25
keloid@meta.data$DF_hi.lo[which(keloid@meta.data$DF_hi.lo == "Doublet" & keloid@meta.data$DF.classifications_0.25_0.3_19 == "Singlet")] <- "Doublet-Low Confidience"
keloid@meta.data$DF_hi.lo[which(keloid@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(keloid@meta.data$DF_hi.lo)
# Doublet-High Confidience  Doublet-Low Confidience                  Singlet 
# 198                       24                     5041 

## 结果展示，分类结果在pbmc@meta.data中
pdf(paste(output, "doubletfinder.pdf",sep = '/'), width=4,height=3)
DimPlot(keloid, reduction = "umap", group.by ="DF.classifications_0.25_0.3_19",cols =c("#DC050C","#1965B0"))+
  ggtitle("DoubletFinder") 
dev.off()





########### 计算双胞率-scDblFinder ########
# 安装并加载所需的R包
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#BiocManager::install("scDblFinder")

# or, to get that latest developments:
# BiocManager::install("plger/scDblFinder")
library(scDblFinder)

output <- paste(outdir,"双胞", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
keloid <- readRDS(file_path)

sce <- as.SingleCellExperiment(keloid)


# 单样本（使用随机方法）
sce <- scDblFinder(sce, dbr=0.1)

# 使用基于集群的方法，只需额外提供clusters参数：
sce <- scDblFinder(sce, clusters="seurat_clusters")


# 多个样本（此处，采用这个）
library(BiocParallel)
#sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(3))
table(sce$scDblFinder.class)
## singlet doublet 
##  122418   18187 

# 查看双胞率统计
doublet_stats <- table(sce$scDblFinder.class)
cat("Doublet statistics:\n")
print(doublet_stats)


# 转换回Seurat对象
keloid <- as.Seurat(sce, counts = "counts", data = "logcounts")

# 将双胞信息添加到原始对象（推荐做法）
# 这样可以保留所有的原始降维信息
keloid_original <- readRDS(file_path)  # 重新读取原始数据
keloid_original$scDblFinder.class <- sce$scDblFinder.class
keloid_original$scDblFinder.score <- sce$scDblFinder.score

# 检查可用的降维方法
cat("Available reductions in final object:\n")
print(Reductions(keloid_original))

# 如果UMAP不存在，重新计算
if (!"umap" %in% Reductions(keloid_original)) {
  cat("UMAP not found, computing UMAP...\n")
  
  # 使用已有的PCA进行计算
  if ("pca" %in% Reductions(keloid_original)) {
    keloid_original <- RunUMAP(keloid_original, reduction = "pca", dims = 1:30)
  } else {
    # 如果没有PCA，先进行标准化和PCA
    keloid_original <- NormalizeData(keloid_original)
    keloid_original <- FindVariableFeatures(keloid_original)
    keloid_original <- ScaleData(keloid_original)
    keloid_original <- RunPCA(keloid_original)
    keloid_original <- RunUMAP(keloid_original, reduction = "pca", dims = 1:30)
  }
}

# 设置颜色方案
doublet_colors <- c("singlet" = "#1965B0", "doublet" = "#DC143C")  # 单胞绿色，双胞红色

# 绘制双胞检测结果图
p1 <- DimPlot(keloid_original, reduction = "umap", 
              group.by = "scDblFinder.class", 
              raster = FALSE,
              cols = doublet_colors) +
  ggtitle("ScDblFinder") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # 增大标题
    axis.title = element_text(size = 16),  # 增大坐标轴标题
    axis.text = element_text(size = 14),   # 增大坐标轴标签
    legend.title = element_text(size = 14, face = "bold"),  # 增大图例标题
    legend.text = element_text(size = 14)  # 增大图例文本
  )

# 查看双胞得分的分布 - 使用自定义颜色梯度
p2 <- FeaturePlot(keloid_original, 
                  features = "scDblFinder.score",
                  reduction = "umap",
                  raster = FALSE,
                  cols = c("lightgrey", "#FF6B6B", "#DC143C")) +  # 从浅灰到红色的渐变
  ggtitle("Doublet Score") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # 增大标题
    axis.title = element_text(size = 16),  # 增大坐标轴标题
    axis.text = element_text(size = 14),   # 增大坐标轴标签
    legend.title = element_text(size = 14, face = "bold"),  # 增大图例标题
    legend.text = element_text(size = 14)  # 增大图例文本
  )


# 合并图形
combined_plot <- p1 | p2

# 打印到屏幕
print(combined_plot)

# 保存为PDF
pdf_file <- file.path(output, "scDblFinder.pdf")
pdf(pdf_file, width = 8, height = 3)  # 设置PDF尺寸
print(combined_plot)
dev.off()


# 可选：分别保存两个图（如果需要更大的单个图）
pdf_file_p1 <- file.path(output, "scDblFinder_classification.pdf")
pdf(pdf_file_p1, width = 4.5, height = 3)
print(p1)
dev.off()

pdf_file_p2 <- file.path(output, "scDblFinder_score_distribution.pdf")
pdf(pdf_file_p2, width = 4, height = 3)
print(p2)
dev.off()


# 保存结果
saveRDS(keloid_original, file.path(output, "scRNA_with_scDblFinder.rds"))

# 输出双胞率统计
doublet_rate <- round(doublet_stats["doublet"] / sum(doublet_stats) * 100, 2)
cat(sprintf("\nDoublet rate: %.2f%%\n", doublet_rate))

# 可选：移除双胞细胞
if (FALSE) {  # 设置为TRUE如果想要移除双胞
  keloid_filtered <- subset(keloid_original, 
                            subset = scDblFinder.class == "singlet")
  cat(sprintf("After filtering: %d cells remaining\n", ncol(keloid_filtered)))
  saveRDS(keloid_filtered, file.path(output, "scRNA_filtered_scDblFinder.rds"))
}







##########查看"tdTomato"和"Epcam"的表达比例###########
#####合并绘制######
genes <- c("EPCAM","KRT7","KRT19")
subset_data <- ScRNA

for (gene in genes) {
  # 获取基因表达数据
  gene_expr <- FetchData(subset_data, vars = gene)
  subset_data[[paste0(gene, "_expr")]] <- gene_expr[[gene]]
  
  # 计算非零表达细胞比例
  expressed_cells <- sum(subset_data[[paste0(gene, "_expr")]] > 0)
  total_cells <- nrow(subset_data@meta.data)
  expression_ratio <- expressed_cells / total_cells * 100
  
  # 从 metadata 中提取表达向量
  expr_vec <- subset_data@meta.data[[paste0(gene, "_expr")]]
  
  # 排序细胞（低表达在下层）
  cells_ordered <- subset_data@meta.data[order(expr_vec, decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)
  
  # 设置标题
  plot_title <- paste0(gene, " (", round(expression_ratio, 2), "%) ")
  
  # PDF 输出路径
  pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge", ".pdf")
  svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge", ".svg")
  
  # PDF 图
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
  
  # SVG 图
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



###########细胞手动注释########

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


setwd("D:/R/GS/WH/20250901-肺癌极限分选/out/")
outdir <- "D:/R/GS/WH/20250901-肺癌极限分选/out/"


output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)

# 定义不同细胞类型的marker基因
cellmarker <- c(
  "SERPINB3", "CEACAM5", "ENO2", "KRT19","GRP", 
  "EPCAM","KRT7","KRT19",   # Cancer cells
  "SCGB1A1", "SCGB3A2",     # Club cells
  #"KCNE1", "FOXJ1", "TPPP3", "TUBB4B", "TUBB", "TP73", "CCDC7", # Ciliated cells 纤毛细胞
  "EPCAM", "KRT18",  "KRT19", "MUC1", "KRT8" ,        # Epithelial Cells /EOC上皮性卵巢癌
  #"STEAP4","CEACAM6","SCGB1A1","MUC5B","MUC5AC","SPDEF","FOXJ1",    # Secretory Cells 分泌性上皮细胞
  "DNAH5", "DNAH9", "TEKT1",  # Ciliated epithelial cells 纤毛上皮细胞
  'AGER', 'CAV1', 'CLIC5' ,'HOPX', 'SEMA3E', 'COL4A3',  #AT1
  "LAMP3", "ABCA3", "SLC34A2", "LPCAT1", "SFTPC", "SFTPA1", "SFTPB", "SFTPD"  # AT2 cells
  
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# 使用DotPlot可视化免疫细胞marker基因的表达
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

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 10, height = 4)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width =10, height = 4)


# 绘制点图
plot <- DotPlot(scedata, features = unique(cellmarker)) + RotatedAxis() +
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18))
# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_2.pdf", sep='/'), plot = plot, width = 16, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_2.svg", sep='/'), plot = plot, width = 16, height = 8)




library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)

col <- c('#E5D2DD',"#FF3366","#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


col<- c(
  # UMAP
  "#CCCCCC","#3C8487", "#D0F199", "#79BC98",  "#094867",'#E59CC4',"#6666CC",
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
  "0"="Unknown", 
  "1"="Unknown", 
  "2"="Lung cancer", 
  "3"="Unknown", 
  "4"="Lung cancer", 
  "5"="Unknown")
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


library(ggsci)
# 绘制细胞类型的umap图
pdf(paste(output, "ann_umap.pdf",sep = '/'), width = 6, height = 4)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        axis.title.x = element_text(size = 12, face = "bold"),  # 增大X轴标题大小
        axis.title.y = element_text(size = 12, face = "bold"),  # 增大Y轴标题大小
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

dev.off() 

#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom"))

# 绘制细胞类型的umap图
svg(paste(output, "ann_umap.svg",sep = '/'), width = 6, height = 4)
DimPlot(object=scedata,group.by = "celltype",reduction='umap',pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        axis.title.x = element_text(size = 12, face = "bold"),  # 增大X轴标题大小
        axis.title.y = element_text(size = 12, face = "bold"),  # 增大Y轴标题大小
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16),
        plot.title = element_blank())

#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=5*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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

svg(paste(output, "ann-diff-umap.svg",sep = '/'),width=5*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=2,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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


###### 添加数量和比例标签 ######
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)

# 计算各类细胞数量与比例
cell_counts <- table(scedata$celltype)
cell_prop <- prop.table(cell_counts)

cell_stats <- data.frame(
  celltype = names(cell_counts),
  count = as.numeric(cell_counts),
  prop = round(100 * as.numeric(cell_prop), 1)
)

# 提取 UMAP 坐标
umap_data <- Embeddings(scedata, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(celltype = scedata$celltype)

# 计算每类的中心点坐标，用于放置标签
centers <- umap_data %>%
  group_by(celltype) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# 合并数量与比例信息
centers <- centers %>%
  left_join(cell_stats, by = "celltype") %>%
  mutate(label = paste0("n = ", count, " cells (", prop, "%)"))

# 绘制 PDF
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

# 绘制 SVG
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






# 定义不同细胞类型的marker基因
cellmarker <- c(
  "EPCAM","MUC1","ERBB2","MET"      # Cancer cells
  
)


# 创建存储所有 RidgePlot、VlnPlot 和 FeaturePlot 的列表
ridge_plots <- list()
vln_plots <- list()
feature_plots <- list()

# 遍历免疫细胞标志基因列表
for (gene in cellmarker) {
  # RidgePlot
  ridge_plots[[gene]] <- RidgePlot(scedata, features = gene, ncol = 1, cols = col) +
    theme(legend.position = "none")  
  
  # VlnPlot
  vln_plots[[gene]] <- VlnPlot(scedata, features = gene, ncol = 1, pt.size = 0, cols = col) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") 
  
  # 获取基因表达数据，确保返回为数值向量
  gene_expr <- as.numeric(FetchData(scedata, vars = gene, assay = "RNA")[[gene]])  # 转换为数值向量
  
  # 检查基因表达是否成功提取
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # 将基因表达数据加入元数据
  scedata@meta.data[[paste0(gene, "_expr")]] <- gene_expr  # 添加表达量信息到元数据
  
  # 使用数值向量进行排序，避免直接索引 ScRNA 对象
  cells_ordered <- scedata@meta.data[order(scedata@meta.data[[paste0(gene, "_expr")]], 
                                           decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # 提取排序后的细胞名称
  
  # FeaturePlot 按排序后的顺序绘制，并使用连续型颜色梯度
  feature_plots[[gene]] <- FeaturePlot(scedata, features = gene, reduction = "umap", 
                                       cells = cell_names_ordered,  # 指定细胞顺序
                                       ncol = 1) + 
    scale_color_gradientn(colors = c(  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
                                       "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
                                       "#F5E24B")) +  # 设置连续型颜色梯度
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
pdf(paste0(output, "/marker_FeaturePlot_umap.pdf"), width = 10, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 2))
dev.off()

svg(paste0(output, "/marker_FeaturePlot_umap.svg"), width = 10, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 2))
dev.off()

pdf(paste0(output, "/marker_VlnPlot_umap.pdf"), width = 10, height = 3)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()
svg(paste0(output, "/marker_VlnPlot_umap.svg"), width = 10, height = 3)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()








# =========================
# 加载必要的包
# =========================
if (!require("ggridges")) install.packages("ggridges")
if (!require("cowplot")) install.packages("cowplot")
if (!require("patchwork")) install.packages("patchwork")
if (!require("grid")) install.packages("grid")
library(ggridges)
library(cowplot)
library(patchwork)
library(grid)
library(Seurat)
library(ggplot2)

# =========================
# 定义 marker 基因
# =========================
cellmarker <- c("EPCAM", "MUC1", "ERBB2", "MET")  # Cancer cells

# 存储组合图
combined_plots <- list()

# =========================
# 遍历基因
# =========================
for (gene in cellmarker) {
  
  # 提取基因表达
  gene_expr <- as.numeric(FetchData(
    scedata,
    vars = gene,
    assay = "RNA"
  )[[gene]])
  
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # UMAP 坐标
  umap_data <- as.data.frame(scedata@reductions$umap@cell.embeddings)
  colnames(umap_data) <- c("UMAP_1", "UMAP_2")
  umap_data$expression <- gene_expr
  
  # =========================
  # 主 UMAP FeaturePlot
  # =========================
  p_main <- FeaturePlot(
    scedata,
    features = gene,
    reduction = "umap",
    pt.size = 0.5,
    order = TRUE
  ) +
    scale_color_gradientn(
      colors = c("#ADD8C0", "#094867"),
      name = "Expression"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12)
    ) +
    labs(
      title = gene,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05))
  
  # =========================
  # 上方：UMAP1 单一山峦
  # =========================
  p_top <- ggplot(umap_data, aes(x = UMAP_1)) +
    
    geom_density(
      aes(y = ..density.. * 5),
      fill = "#8FD3C8",
      alpha = 0.6,
      color = NA
    ) +
    
    geom_smooth(
      aes(y = expression),
      method = "loess",
      span = 0.3,
      color = "#094867",
      fill = "#094867",
      alpha = 0.2,
      linewidth = 0.8
    ) +
    
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL)
  
  # =========================
  # 右侧：UMAP2 单一山峦
  # =========================
  p_right <- ggplot(umap_data, aes(x = UMAP_2)) +
    
    geom_density(
      aes(y = ..density.. * 5),
      fill = "#8FD3C8",
      alpha = 0.6,
      color = NA
    ) +
    
    geom_smooth(
      aes(y = expression),
      method = "loess",
      span = 0.3,
      color = "#094867",
      fill = "#094867",
      alpha = 0.2,
      linewidth = 0.8
    ) +
    
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    coord_flip()
  
  # =========================
  # 空白占位
  # =========================
  p_empty <- ggplot() +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  # =========================
  # patchwork 组合
  # =========================
  combined_plots[[gene]] <- p_top + p_empty + p_main + p_right +
    plot_layout(
      ncol = 2,
      nrow = 2,
      widths = c(4, 1),
      heights = c(1, 4)
    )
}

# =========================
# 输出 PDF
# =========================
if (length(cellmarker) > 0) {
  
  # 每个基因单独一页
  pdf(paste0(output, "/marker_FeaturePlot_individual.pdf"),
      width = 5,
      height = 6)
  for (gene in cellmarker) {
    print(combined_plots[[gene]])
  }
  dev.off()
  
  # 所有基因合并一页
  combined_grid <- cowplot::plot_grid(
    plotlist = combined_plots,
    ncol = 2,
    nrow = ceiling(length(cellmarker) / 2)
  )
  
  pdf(
    paste0(output, "/marker_FeaturePlot_combined.pdf"),
    width = 10,
    height = 6 * ceiling(length(cellmarker) / 2)
  )
  print(combined_grid)
  dev.off()
}











# 加载必要的包
if (!require("ggridges")) install.packages("ggridges")
if (!require("cowplot")) install.packages("cowplot")
if (!require("patchwork")) install.packages("patchwork")
if (!require("grid")) install.packages("grid")
library(ggridges)
library(cowplot)
library(patchwork)
library(grid)

# 定义不同细胞类型的marker基因
cellmarker <- c("EPCAM", "MUC1", "ERBB2", "MET")  # Cancer cells

# 创建存储所有组合图的列表
combined_plots <- list()

# 遍历免疫细胞标志基因列表
for (gene in cellmarker) {
  # 获取基因表达数据
  gene_expr <- as.numeric(FetchData(scedata, vars = gene, assay = "RNA")[[gene]])
  
  # 检查基因表达是否成功提取
  if (is.null(gene_expr) || length(gene_expr) == 0) {
    stop(paste("Failed to fetch data for gene:", gene))
  }
  
  # 获取UMAP坐标
  umap_data <- as.data.frame(scedata@reductions$umap@cell.embeddings)
  colnames(umap_data) <- c("UMAP_1", "UMAP_2")
  umap_data$expression <- gene_expr
  
  # 创建UMAP散点图 - 添加边框和坐标轴标题
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
      # 添加面板边框
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      # 坐标轴设置
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
    # 添加坐标轴刻度
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05))
  
  # 创建表达水平分类
  umap_data$expr_level <- "Low"
  if (any(umap_data$expression > 0)) {
    expr_positive <- umap_data$expression[umap_data$expression > 0]
    median_expr <- quantile(expr_positive, 0.5)
    
    umap_data$expr_level[umap_data$expression > 0 & 
                           umap_data$expression <= median_expr] <- "Medium"
    umap_data$expr_level[umap_data$expression > median_expr] <- "High"
  }
  
  # 设置因子顺序，确保图例顺序正确
  umap_data$expr_level <- factor(umap_data$expr_level, 
                                 levels = c("Low", "Medium", "High"))
  
  # 方法1: 使用密度阴影面积图展示UMAP1轴上的表达分布（带图例）
  p_top <- ggplot(umap_data, aes(x = UMAP_1, fill = expr_level)) +
    # 添加不同表达水平细胞的密度阴影
    geom_density(data = subset(umap_data, expr_level == "Low"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.5, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "Medium"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "High"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    # 添加表达量平滑曲线
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
  
  # 使用密度阴影面积图展示UMAP2轴上的表达分布（带图例）
  p_right <- ggplot(umap_data, aes(x = UMAP_2, fill = expr_level)) +
    # 添加不同表达水平细胞的密度阴影
    geom_density(data = subset(umap_data, expr_level == "Low"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.5, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "Medium"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    geom_density(data = subset(umap_data, expr_level == "High"), 
                 aes(y = ..density.. * 5), 
                 alpha = 0.4, color = NA) +
    # 添加表达量平滑曲线
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
  
  # 创建空白图作为占位符
  p_empty <- ggplot() + 
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  # 使用patchwork组合图形
  combined_plots[[gene]] <- p_top + p_empty + p_main + p_right +
    plot_layout(ncol = 2, nrow = 2, 
                widths = c(4, 1), 
                heights = c(1, 4))
}



# 另外提供一种更紧凑的版本，将所有基因组合在一个图中
if (length(cellmarker) > 0) {
  # 创建一个多页PDF，每页一个基因
  pdf(paste0(output, "/marker_FeaturePlot_individual.pdf"), width = 5, height = 6)
  for (gene in cellmarker) {
    print(combined_plots[[gene]])
  }
  dev.off()
  
  # 创建组合图（所有基因在一页）
  combined_grid <- cowplot::plot_grid(plotlist = combined_plots, ncol = 2, nrow = ceiling(length(cellmarker)/2))
  
  pdf(paste0(output, "/marker_FeaturePlot_combined.pdf"), width = 10, height = 6 * ceiling(length(cellmarker)/2))
  print(combined_grid)
  dev.off()

}











####### 计算细胞比例 ###########
col <- c('#E5D2DD',"#FF3366","#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

col<- c(
  # UMAP
  "#3C8487", "#D0F199", "#79BC98",  "#094867",'#E59CC4',"#6666CC",
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

table(scedata$seurat_clusters)

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

p <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Counts') +
  scale_fill_manual(values = col) +
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # 修改X轴标签
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 20),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小
# 添加细胞数文本标签
p <- p + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)


p <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Ratio') +
  scale_fill_manual(values = col) +
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # 修改X轴标签
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 20),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小
# 添加细胞比例文本标签
p <- p + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 5*length(unique(scedata$orig.ident)), height = 8, dpi = 800)



############分组############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# 计算每个处理组中每种细胞类型的比例
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

########## 绘制细胞数堆叠图 ##########
p1 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 24, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 24),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 26), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 24))  # 修改图例文本大小

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)

########## 绘制细胞比例堆叠图 ##########
p2 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 24, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 24),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 26), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 24))  # 修改图例文本大小

file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 5*length(unique(scedata$treatment)), height = 7, dpi = 800)






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
  "SERPINB3", "CEACAM5", "ENO2", "KRT19","GRP", 
  "EPCAM","KRT7","KRT19"    # Cancer cells
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
    scores >= q80 ~ "cancer_high",
    scores <= q20 ~ "cancer_low",
    TRUE ~ "cancer_medium"
  ),
  levels = c("cancer_high", "cancer_medium", "cancer_low")
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
  labs(title = "Cancer cell score (16-quantile bins)", color = "Bin")

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
  labs(title = "Cancer cell score", color = "Class")


ggsave(filename = file.path(output, "Meta_class_DimPlot.pdf"),
       plot = p_hml, width = 5, height = 4)

# 7. 如仍想保留连续分数的 FeaturePlot（原图），可选保留：
p_cont <- FeaturePlot(
  sce,
  features = score_col,
  order = TRUE,
  ncol = 1,
  cols = c("#CCCCCC","#FF3366")
) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(title = "Cancer cell score", color = "Met_score")

ggsave(filename = file.path(output, "Meta_FeaturePlot.pdf"),
       plot = p_cont, width = 4, height = 3)


### 分组 ###
# 7. 如仍想保留连续分数的 FeaturePlot（原图），可选保留：
p_cont <- FeaturePlot(
  sce,
  features = score_col,
  split.by = 'celltype',
  order = TRUE,
  ncol = 2,
  cols = c("lightgrey", "#FF3366")
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
  colname.patient = "celltype",
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


























######## 非负矩阵（NMF）#######
library(NMF)
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library(gtools)
library(ggplot2)
library(cowplot)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
library(pheatmap)

col <- c('#E5D2DD',"#FF3366","#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'非负矩阵NMF', sep='/')
dir.create(output, showWarnings = FALSE, recursive = TRUE)

file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

# 检查数据对象类型并提取表达矩阵
if (class(scedata)[1] == "Seurat") {
  # 处理Seurat对象
  cat("检测到Seurat对象\n")
  # 使用默认的assay，通常是"RNA"
  default_assay <- DefaultAssay(scedata)
  cat("使用assay:", default_assay, "\n")
  
  # 尝试获取counts矩阵，如果没有则使用data矩阵
  if (!is.null(GetAssayData(scedata, slot = "counts"))) {
    V <- as.matrix(GetAssayData(scedata, slot = "counts"))
    cat("使用counts矩阵\n")
  } else if (!is.null(GetAssayData(scedata, slot = "data"))) {
    V <- as.matrix(GetAssayData(scedata, slot = "data"))
    cat("使用data矩阵（log归一化数据）\n")
  } else {
    stop("无法找到合适的表达矩阵")
  }
  
} else if (class(scedata)[1] == "SingleCellExperiment") {
  # 处理SingleCellExperiment对象
  cat("检测到SingleCellExperiment对象\n")
  
  # 获取可用的assay名称
  assay_names <- assayNames(scedata)
  cat("可用的assay:", paste(assay_names, collapse = ", "), "\n")
  
  # 优先使用counts，然后是logcounts
  if ("counts" %in% assay_names) {
    V <- as.matrix(assay(scedata, "counts"))
    cat("使用counts矩阵\n")
  } else if ("logcounts" %in% assay_names) {
    V <- as.matrix(assay(scedata, "logcounts"))
    cat("使用logcounts矩阵\n")
  } else if (length(assay_names) > 0) {
    V <- as.matrix(assay(scedata, 1))
    cat("使用第一个assay:", assay_names[1], "\n")
  } else {
    stop("没有可用的assay")
  }
  
} else {
  # 如果是普通矩阵或数据框
  cat("检测到普通矩阵/数据框对象\n")
  V <- as.matrix(scedata)
}

# 确保数据为非负
if (any(V < 0)) {
  cat("数据包含负值，进行非负化处理...\n")
  V[V < 0] <- 0
}

# 移除在所有细胞中都为0的基因
gene_sums <- rowSums(V)
V <- V[gene_sums > 0, ]
cat("过滤后数据维度:", dim(V), "\n")

# 定义K值列表
K_values <- 2:8

# 存储所有K值的W和H矩阵
W_list <- list()
H_list <- list()

# 执行不同K值的NMF分解并计算W和H矩阵
set.seed(123)
for (K in K_values) {
  cat("Running NMF for K =", K, "\n")
  
  # 执行NMF分解
  nmf_result <- nmf(V, rank = K, nrun = 10, .options = "v")
  
  # 提取W和H矩阵
  W <- basis(nmf_result)  # 基因模式
  H <- coef(nmf_result)   # 细胞模式
  
  # 保存到列表
  W_list[[as.character(K)]] <- W
  H_list[[as.character(K)]] <- H
}

# 合并所有K值的W和H矩阵，计算平均归一化系数
W_avg <- matrix(0, nrow = nrow(V), ncol = length(K_values))
H_avg <- matrix(0, nrow = ncol(V), ncol = length(K_values))

# 计算每个基因和细胞的平均系数
for (i in 1:length(K_values)) {
  W_avg[, i] <- rowMeans(W_list[[i]], na.rm = TRUE)
  H_avg[, i] <- rowMeans(H_list[[i]], na.rm = TRUE)
}

# 对W和H矩阵进行归一化
W_avg_norm <- scale(W_avg, center = FALSE, scale = apply(W_avg, 2, max))
H_avg_norm <- scale(H_avg, center = FALSE, scale = apply(H_avg, 2, max))

# 设置行名和列名
rownames(W_avg_norm) <- rownames(V)
colnames(W_avg_norm) <- paste0("K", K_values)
rownames(H_avg_norm) <- colnames(V)
colnames(H_avg_norm) <- paste0("K", K_values)

# 计算每个基因的得分，并选择前150个基因
gene_scores <- rowSums(W_avg_norm)
top_genes_idx <- order(gene_scores, decreasing = TRUE)[1:min(150, length(gene_scores))]
top_genes <- rownames(V)[top_genes_idx]

cat("找到前", length(top_genes), "个特征基因\n")

# 保存结果
write.csv(data.frame(Gene = top_genes, Score = gene_scores[top_genes_idx]), 
          file.path(output, "top_genes_NMF.csv"), row.names = FALSE)

# 基因特征矩阵（W_avg_norm）的热图
p1 <- pheatmap(W_avg_norm, 
               main = "Gene Feature Matrix (W) - Average Normalized Coefficients",
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               show_rownames = FALSE,
               show_colnames = TRUE,
               color = colorRampPalette(brewer.pal(9, "Blues"))(50))

# 细胞特征矩阵（H_avg_norm）的热图
p2 <- pheatmap(H_avg_norm, 
               main = "Cell Feature Matrix (H) - Average Normalized Coefficients",
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               show_rownames = FALSE,
               show_colnames = TRUE,
               color = colorRampPalette(brewer.pal(9, "Blues"))(50))

# 选择前150个基因并进行热图展示
W_top150 <- W_avg_norm[top_genes_idx, ]
rownames(W_top150) <- top_genes

# 绘制选出的前150个基因的热图
p3 <- pheatmap(W_top150, 
               main = "Top 150 Genes - Average Normalized Coefficients",
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize_row = 6,
               color = colorRampPalette(brewer.pal(9, "Blues"))(50))

# 保存热图
png(file.path(output, "W_matrix_heatmap.png"), width = 1000, height = 800)
print(p1)
dev.off()

png(file.path(output, "H_matrix_heatmap.png"), width = 1000, height = 800)
print(p2)
dev.off()

png(file.path(output, "top150_genes_heatmap.png"), width = 1000, height = 1200)
print(p3)
dev.off()

# 保存NMF结果
saveRDS(list(W_list = W_list, H_list = H_list, 
             W_avg_norm = W_avg_norm, H_avg_norm = H_avg_norm,
             top_genes = top_genes, gene_scores = gene_scores), 
        file.path(output, "NMF_results.rds"))

# 绘制特征基因得分分布
gene_score_df <- data.frame(Gene = rownames(V), Score = gene_scores)
gene_score_df <- gene_score_df[order(gene_score_df$Score, decreasing = TRUE), ]
gene_score_df$Rank <- 1:nrow(gene_score_df)

# 前150个基因标记
gene_score_df$Top150 <- ifelse(gene_score_df$Rank <= 150, "Top150", "Other")

p4 <- ggplot(gene_score_df, aes(x = Rank, y = Score, color = Top150)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("Top150" = "#FF3366", "Other" = "#3399CC")) +
  labs(title = "Gene Scores from NMF Analysis",
       x = "Gene Rank", y = "NMF Score") +
  theme_minimal()

ggsave(file.path(output, "gene_scores_plot.png"), p4, width = 10, height = 6)

cat("NMF分析完成！结果保存在:", output, "\n")









