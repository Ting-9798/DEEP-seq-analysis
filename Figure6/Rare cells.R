



# 清空环境
rm(list = ls())

library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
#install.packages("tidydr")

# 设置工作目录
setwd("D:/R/GS/WH/20250903-稀有细胞/data/")
# 定义文件夹路径
data_dir <- "D:/R/GS/WH/20250903-稀有细胞/data/"

outdir <- "D:/R/GS/WH/20250903-稀有细胞/data/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

folders <- c("0606-WF-7/","0620-WF-14/","0829-0825-WF-PNEC-1/","0627-WF-PNEC-3/")

# 红细胞marker基因列表
HB.genes <- c("Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-y","Hbb-bt","Hbb-bs", "Hbb-bh1", "Hbb-bh2")

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
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # 匹配红细胞基因并计算百分比
  HB_m <- match(HB.genes, rownames(seurat_obj@assays$RNA))
  HB.genes.filtered <- rownames(seurat_obj@assays$RNA)[HB_m]
  HB.genes.filtered <- HB.genes.filtered[!is.na(HB.genes.filtered)]
  
  seurat_obj[["percent.HB"]] <- PercentageFeatureSet(seurat_obj, features = HB.genes.filtered)
  
  # 画出红细胞基因比例分布图
  pdf(paste0(OUTPUT,sample_name, "_percent_HB_genes.pdf"), w=5, h=5)
  print(VlnPlot(seurat_obj, features = "percent.HB", pt.size = 0) + 
          ggtitle(paste(sample_name, "Hemoglobin gene percentage")) +
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.5)))
  dev.off()
  
  # 去除红细胞基因比例过高的细胞
  #seurat_obj <- subset(seurat_obj, subset = percent.HB < 1)
  
  
  # 质控过滤
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
  
  # 数据标准化与高变基因寻找
  #seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 根据文件夹名称分配 orig.ident 分组
  if (str_detect(sample_name, "WF-7")) {
    seurat_obj$orig.ident <- "7d"
    seurat_obj$treatment <- "7d"
  } else if (str_detect(sample_name, "WF-14")) {
    seurat_obj$orig.ident <- "14d"
    seurat_obj$treatment <- "14d"
  } else if (str_detect(sample_name, "WF-PNEC-1")) {
    seurat_obj$orig.ident <- "0d"
    seurat_obj$treatment <- "0d"
  } else if (str_detect(sample_name, "WF-PNEC-3")) {
    seurat_obj$orig.ident <- "3d"
    seurat_obj$treatment <- "3d"
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
#summary(combined_seurat@meta.data$treatment)

# 设置为工作目录并保存整合后的 Seurat 对象
saveRDS(combined_seurat, file = "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/combined_seurat.rds")



#######################Seurat分析#####################
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 设置输出目录
setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("0d", "3d","7d", "14d"))

View(ScRNA@meta.data)
summary(ScRNA$`treatment`)

# 根据Cbr2表达量定义exp列
Cbr2_threshold <- 0  # 假设表达阈值为0，可根据实际需求调整
ScRNA@meta.data$exp <- ifelse(ScRNA@assays$RNA@data["tdTomato", ] > Cbr2_threshold, "tdTomato+", "tdTomato-")

# 只保留tdTomato+细胞
ScRNA <- subset(ScRNA, subset = exp == "tdTomato+")
#View(ScRNA@meta.data)
summary(ScRNA$`treatment`)


# 生成小提琴图，显示质控指标
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 9, height = 6)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "treatment", pt.size = 0,cols = col)
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


svg(paste(OUTPUT, "QC-ViolinPlot.svg"), width = 8, height = 6)
# 小提琴图1：nFeature_RNA
p1 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nFeature_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values =col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 小提琴图2：nCount_RNA
p2 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nCount_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 合并图像
library(patchwork)
(p1 | p2) + plot_layout(ncol = 2)
dev.off()

#QC:基因数与线粒体基因以及RNA数量的相关性
pdf(paste(OUTPUT,"cor-plot.pdf"),width = 15,height = 6)
plot1 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = col)
plot2 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols = col)
CombinePlots(plots = list(plot1, plot2),legend = "right")
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
all.genes <- rownames(ScRNA)
ScRNA<-ScaleData(ScRNA,features = all.genes)

#运行PCA
ScRNA<-RunPCA(ScRNA,npcs = 60)

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


save(ScRNA, file = "ScRNA（批次矫正后分群前）.RData")



#### 7.细胞分群与注释 ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")

load("ScRNA（批次矫正后分群前）.RData")
#细胞分群
#ScRNA <- ScRNA %>% 
#  RunUMAP(dims = 1:20) %>% 
#  RunTSNE(dims = 1:20) %>%
#  FindNeighbors(dims = 1:20)


#细胞分群
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:60,reduction = "harmony",
          n.neighbors = 15,     # 默认30，减小后类群间拉开
          min.dist = 0.1,       # 默认0.3，减小后类群更紧凑，类群之间更分离
          spread = 2,
          seed.use = 42) %>%       # 增大 spread 会扩大类群之间距离
  RunTSNE(dims = 1:60, reduction = "harmony",         # 使用更多的主成分
          perplexity = 50,      # 默认是 30，可以适当提高
          theta = 0.3,          # 更精确的计算
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
ScRNA$seurat_clusters <- ScRNA@active.ident##根据聚类树选择你要的resolution
table(Idents(ScRNA))

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_tsne.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "tsne", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# 单独生成umap图
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

saveRDS(ScRNA, "ScRNA（分群后）.rds")


##########查看"tdTomato"和"Epcam"的表达比例###########
#####合并绘制######
genes <- c("tdTomato", "Epcam")
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
  pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge(filter)", ".pdf")
  svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_merge(filter)", ".svg")
  
  # PDF 图
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
  
  # SVG 图
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


# 需要绘图的分组和基因
treatment_groups <- c("0d", "3d","7d", "14d")
genes <- c("tdTomato", "Epcam")

# 遍历每个 treatment 和 gene 绘图
for (treat in treatment_groups) {
  # Subset 对应 treatment 的细胞
  subset_data <- subset(ScRNA, subset = treatment == treat)
  
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
    plot_title <- paste0(gene, " (", round(expression_ratio, 2), "%) - ", treat)
    
    # PDF 输出路径
    pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, "_filter.pdf")
    svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, "_filter.svg")
    
    # PDF 图
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
    
    # SVG 图
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


col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 设置输出目录
setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)
scedata$`treatment` <- factor(scedata$`treatment`, levels = c("0d", "3d","7d", "14d"))

# 定义不同细胞类型的marker基因
cellmarker <- c(
  "Epcam", "Krt18", "Cd24", "Krt19",           # Epithelial Cells /EOC上皮性卵巢癌
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
  #"Ctsd","Ctss",    # Monocyte 单核细胞
  #"Ccl22", "Gm2a", "H2-Ab1",  # Dendritic cells  DC
  #"Col1a1", "Col1a2", "Dcn",  # Fibroblasts / Stromal 基质细胞
  "Bgn", "Mgp",  "Sparc",        #Stromal 基质细胞 "Pdgfra",
  #"Myh11", "Cnn1", "Smtn",               # Smooth Muscle Cells

  "Pcna", "Top2a", "Cdk1",            # Proliferating cells 增殖细胞
  #"Birc5",    # Proliferative endothelial cell
  "Pdgfrb", "Rgs5",   # Pericytes 周细胞
  "Trp63", "Krt5", "Krt14",    # Basal cells
  #"Muc5ac", "Spdef", "Clca1",  # Goblet cells  杯状细胞
  "Sulf1",                 # Tuft cells
  "Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells 纤毛上皮细胞
  "Tppp3","Tubb4b","Tubb1",   #Cilliated纤毛细胞
  "Scgb1a1","Scgb3a2",     #club cell
  "Ascl1","Mash1","Calca","Calcb","Uchl1","Syp","Resp18","Pcsk1","Scg5","Chgb",    #PNEC (肺神经内分泌细胞)
  #"Igfbp5","Fmnl2","Pcsk2","Cacna2d1","Chgb","Sez6l2", ### 神经元内分泌细胞 Neuroendocrine cells
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
  #"Gzma", "Gzmb",  "Ifng", "Ccl5", "Il2", "Tbx21", # Cytotoxic T cells 细胞毒性T细胞
  #"Icos","Grap2", "Csf2", "Gata3" , "Pdcd1" ,               #T helper cell(Th cell)
  #"S100a8", "S100a9",    # Neutrophils
  #"Hba-a1", "Hbb-bs", "Hba-a2", "Hbb-bt",    #Erythrocytes 红细胞
  #"Tpsab1", "Tpsb2",                                   # Mast cells 肥大细胞
  
  
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
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 12, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 12, height = 8)




############################################################################################
########################### 去除 Endothelial Cells 和 Stromal #############################
# 读取已有文件
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata1 <- readRDS(file_path)

# 去除 Endothelial Cells 和 Stromal
meta_filtered <- subset(scedata1@meta.data, !seurat_clusters %in% c("3","6","9","16","18","21","23","25","27","28"))
scedata <- subset(scedata1, cells = rownames(meta_filtered))

# 绘制 UMAP 图（按celltype）
pdf(file.path(output, "ann_umap_filtered.pdf"), width = 5.5, height = 5)
DimPlot(scedata, group.by = "seurat_clusters", reduction = 'umap', pt.size = 0.1, label = TRUE, label.size = 5, repel = TRUE, cols = col, label.box = TRUE) +
  theme_dr(xlength = 0.15, ylength = 0.15, arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        plot.title = element_blank())
dev.off()

# 绘制分组 UMAP 图（按treatment分组）
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



# 设置输出目录
setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/"

# 保存过滤后的数据
saveRDS(scedata, file.path(outdir, "celltype（去杂）.rds"))

##############重新降维聚类####################

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

file_path <- file.path(outdir, "celltype（去杂）.rds")
ScRNA <- readRDS(file_path)
View(ScRNA@meta.data)


#### 6.归一化与PCA降维 ####
#归一化
ScRNA<-ScaleData(ScRNA)

#运行PCA
ScRNA<-RunPCA(ScRNA,npcs = 60)

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


save(ScRNA, file = "ScRNA（批次矫正后分群前）.RData")



#### 7.细胞分群与注释 ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")

load("ScRNA（批次矫正后分群前）.RData")
#细胞分群
#ScRNA <- ScRNA %>% 
#  RunUMAP(dims = 1:20) %>% 
#  RunTSNE(dims = 1:20) %>%
#  FindNeighbors(dims = 1:20)


#细胞分群
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:30,reduction = "harmony") %>%       # 增大 spread 会扩大类群之间距离
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
ScRNA$seurat_clusters <- ScRNA@active.ident##根据聚类树选择你要的resolution
table(Idents(ScRNA))

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_tsne.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "tsne", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "treatment", cols = col)
dev.off()

# 单独生成umap图
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

saveRDS(ScRNA, "ScRNA（分群后）.rds")


##########查看"tdTomato"和"Epcam"的表达比例###########
#####合并绘制######
genes <- c("tdTomato", "Epcam")
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


# 需要绘图的分组和基因
treatment_groups <- c("0d", "3d","7d", "14d")
genes <- c("tdTomato", "Epcam")

# 遍历每个 treatment 和 gene 绘图
for (treat in treatment_groups) {
  # Subset 对应 treatment 的细胞
  subset_data <- subset(ScRNA, subset = treatment == treat)
  
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
    plot_title <- paste0(gene, " (", round(expression_ratio, 2), "%) - ", treat)
    
    # PDF 输出路径
    pdf_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, ".pdf")
    svg_path <- paste0(OUTPUT, gene, "_FeaturePlot_umap_", treat, ".svg")
    
    # PDF 图
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
    
    # SVG 图
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


col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 设置输出目录
setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)

scedata$`treatment` <- factor(scedata$`treatment`, levels = c("0d", "3d","7d", "14d"))

# 定义不同细胞类型的marker基因
cellmarker <- c(
  "Epcam", "Krt18", "Cd24", "Krt19",           # Epithelial Cells /EOC上皮性卵巢癌
  #"Pecam1", "Cdh5","Gpihbp1",  # Endothelial Cells
  # Endothelial
  #"Pecam1","Kdr","Klf2","Car4",
  # Lymphatic Endothelium
  #"Prox1","Lyve1","Pdpn",
  # Fibroblast
  #"Col1a1","Col1a2","Pdgfra","Dcn","Acta2",
  # Smooth muscle / Pericytes
  #"Acta2","Myh11","Tagln","Pdgfrb","Rgs5",
  #"Lsr","Cxcl10", "Clec4g", "Igfbp7", "Adamts1", "Plpp3", "Iigp1", "Kdr", "Nrp1", "Cyp4b1", "Socs3",  # Endothelial Cell
  #"Ctsd","Ctss",    # Monocyte 单核细胞
  #"Ccl22", "Gm2a", "H2-Ab1",  # Dendritic cells  DC
  #"Col1a1", "Col1a2", "Dcn",  # Fibroblasts / Stromal 基质细胞
  #"Bgn", "Mgp",  "Sparc",        #Stromal 基质细胞 "Pdgfra",
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
  #"Gzma", "Gzmb",  "Ifng", "Ccl5", "Il2", "Tbx21", # Cytotoxic T cells 细胞毒性T细胞
  #"Icos","Grap2", "Csf2", "Gata3" , "Pdcd1" ,               #T helper cell(Th cell)
  #"S100a8", "S100a9",    # Neutrophils
  #"Hba-a1", "Hbb-bs", "Hba-a2", "Hbb-bt",    #Erythrocytes 红细胞
  #"Tpsab1", "Tpsb2",                                   # Mast cells 肥大细胞
  
  "Pcna", "Top2a", "Cdk1",            # Proliferating cells 增殖细胞
  #"Birc5",    # Proliferative endothelial cell
  "Trp63", "Krt5", "Krt14",    # Basal cells
  #"Muc5ac","Muc5b", "Spdef", "Clca1","Gsto1","Tff3",  # Goblet cells  杯状细胞
  #"Sulf1",                 # Tuft cells
  "Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells 纤毛上皮细胞
  "Foxj1","Tppp3","Tubb4b",  #Cilliated纤毛细胞
  "Cldn10","Hp","Aldh1a7","Aox3",   # Clara cells	
  "Upk3a","Cyp2f2","N1icd","Scgb1a1","Scgb3a2","Chad",     # club cells   "Upk3a","Cyp2f2","N1icd",
  "Ascl1","Mash1","Calca","Calcb","Uchl1","Syp","Resp18","Pcsk1","Scg5","Chgb",    #PNEC (肺神经内分泌细胞)
  #"Igfbp5","Fmnl2","Pcsk2","Cacna2d1","Chgb","Sez6l2", ### 神经元内分泌细胞 Neuroendocrine cells
  'Ager', 'Hopx', 'Rtkn2', 'Aqp5',"Cav1 ","Spock2",   #AT1
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  #"Lrg1","Lcn2","Retnla","Il33","Car8","Ank3","Cftr",         # Activated AT2
  #"Birc5","Top2a",       # Proliferating AT2s
  "Cldn4","Sfn","Clu","Krt19","Krt8"     # PATs [pre-AT1 过渡状态]  Epithelial transitional states 
  #"Ndrg1","Cldn4","Krt8","Sprr1a","AW112010"   # DATPs (Damage-associated transition progenitors)
  
  ####肺癌######
  #"Pdgfrb", "Rgs5",   # Pericytes 周细胞
  #"Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells 纤毛上皮细胞
  #"Foxj1", "Tppp3","Tubb4b","Tubb1", "Tp73", "Ccdc7", #Cilliated纤毛细胞
  #"Scgb1a1","Scgb3a2","Calcb","Notch2","Hes1", #club cell 
  #"Ascl1","Calca","Syp","Resp18","Pcsk1","Scg5","Chgb","Sez6l2",    #PNEC
  #'Ager', 'Hopx', 'Rtkn2', 'Aqp5', 'Cav1', 'Cav2', 'Wasl', 'Igfbp2', 'Itgb6', 'Sema3a', 'Sema3e', 'Col4a3',  #AT1
  #'Sftpc', 'Sftpb', 'Sftpd', 'Sftpa1', 'Lamp3', 'Abca3', 'Slc34a2', 'Lpcat1', 'Etv5', 'Il33', 'Cpm'  #AT2
  #"Cd68", "Tspo", "Cd163", "Sepp1","Lgals3","Cd11b", "Apoe",   # Macrophage
  #"Mrc1","Ccl3", "Tmem176a", "Tmem176b", "Arg1",  #Alveolar Macrophages
  #"Il2rb","Nkg7","Adamts14","Klra4","Klrf1", "Prf1",        # Natural Killer cell,NK 
  #"B4galt1", "Plcb1", "Brca1", "Mcm6", "Hells", "Cdt1", "Dtl", "Ung", "Rmi2", # NK T cells
  #"Cd8a","Cd8b", # CD8+ T
  #"Trbc2","Icos", # T cells
  #"Gzma", "Gzmb",  "Ifng", "Ccl5", "Il2", "Tbx21", # Cytotoxic T cells 细胞毒性T细胞
  #"Icos","Grap2", "Csf2", "Gata3" , "Pdcd1" ,               #T helper cell(Th cell)
  #"S100a8", "S100a9", "G0s2", # Neutrophils
  
  
  #"Hba-a1", "Hbb-b1", "Gypa", "Epor", "Alas2",   #Erythrocytes 红细胞
  #"Acta2", "Bgn", "Col1a1", "Lgr6", "Mgp", "Pdgfra", "Pdgfrb", "Sparc", "Thy1",       #Stromal 基质细胞
  #"Acta2", "Myh11", "Cnn1", "Smtn", "Tagln",               # Smooth Muscle Cells
  #"Cst3","Ccr7", "H2-Eb1", "Ccl22", "H2-Aa", "Naaa", "Gm2a", "H2-Ab1", "Ppt1", "Cytyp", "S100a8","S100a9","Cd74", # DC
  
  #"Pax8","Wt1","Ca125","Tp53","Cdh1",                     # Ovarian Cancer Cells
  #"Ctss","Cx3cr1", "Cxcr4", "P2ry12","Sall1","Hmox1","Ctsd" , "Iba1", "Cd11b", "Cd45", "P2ry12", "Tmem119","Fcrls", "Tlr4", "Mafb",    #Microglial cell(小胶质细胞)
  
  
  #"Pecam1", "Vegfr2", "Vwf", "Cdh5",  # Endothelial Cells
  #"Lyz","Cd68","Ms4a6a","Cd1e","Il3ra","Lamp3",            # 髓系细胞Myeloid
  #"Tpsab1","Tpsb2",                                        # Mast cell 肥大细胞
  #"Ctss","Ctsd","Itgal","Cd14" ,"Lyz2","Ccr2","Fcgr3",  "Cxcr3","Cx3cr1","Cd14" ,"Ccr2","Cd62l",        # Monocyte 单核细胞 
  #"Ctss","Cd14", "Ly6c1", "Itgam", "Fcgr1", "Ccr2", "Csfr1", "Mpo", "Lyz2","Cd68",        # Monocyte 单核细胞
  #"Ctss","Cd11c", "Cd80", "Cd86", "Cd14", "Flt3", "Irf8", "Cx3cr1", "Ccr2", "Itgam", "MhcII",  #单核细胞-树突状细胞（Monocyte-dendritic cell）
  
  #"Cd68", "Tspo", "Cd163", "Sepp1", "C1qa",    # Macrophage
  #"Chil3", "Ccl9", "S100a4", "Lyz", "Thbs1", "Ms4a4c", "F10", "Ly6c", "Gda", "Lgals3","Cd11b", "C1qb", "Apoe", "Cd14", "Rnase1", "Mrc1", "Marco",   # Macrophage
  #"Cd105", "Cd73", "Cd90", "Stro-1", "Cd44",               # MSC (Mesenchymal Stem Cells)
  
  #"Map2", "Tubb3", "Nefl", "Syp", "Snap25", "Rbfox3", "Dcx","Syn1" ,"Slc17a7","Setd2", # Neurons
  #"Cd2","Ccr5","Il2rb","Nkg7","Gzma","Adamts14","Klra4","Klrd1","Klrb1", "Klrf1","Il18rap","Ncam1","Ncr1","Gnly","Fcgr3a", "Prf1", "Gzmb", "Ifng",  ## Natural Killer cell,NK 
  #"Cd2","Cd3d","Cd3e", "Cd3g","Cd4", "Cd8a", "Cd8", "Cd28", "Ccr7","Ptprc","Trbc2", "Trac", "Icos","Grap2", # T Cells
  #"Cd8a", "Cd8b","Cd3e", "Tbx21","Pdcd1" ,"Sidt1",                                 # CD8 T Cell
  #"Rora","Csf2", "Gata3" , "Pdcd1"                #T helper cell(Th cell)
  #"Il2ra","Il32","Tnfrsf18","Ikzf2", "Areg","Il5"         # Treg
  
  #"Cd19","Cd79a","Cd79b", "Ms4a1","Mzb1","Ms4a1",          # B cell
  #"Bank1", "Bcl11a", "Cd19", "Cd22", "Cd37", "Cd74", "Cd79a", "Cd79b", "Cxcr4", "Ebf1",  # B Cell 1
  #"Creld2", "Crip1", "Derl3", "Dnajc3", "Eaf2", "Edem1", "Edem2", "Fam46c", "Glipr1", "Gm43291", "H13", "Herpud1", "Hsp90b1", "Igha", "Igkc", "Iglc2", "Jchain",  # B Cell Progenitor
  
  
  #"S100a8", "S100a9", "Ly6g", "Mpo", "Csf3r" ,            #Neutrophil 中性粒细胞
  #"Foxj1", "Tubb4b","Tubb1", "Tp73", "Ccdc7", "Dnah5", "Dnah9", "Cfap298", "Ccdc39", "Ccdc40", "Rsph4a", "Spef2", "Tekt1",  #Cilliated纤毛细胞
  #"Alb", "Apoa1", "Fgb", "Gc", "Ahsg", "Kng1", "Mup3", "Car3", "Gsta3", "Hpd", "Ass1", "Mat1a", "Bhmt", "Fabp1", "Aldob", "Wfdc21",  # Hepatocyte
  #"Alcam", "Ambp", "Ankrd1", "Anxa5", "Atp1b1", "Bicc1", "Ces1d", "Cldn3", "Cldn7", "Clu", "Cp", "Cyr61", "Cystm1", "Dbi", "Ddit4l", "Dsg2",  # Cholangiocyte
  #"Col1a1", "Col1a2", "Col3a1", "Colec11", "Cxcl12", "Cygb", "Dcn", "Dpt", "Ecm1", "Efemp1", "Gsn", "Ifitm1", "Igfbp5", "Igfbp6"  # HSC
  
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

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")


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

# 将细胞类型添加到meta数据中
scedata$celltype <- scedata@active.ident
head(scedata@meta.data)

saveRDS(scedata, "celltype.rds")


library(ggsci)
# 绘制细胞类型的umap图
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

#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom"))


# 绘制细胞类型的umap图
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




##### 添加文本标签的边框 ####
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
umap_data$label <- rownames(umap_data)  # 使用行名作为标签（可以替换为其他标签）

# 创建一个含有每个celltype的中心点坐标的数据框
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# 计算每个celltype的中心点坐标
celltype_centers <- umap_df %>%
  group_by(celltype) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# 绘制细胞类型的umap图并保存为PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype, label = label)) +
  geom_point(size = 0.1) +
  ggrepel::geom_label_repel(
    data = celltype_centers, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),  # 使用对应的颜色
    size = 7, fontface = "bold",   # 设置标签为加粗
    box.padding = 0.5,             
    point.padding = 0.5,           
    segment.color = "grey30",
    label.size = 0.5,             # 设置文本框的边框宽度
    label.r = 0.3,                # 设置圆角
    label.border = "black",       # 设置文本框边框为黑色
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +  # 设置细胞类型的颜色
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20)) + # 居中显示标题
  ggtitle(title_label)  # 添加细胞总数作为标题
dev.off()

# 绘制细胞类型的umap图并保存为SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 7, height = 6)
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype, label = label)) +
  geom_point(size = 0.1) +
  ggrepel::geom_label_repel(
    data = celltype_centers, 
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),  # 使用对应的颜色
    size = 7, fontface = "bold",   # 设置标签为加粗
    box.padding = 0.5,             
    point.padding = 0.5,           
    segment.color = "grey30",
    label.size = 0.5,             # 设置文本框的边框宽度
    label.r = 0.3,                # 设置圆角
    label.border = "black",       # 设置文本框边框为黑色
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +  # 设置细胞类型的颜色
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        legend.position = "none",
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 20)) + # 居中显示标题
  ggtitle(title_label)  # 添加细胞总数作为标题
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
    width = 3.5*length(unique(scedata$treatment)), height = 10.5)

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





###########绘制分群注释点图
library(ggh4x)

# 定义不同细胞类型的marker基因
cellmarker <- c(
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  "Cldn10","Hp","Aldh1a7","Aox3",   # Clara cells	
  "Scgb1a1","Scgb3a2","Chad",     # club cells
  "Dnah5", "Dnah9", "Tekt1",  # Ciliated epithelial cells 纤毛上皮细胞
  "Foxj1","Tppp3","Tubb4b",  #Cilliated纤毛细胞
  "Pcna", "Top2a", "Cdk1",            # Proliferating cells 增殖细胞
  "Trp63", "Krt5", "Krt14",    # Basal cells
  'Ager', 'Hopx', 'Rtkn2', "Cav1 ","Spock2",   #AT1
  "Ascl1","Mash1","Calca","Calcb","Uchl1","Syp","Resp18","Pcsk1","Scg5","Chgb"    #PNEC (肺神经内分泌细胞)
  
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
    legend.title = element_text(size = 18),  # 增大图例标题大小
    legend.text = element_text(size = 16),    # 增大图例文本大小
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

# 保存图像
ggsave(plot = p, filename = paste(output, "ann_DotPlot.pdf",sep = '/'), width = 12, height = 6)
ggsave(plot = p, filename = paste(output, "ann_DotPlot.svg",sep = '/'), width = 12, height = 6)


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

pdf(file = paste(output, "ann_Heatmap.pdf",sep = '/'), width = 7, height = 6)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # 自定义高值颜色
dev.off()

svg(file = paste(output, "ann_Heatmap.svg",sep = '/'), width = 7, height = 6)
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

ScRNA_3d <- subset(ScRNA, subset = treatment == "3d")

# 保存结果
saveRDS(ScRNA_3d, file = file.path(outdir, "celltype_3d.rds"))


# 定义不同细胞类型的marker基因
cellmarker <- c(
  
  "Scgb1a1",    # club cells
  'Lamp3',    #AT2
  "Dnah5",
  "Pcna",            # Proliferating cells 增殖细胞
  "Krt5",  # Basal cells
  'Ager'   #AT1
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
# 设置T细胞激活相关基因
cellmarker <- c(
  "Ascl1","Mash1","Calca","Calcb","Syp","Resp18","Pcsk1","Scg5","Chgb"      #PNEC (肺神经内分泌细胞)
  
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
    scale_color_gradientn(colors = c('#E5D2DD',  "#FF3366")) +  # 设置连续型颜色梯度
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
# 设置T细胞激活相关基因
cellmarker <- c(
  "Trp63","Calca","Mki67"             #"tdTomato","Epcam"   
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











###### 阳性细胞 #########
cellmarker <- c("tdTomato","Epcam")
cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]

# 创建数据框存储平均表达量
avg_expr_data <- data.frame(Gene = character(), Avg_Expression = numeric(), stringsAsFactors = FALSE)

# 循环计算平均表达量并绘图
for (gene in cellmarker) {
  
  # 获取基因表达数据
  gene_expr <- FetchData(ScRNA, vars = gene)
  
  # 确保 gene_expr 是向量而非数据框
  gene_expr_vec <- gene_expr[[gene]]
  
  # 计算该基因的平均表达量
  avg_expression <- mean(gene_expr_vec)
  
  # 依据固定阈值 1.5 计算表达细胞数
  expressed_cells <- sum(gene_expr_vec > 0)
  total_cells <- nrow(ScRNA@meta.data)
  expression_ratio <- expressed_cells / total_cells * 100
  
  # 记录平均表达量
  avg_expr_data <- rbind(avg_expr_data, data.frame(Gene = gene, Avg_Expression = avg_expression))
  
  # 按表达量排序
  ordered_cells <- order(gene_expr_vec, decreasing = FALSE)
  cell_names_ordered <- rownames(ScRNA@meta.data)[ordered_cells]
  
  # 设置标题，标注表达比例
  plot_title <- paste0(gene, " Expression (", round(expression_ratio, 2), "%)")
  
  # 添加基因表达量信息到 Seurat 对象
  ScRNA[[paste0(gene, "_expr")]] <- gene_expr_vec
  
  # 绘制 FeaturePlot
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
  
  # 保存 PDF
  ggsave(paste0(output, "/", gene, "_FeaturePlot_umap.pdf"), plot = p, width = 4, height = 4)
  
  # 保存 SVG
  ggsave(paste0(output, "/", gene, "_FeaturePlot_umap.svg"), plot = p, width = 4, height = 4)
}

# 保存平均表达量到 TXT 文件
write.table(avg_expr_data, file = paste0(output, "/Tcell_activation_genes_avg_expression.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




####### 计算细胞比例 ###########
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


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

# 排序
cell_counts_group$Sample <- factor(cell_counts_group$Sample, levels = c("0d", "3d","7d", "14d"))

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
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


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
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 2*length(unique(scedata$orig.ident)), height = 6, dpi = 800)



############分组############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# 计算每个处理组中每种细胞类型的比例
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

# 排序
cell_counts_treatment$Treatment <- factor(cell_counts_treatment$Treatment, levels = c("0d", "3d","7d", "14d"))

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
ggsave(file_path, plot = p1, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)

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
ggsave(file_path, plot = p2, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 2*length(unique(scedata$treatment)), height = 6, dpi = 800)





#################绘制细胞比例折线面积图###############
# 替换为你提供的颜色方案
cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

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
        axis.text = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  geom_vline(aes(xintercept = "3d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "7d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "14d"), linetype = "dashed", size = 1, colour = "white")

ggsave(paste0(output, "/genecount_treatment_area.pdf"), plot = p_counts_area, width = 6, height = 4, dpi = 800)


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
        axis.text = element_text(color = "black", size = 16),
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  geom_vline(aes(xintercept = "3d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "7d"), linetype = "dashed", size = 1, colour = "white") +
  geom_vline(aes(xintercept = "14d"), linetype = "dashed", size = 1, colour = "white")

ggsave(paste0(output, "/geneRatio_treatment_area.pdf"), plot = p_ratio_area, width = 6, height = 4, dpi = 800)







###############分别绘制每个样本########################
library(ggplot2)
library(dplyr)

# 定义需要绘图的样本
samples_to_plot <- c("0d", "3d","7d", "14d")

# 创建输出目录（如果不存在）
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

for (sample_id in samples_to_plot) {
  
  # 筛选当前样本数据
  cell_counts_group_filtered <- subset(cell_counts_group, Sample == sample_id)
  
  ### —— 细胞数量图 —— ###
  # 聚合并生成图例标签（细胞数）
  cell_counts_group_agg <- aggregate(Counts ~ CellType, cell_counts_group_filtered, sum)
  cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", cell_counts_group_agg$Counts, ")")
  legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)
  
  # 绘制数量柱状图
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
  
  # 保存细胞数量图
  ggsave(paste0(output, "/genecount_", sample_id, ".pdf"), plot = p_counts, width = 6, height = 8, dpi = 800)
  ggsave(paste0(output, "/genecount_", sample_id, ".svg"), plot = p_counts, width = 6, height = 8, dpi = 800)
  
  ### —— 细胞比例图 —— ###
  # 聚合并生成图例标签（比例）
  cell_counts_group_agg <- aggregate(Ratio ~ CellType, cell_counts_group_filtered, sum)
  cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", round(cell_counts_group_agg$Ratio * 100, 2), "%)")
  legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)
  
  # 绘制比例柱状图
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
  
  # 保存细胞比例图
  ggsave(paste0(output, "/geneRatio_", sample_id, ".pdf"), plot = p_ratio, width = 6.5, height = 8, dpi = 800)
  ggsave(paste0(output, "/geneRatio_", sample_id, ".svg"), plot = p_ratio, width = 6.5, height = 8, dpi = 800)
}






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
output <- paste(outdir,'差异分析', sep='/')
dir.create(output, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)

# 寻找 Res 和 Sen 组之间的差异基因
logFCfilter <- 0.25        # 定义 log2FC 过滤值
adjPvalFilter <- 0.05   # 定义矫正后 P 值过滤值

# 寻找 Epi_cisplatin_res 和 Epi_other 组之间的差异基因
scRNAsub.cluster.markers <- FindMarkers(object = scRNAsub, 
                                        ident.1 = "3d",
                                        ident.2 = "0d", 
                                        group.by = "treatment", 
                                        logfc.threshold = 0, 
                                        min.pct = 0.25, 
                                        test.use = "wilcox")
scRNAsub.cluster.markers$gene <- rownames(scRNAsub.cluster.markers)

# 添加显著性标注
scRNAsub.cluster.markers <- scRNAsub.cluster.markers %>%
  mutate(Significance = ifelse(p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter, 
                               ifelse(avg_log2FC > 0, "Up", "Down"), "Normal"))
write.table(scRNAsub.cluster.markers, file = file.path(output,"sig.markers_ann_Tol-aLN_vs_PBS.txt"), sep = "\t",row.names = T, quote = FALSE)

saveRDS(scRNAsub.cluster.markers, file = file.path(output, "ScRNA.sig.markers.rds"))

# 分别保存上调基因和下调基因
upregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")
downregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")
write.csv(upregulated_genes, file = file.path(output, "upregulated_genes_Tol-aLN_vs_PBS.csv"), row.names = TRUE, quote = FALSE)
write.csv(downregulated_genes, file = file.path(output, "downregulated_genes_Tol-aLN_vs_PBS.csv"), row.names = TRUE, quote = FALSE)

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
  filter(p_val_adj < 0.01 & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  head(15)
top_genes_downregulated <- downregulated_genes_df %>%
  filter(p_val_adj < 0.01 & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  head(10)

# 感兴趣的基因列表
genes <- c("Gpnmb","Spp1","Ctsd","Nfkb1","Ifngr1","Camk4","Zeb1","Rora","Cd14","Il1b","Il12rb2","Cxcl2","Fth1","Icos","Stat4","Cd63","Thbs1","Tyrobp")

# 筛选出感兴趣的基因
interested_genes <- scRNAsub.cluster.markers %>%
  filter(gene %in% genes)

# 绘制火山图
p <- ggplot(scRNAsub.cluster.markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), size = 2, shape = 18) +
  scale_color_manual(values = c("#0099CC", "#FFCCCC", "#ca0020")) +
  geom_hline(yintercept = -log10(adjPvalFilter), linetype = "dashed") +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed") +
  #geom_text_repel(data = top_genes_upregulated, aes(label = top_genes_upregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  # 添加感兴趣的基因标注
  geom_text_repel(data = interested_genes, aes(label = gene), size = 5, fontface = "bold", box.padding = 0.6, max.overlaps = 50) +
  theme_classic() +
  labs(title = "Tol-aLN vs PBS", 
       x = "log2 Fold Change", y = "-log10 Adjusted P-value", color = "Significance") +
  #annotate("text", x = 0.35, y = 150, label = paste("Up-regulated:", upregulated_genes), 
  #         hjust = 0, size = 6, color = "#FF3300", fontface = "bold") +
  #annotate("text", x = 0.35, y = 135, label = paste("Down-regulated:", downregulated_genes), 
  #         hjust = 0, size = 6, color = "#6699CC", fontface = "bold") +
  #annotate("text", x = 0.35, y = 120, label = paste("DEGs:", total_diff_genes), 
  #         hjust = 0, size = 6, color = "black", fontface = "bold") +
  scale_x_continuous(
    limits = c(-3, 2),                      # 设置 X 轴范围
    breaks = seq(-3, 2, by = 1),            # 设置刻度
    expand = expansion(mult = c(0.05, 0.05)) # 增加中间区域扩展
  ) +
  scale_y_continuous(
    limits = c(0, 200),                      # 设置 X 轴范围
    breaks = seq(0, 200, by = 50),            # 设置刻度
    expand = expansion(mult = c(0.05, 0.05)) # 增加中间区域扩展
  ) +
  theme(plot.title = element_text(size = 22, face = "bold", hjust = 0), 
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 18))

# 保存图片
ggsave(file.path(output, "Tol-aLN_vs_PBS_volcano_plot.svg"), p, width = 8, height = 7, dpi = 300)
ggsave(file.path(output, "Tol-aLN_vs_PBS_volcano_plot.pdf"), p, width = 8, height = 7, dpi = 300)




#####12.展示已知的细胞marker基因的表达情况####
# 定义要绘制的基因列表
genes <- c("Gpnmb","Spp1","Ctsd","Nfkb1","Ifngr1","Camk4","Zeb1","Rora","Cd14","Il1b","Il12rb2","Cxcl2","Fth1","Icos","Stat4","Cd63","Thbs1","Tyrobp")


#"MYO1E","MYC","KLF6","SAT1","JUN",

# 循环绘制并保存每个基因的特征图
for (gene in genes) {
  # 提取表达数据
  feature_data <- FetchData(scRNAsub, vars = c(gene, "UMAP_1", "UMAP_2", "treatment"))
  
  # 按表达量排序
  feature_data <- feature_data %>% arrange(!!sym(gene))
  
  # 分组绘图
  p1 <- ggplot(feature_data, aes(x = UMAP_1, y = UMAP_2, color = !!sym(gene))) +
    geom_point(size = 0.2, alpha = 0.8) +
    scale_color_gradientn(colors = c("#660066","#00CC99","#FFFF99")) +
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



################批量进行差异分析################
library(scRNAtoolVis)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggrepel)
library(org.Mm.eg.db) # 小鼠数据库
#library(org.Hs.eg.db) # 小鼠数据库
library(clusterProfiler)
library(enrichplot)
library(DOSE)


col <- c('#437eb8','#FF6666',"#FFFFCC",'#FFCC99','#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFCCCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")


# 设置输出目录
setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/"

# 创建输出目录
output <- file.path(outdir, "差异分析")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# 读取数据
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
  
  
  # 计算上调和下调基因数目
  upregulated_genes <- sum(markers$Significance == "Up")
  downregulated_genes <- sum(markers$Significance == "Down")
  total_diff_genes <- upregulated_genes + downregulated_genes
  
  # 分别保存上调基因和下调基因的数据框
  upregulated_genes_df <- markers %>%
    filter(Significance == "Up")
  downregulated_genes_df <- markers %>%
    filter(Significance == "Down")
  
  # 筛选出要显示标签的前10个上调和下调基因
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
  
  # 将 SYMBOL 转换为 ENTREZID
  gene_up_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                               keys = gene_up, 
                                                               columns = 'ENTREZID', 
                                                               keytype = 'SYMBOL')[,2]))
  gene_down_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                                 keys = gene_down, 
                                                                 columns = 'ENTREZID', 
                                                                 keytype = 'SYMBOL')[,2]))
  
  # 进行GO富集分析
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
  
  # 删除 KEGG 分析中 Description 中的物种后缀
  kegg_up@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_up@result$Description)
  kegg_down@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_down@result$Description)
  
  # 将geneID从ENTREZID转为SYMBOL
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






###############差异基因表达##################
library(Seurat)
library(tidyverse)
library(ggsci)

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
ScRNA <- readRDS("celltype.rds")

# 设置待分析的T细胞激活相关基因
genes <- c("Gpnmb","Spp1","Ctsd","Nfkb1","Ifngr1","Camk4","Zeb1","Rora","Cd14","Il1b","Il12rb2","Cxcl2","Fth1","Icos","Stat4","Cd63","Thbs1","Tyrobp")


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
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.pdf", sep='/'), plot = plot, width = 5, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_by_treatment.svg", sep='/'), plot = plot, width = 5, height = 5)









## 推断拟时序起点-CytoTRACE ##

#install
#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") 
library(CytoTRACE2)
library(tidyverse)
library(Seurat)


col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272',"#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")


output <- paste(outdir,'推断拟时序起点CytoTRACE', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data1 <- readRDS(file_path)
summary(data1$celltype)


data1@meta.data$CB <- rownames(data1@meta.data)
sample_CB <- data1@meta.data %>% 
  group_by(celltype) %>% 
  sample_frac(0.3)
sce3 <- subset(data1,CB %in% sample_CB$CB) 
sce3 <- data1


#######输入seurat 对象###########
cytotrace2_result_sce <- cytotrace2(sce3, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'mouse',
                                    seed = 1234)
cytotrace2_result_sce


# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = sce3@meta.data$celltype) %>% 
  set_rownames(., colnames(sce3))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)
# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno

# 使用ggsave将每张图保存为PDF
ggsave(file.path(output, "CytoTRACE2_UMAP.pdf"), plot = plots$CytoTRACE2_UMAP, w=4,h=3)
ggsave(file.path(output, "CytoTRACE2_Potency_UMAP.pdf"), plot = plots$CytoTRACE2_Potency_UMAP,,w=4,h=3)
ggsave(file.path(output, "CytoTRACE2_Relative_UMAP.pdf"), plot = plots$CytoTRACE2_Relative_UMAP,,w=4,h=3)
ggsave(file.path(output, "CytoTRACE2_Boxplot_byPheno.pdf"), plot = plots$CytoTRACE2_Boxplot_byPheno,,w=4,h=3)



############# [CytoTRACE评分最高，代表分化程度低，推断为这个细胞数据集中细胞的起点] ###########

# 绘制FeaturePlot并调整风格
p1 <- FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative", pt.size = 0.5) + 
  scale_colour_gradientn(colours = 
                           c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                             "#66C2A5", "#5E4FA2"), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)

# 使用ggsave将图形保存为PDF
ggsave(file.path(output, "CytoTRACE2_FeaturePlot.pdf"), plot = p1,w=4,h=3)




library(ggpubr)
p1 <- ggboxplot(cytotrace2_result_sce@meta.data, x="celltype", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",#轮廓颜色
                fill="celltype",#填充
                palette = "npg",
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "none")+ #图例放右边 
  theme(
    axis.title.x = element_text(size = 14),  # 增大X轴标题字体
    axis.title.y = element_text(size = 14),  # 增大Y轴标题字体
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),  # 增大X轴文本并倾斜30度
    axis.text.y = element_text(size = 12),  # 增大Y轴文本
    legend.title = element_text(size = 12),  # 增大图例标题字体
    legend.text = element_text(size = 10)  # 增大图例文本字体
  )

# 使用ggsave将图形保存为PDF
ggsave(file.path(output, "CytoTRACE2_Boxplot_byPheno.pdf"), plot = p1,w=4,h=3)






library(ggpubr)

# 对celltype根据CytoTRACE2_Score进行排序
cytotrace2_result_sce$celltype <- with(cytotrace2_result_sce@meta.data, 
                                       reorder(celltype, CytoTRACE2_Score, FUN = median))

view(cytotrace2_result_sce)

p1 <- ggboxplot(cytotrace2_result_sce@meta.data, x = "celltype", y = "CytoTRACE2_Score", 
                width = 0.6, 
                color = "black",  # 轮廓颜色
                fill = "celltype",  # 填充
                palette = "npg",
                xlab = F,  # 不显示x轴的标签
                bxp.errorbar = T,  # 显示误差条
                bxp.errorbar.width = 0.5,  # 误差条大小
                size = 1,  # 箱型图边线的粗细
                outlier.shape = NA,  # 不显示outlier
                legend = "none") +  # 图例放右边 
  theme(
    axis.title.x = element_text(size = 14),  # 增大X轴标题字体
    axis.title.y = element_text(size = 14),  # 增大Y轴标题字体
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),  # 增大X轴文本并倾斜30度
    axis.text.y = element_text(size = 12),  # 增大Y轴文本
    legend.title = element_text(size = 12),  # 增大图例标题字体
    legend.text = element_text(size = 10)  # 增大图例文本字体
  )

# 指定组比较
my_comparisons <- list(c("PNEC", "AT1"), c("PNEC", "AT2"), c("PNEC", "Ciliated"),
                       c("PNEC", "Club cells"), c("PNEC", "Basal cells"), c("PNEC", "Proliferating cells"))

# 添加统计比较
p1 <- p1 + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")  # 添加显著性标签

# 使用ggsave将图形保存为PDF
ggsave(file.path(output, "CytoTRACE2_Boxplot_byPheno.pdf"), plot = p1, w = 4, h = 3)











#####################拟时序分析 monocle2 ###############################
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

output <- paste(outdir,'monocle2', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data1 <- readRDS(file_path)
summary(data1$treatment)

# 去除 0d
data <- subset(data1, subset = treatment != "0d")
summary(data$treatment)




##提取原始的表达矩阵并稀疏化：UMI count
expr_matrix<-as(as.matrix(data@assays$RNA@data), 'sparseMatrix')
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
                      lowerDetectionLimit = 0.1,
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
diff <- differentialGeneTest(cds[express_genes,],fullModelFormulaStr="~celltype",cores=2) 
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
types <- c("Pseudotime", "State", "celltype", "treatment")

# 自定义颜色
custom_colors <- c(
  "Tumor" = '#FF6666',       # 红色
  "Normal" = '#E5D2DD',     # 绿色
  "DCs-cancer" = '#FF6666',  # 蓝色
  "DCs-no_cancer" = '#E5D2DD' # 橙色
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
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
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
         plot = plot_cell_traj, width = 10, height = 5)
  
  # 保存为 SVG 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 10, height = 5)
}


#saveRDS(cds,  file = file.path(output, "monocle2.rds"))
saveRDS(cds,"monocle2.rds")


#####拟时序后计算细胞比例#####

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



file_path <- file.path(outdir, "monocle2.rds")
cds <- readRDS(file_path)

##根节点的确认
cds <- orderCells(cds,root_state=1)
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
         plot = plot_cell_traj, width = 9, height = 5)
  
  # 保存为 SVG 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 9, height = 5)
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
plot_pseu_heatmap <- plot_pseudotime_heatmap(cds[topgene,],num_clusters = 5,cores = 1,
                                             show_rownames = T,return_heatmap=T,hmcols = colorRampPalette(c("#1f77b4", "#ffffff", "#FF3366"))(100))
pdf(paste0(output,"/monocle_pheatmap.pdf"))
print(plot_pseu_heatmap)
dev.off()
ggsave(paste0(output,"/monocle_pheatmap.svg"),plot_pseu_heatmap,width = 6, height = 7)



## 关键驱动基因的表达变化图
#选择前4个基因
keygenes <- ordergene[1:8] 
#cds_subset <- cds[keygenes,]

print(ordergene)

# 挑选感兴趣的基因（用户需要提供感兴趣的基因名称或索引）
interested_genes <- c( "PLOD1", "SLC2A5", "TNFRSF14","TNFRSF9", "ERRFI1", "ISG15", "AURKAIP1", "ENO1")

keygenes <- intersect(interested_genes, rownames(cds))  # 确保基因在数据集中存在

# 检查选中的基因
if (length(keygenes) == 0) {
  stop("未找到感兴趣的基因，请检查基因名称是否正确。")
} else {
  message("找到以下感兴趣的基因：", paste(keygenes, collapse = ", "))
}

# 子集数据
cds_subset <- cds[keygenes,]

# 自定义颜色
custom_colors <- c(
  "Tumor" = '#FF6666',       # 红色
  "Normal" = '#E5D2DD',     # 绿色
  "BC-cancer" = '#FF6666',  # 蓝色
  "BC-no_cancer" = '#E5D2DD' # 橙色
)

for (type in types) {
  if (type == "Pseudotime") {
    # 对于连续型数据使用渐变色
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16),
            legend.position = "top",
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            strip.text = element_text(size = 14)) +
      facet_wrap(~ gene_short_name, ncol = 2)
  } else if (type %in%c("State","celltype")) {
    # 对于离散型数据使用自定义颜色
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_manual(values = col) +  # 应用离散颜色
      theme(legend.text = element_text(size = 14),  
            legend.title = element_text(size = 16),
            legend.position = "top",
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            strip.text = element_text(size = 14)) +
      facet_wrap(~ gene_short_name, ncol = 2)+
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # 对于离散型数据使用自定义颜色
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_manual(values = col) +  # 应用离散颜色
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














#####################拟时序分析 monocle2 ###############################
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

# 设置输出目录
setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/"

output <- paste(outdir,'monocle2(notch)', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data1 <- readRDS(file_path)
summary(data1$treatment)

# 去除 0d
data <- subset(data1, subset = treatment != "0d")
summary(data$treatment)

saveRDS(data, "celltype(no-0d).rds")


##提取原始的表达矩阵并稀疏化：UMI count
expr_matrix<-as(as.matrix(data@assays$RNA@data), 'sparseMatrix')
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
                      lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())

#估计文库大小及分散度--归一化处理
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 设置目标基因列表（Notch 通路相关）
target_genes <- c("Rest", "Numbl", "Numb", "Dlk2", "Dlk1", "Dll3", "Rbpj", "Mamld1", "Maml3", "Maml2", "Maml1",
                  "Rfng", "Mfng", "Lfng", "Pofut2", "Pofut1", "Jag2", "Jag1", "Dll4", "Dll1", 
                  "Nrarp", "Heyl", "Hey2", "Hey1", "Hes7", "Hes6", "Hes5", "Hes3", "Hes2", "Hes1",
                  "Notch4", "Notch3", "Notch2", "Notch1")

# 筛选目标基因中实际存在于表达矩阵中的部分
ordering_genes <- intersect(target_genes, rownames(expr_matrix))

cds <- setOrderingFilter(cds, ordering_genes)

## 由于基因数较少，直接用全部基因进行排序基因设定
#cds <- setOrderingFilter(cds, rownames(expr_matrix))

# 可视化排序基因
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
types <- c("Pseudotime", "State", "celltype", "treatment")

# 自定义颜色
custom_colors <- c(
  "Tumor" = '#FF6666',       # 红色
  "Normal" = '#E5D2DD',     # 绿色
  "DCs-cancer" = '#FF6666',  # 蓝色
  "DCs-no_cancer" = '#E5D2DD' # 橙色
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
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
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
         plot = plot_cell_traj, width = 10, height = 5)
  
  # 保存为 SVG 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 10, height = 5)
}


#saveRDS(cds,  file = file.path(output, "monocle2.rds"))
saveRDS(cds,"monocle2.rds")


#####拟时序后计算细胞比例#####

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











############## 返回多个亚群到大群 ###############

col <- c("#CC0066","#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272',"#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")

setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/"

output <- paste(outdir,"合并多亚群", sep='/')
dir.create(output)

# 读取数据
#DCs <- readRDS('D:/R/GS/HZ/20250324-脾/out1/DCs/celltype.rds')  # 包含亚群注释的细分群
PNEC <- readRDS('D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/PNEC/celltype(Notch2).rds')  # 包含亚群注释的细分群
ScRNA <- readRDS("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/celltype.rds")  # 包含所有细胞的大群

# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
#Cells.sub <- subset(ScRNA@meta.data, celltype == c("T cells","CTL","Macrophages","Monocytes","DCs","Neutrophils","B Cells","Plasma Cells","NK T"))
#summary(Cells.sub$celltype)
#all <- subset(ScRNA, cells=row.names(Cells.sub))
#table(all@meta.data$celltype)
#View(all@meta.data)

# 设定每个亚群的标识
#Idents(DCs) <- "celltype"  # DCs细胞的亚群信息
Idents(PNEC) <- "exp"  # T细胞的亚群信息
Idents(ScRNA) <- "celltype"  # 所有细胞的大群信息

#view(DCs$celltype)

# 提取亚群
#selected_celltypes <- c("DCs","PNEC","T")
#all <- subset(ScRNA, idents = selected_celltypes)
all <- subset(ScRNA)
table(all@meta.data$celltype)

# 将多个亚群的标识合并到大群
Idents(all, cells = colnames(PNEC)) <- Idents(PNEC)
#Idents(all, cells = colnames(Oligo)) <- Idents(Oligo)

# 创建一个新的meta信息列来存储合并后的细胞群体信息
all$celltype_merged <- Idents(all)

# 更新 Idents
Idents(all) <- "celltype_merged"

# 查看合并后的分群统计
table(all@meta.data$celltype_merged)

# 保存合并后的数据为rds文件
saveRDS(all, file = "celltype（合并亚群）.rds")

# 绘制细胞类型的UMAP图并保存为PDF
pdf(paste(output, "ann_umap.pdf", sep = '/'), width = 10, height = 6)
DimPlot(object = all, group.by = "celltype_merged", reduction = 'umap', pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        axis.title.x = element_text(size = 16, face = "bold"),  # 增大X轴标题大小
        axis.title.y = element_text(size = 16, face = "bold"),  # 增大Y轴标题大小
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18),
        plot.title = element_blank())
dev.off()

# 绘制细胞类型的UMAP图并保存为SVG
svg(paste(output, "ann_umap.svg", sep = '/'), width = 10, height = 6)
DimPlot(object = all, group.by = "celltype_merged", reduction = 'umap', pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03),
        axis.title.x = element_text(size = 16, face = "bold"),  # 增大X轴标题大小
        axis.title.y = element_text(size = 16, face = "bold"),  # 增大Y轴标题大小
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18),
        plot.title = element_blank())
dev.off()


# 绘制细胞类型的umap图
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

#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom"))


# 绘制细胞类型的umap图
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
#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom")) +
dev.off()


# 绘制根据处理方式分组的UMAP图并保存为PDF
pdf(paste(output, "ann-diff-umap.pdf", sep = '/'), width = 13, height = 5)
DimPlot(all, reduction = "umap", split.by = "treatment", pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col)+
  theme(
    # 增大标签文本大小
    axis.text.x = element_text(size = 16),  # X轴标签大小
    axis.text.y = element_text(size = 16),  # Y轴标签大小
    axis.title.x = element_text(size = 18, face = "bold"),  # 增大X轴标题大小
    axis.title.y = element_text(size = 18, face = "bold"),  # 增大Y轴标题大小
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # 增大标题大小
    legend.title = element_text(size = 18),  # 增大图例标题大小
    legend.text = element_text(size = 18)    # 增大图例文本大小
  )
dev.off()

# 绘制根据处理方式分组的UMAP图并保存为SVG
svg(paste(output, "ann-diff-umap.svg", sep = '/'), width = 13, height = 5)
DimPlot(all, reduction = "umap", split.by = "treatment", pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE,cols = col)+
  theme(
    # 增大标签文本大小
    axis.text.x = element_text(size = 16),  # X轴标签大小
    axis.text.y = element_text(size = 16),  # Y轴标签大小
    axis.title.x = element_text(size = 18, face = "bold"),  # 增大X轴标题大小
    axis.title.y = element_text(size = 18, face = "bold"),  # 增大Y轴标题大小
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # 增大标题大小
    legend.title = element_text(size = 18),  # 增大图例标题大小
    legend.text = element_text(size = 18)    # 增大图例文本大小
  )
dev.off()









#############################不同时间点间的细胞通讯比较（0d, 3d, 7d, 14d）###########################
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

# 定义颜色
col <- c('#437eb8','#FF6666','#FFCC99','#FF9999',"#FFCCCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFFFCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")

col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272',"#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")


# 读取数据
file_path <- file.path(outdir, "celltype.rds")
seuratdata <- readRDS(file_path)
head(seuratdata@meta.data)


# 设置输出目录
output <- paste(outdir,"cellchat", sep='/')
dir.create(output, recursive = TRUE)
setwd(output)


# -----------------------------
#seuratdata <- subset(seuratdata,subset = celltype != "Proliferating cells")

# 清理因子水平（非常重要）
seuratdata@meta.data$celltype <- droplevels(seuratdata@meta.data$celltype)

# 检查是否已删除
table(seuratdata@meta.data$celltype)

# 检查treatment中的时间点
print("Available treatment levels:")
print(unique(seuratdata@meta.data$treatment))

# 提取各时间点的数据
time_points <- c("3d", "7d", "14d")
seurat_subset <- list()

for (tp in time_points) {
  seurat_subset[[tp]] <- subset(seuratdata, subset = treatment %in% tp)
  print(paste("Number of cells in", tp, ":", ncol(seurat_subset[[tp]])))
}

# 创建CellChat对象列表
cellchat_list <- list()

for (tp in time_points) {
  cat("Processing", tp, "...\n")
  
  # 创建CellChat对象
  cellchat_obj <- createCellChat(
    object = seurat_subset[[tp]]@assays$RNA@data, 
    meta = seurat_subset[[tp]]@meta.data, 
    group.by = "celltype"
  )
  
  # 检查和清理未使用的因子水平
  cellchat_obj@idents <- droplevels(cellchat_obj@idents)
  print(paste("Cell types in", tp, ":"))
  print(levels(cellchat_obj@idents))
  
  # 设置数据库
  cellchat_obj@DB <- CellChatDB.mouse
  
  # 分析流程
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- computeCommunProb(cellchat_obj, raw.use = TRUE, population.size = TRUE)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
  
  # 保存单个对象
  saveRDS(cellchat_obj, file = paste0("cellchat_", tp, ".rds"))
  
  cellchat_list[[tp]] <- cellchat_obj
}

# 合并所有时间点的CellChat对象
cellchat <- mergeCellChat(
  cellchat_list, 
  add.names = names(cellchat_list), 
  cell.prefix = TRUE
)
saveRDS(cellchat, "cellchat.rds")

# 如果需要读取已保存的对象
# cellchat_list <- list()
# for (tp in time_points) {
#   cellchat_list[[tp]] <- readRDS(paste0("cellchat_", tp, ".rds"))
# }
# cellchat_merged <- readRDS("cellchat_merged.rds")

# 对比分析 - 相互作用数量和强度
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



################### 配体-受体对比分析##########################
## 展示所有配体受体对的差异
levels(cellchat@idents$joint)
#levels(cellchat_obj@idents)
#levels(cellchat_normal@idents)


##### DCs #####
p <- netVisual_bubble(cellchat, sources.use = c(7), targets.use = c(1,2,3,4,5,6),comparison = c(1, 2,3),angle.x = 45)
# 增大图例文本和标题大小
p <- p + theme(
  legend.title = element_text(size = 18),  # 图例标题大小
  legend.text = element_text(size = 16),   # 图例文本大小
  plot.title = element_text(size = 18),    # 标题大小
  axis.title = element_text(size = 16),    # 坐标轴标题大小
  axis.text = element_text(size = 16)      # 坐标轴刻度文本大小
)
ggsave("Compare_LR_bubble(PNEC).pdf", p, width = 12, height = 60,limitsize = FALSE)




sig.interactions <- subsetCommunication(cellchat)
unique(sig.interactions$interaction_name)

pairLR.use <- as.data.frame(c("TGFB1_TGFBR1_TGFBR2", "TGFB1_ACVR1B_TGFBR2", "BMP2_BMPR1A_BMPR2", "BMP4_BMPR1A_BMPR2", "GDF15_TGFBR2", "AREG_EGFR",
                              "AREG_EGFR_ERBB2", "HBEGF_EGFR", "HBEGF_EGFR_ERBB2", "NRG1_ERBB2_ERBB3", "NRG1_ERBB3", "FGF1_FGFR1", "FGF2_FGFR1",
                              "FGF9_FGFR3", "GAS6_AXL", "GAS6_MERTK", "GAS6_TYRO3", "MIF_CD74_CD44", "IL34_CSF1R", "CSF1_CSF1R", "FN1_ITGA5_ITGB1", 
                              "LAMA5_ITGA6_ITGB1", "THBS1_ITGA3_ITGB1", "COL4A1_ITGA1_ITGB1", "THBS1_CD47", "THBS1_CD36", "CDH1_CDH1", "WNT5B_FZD5",
                              "WNT3A_FZD1_LRP5", "WNT3A_FZD1_LRP6", "SHH_PTCH1", "JAG1_NOTCH1", "SEMA3A_NRP1_PLXNA2", "SEMA4D_PLXNB1", "VEGFA_VEGFR1", "VEGFA_VEGFR2",
                              
                              "TGFB1_TGFBR1_TGFBR2", "BMP2_BMPR1A_ACVR2A", "BMP4_BMPR1A_ACVR2A", "WNT3A_FZD1_LRP5", "VEGFA_VEGFR1", "FGF1_FGFR1", "EGFR_TGFA", 
                              "NRG1_ERBB3", "CX3CL1_CX3CR1", "IL34_CSF1R", "GAS6_AXL", "GAS6_MERTK", "GAS6_TYRO3", "CD80_CTLA4", "CD86_CTLA4", "ICOSL_CTLA4",
                              "CD44_CD44", "CD44_SDC1", "CD44_SDC4", "SDC1_SDC1", "SDC4_SDC4", "CD44_CD36", "CD44_CD47", "SDC1_CD36", "SDC1_CD47", "SDC4_CD36", 
                              "SDC4_CD47", "THBS1_CD36", "THBS3_CD36", "THBS1_CD47", "THBS3_CD47", "THBS1_SDC1", "THBS3_SDC1", "THBS1_SDC4", "THBS3_SDC4", 
                              "THBS1_CD44", "THBS3_CD44", "THBS1_ITGA1_ITGB1", "THBS3_ITGA1_ITGB1", "THBS1_ITGA3_ITGB1", "THBS3_ITGA3_ITGB1", "THBS1_ITGA6_ITGB1", 
                              "THBS3_ITGA6_ITGB1", "THBS1_ITGA7_ITGB1", "THBS3_ITGA7_ITGB1", "THBS1_ITGA9_ITGB1", "THBS3_ITGA9_ITGB1", "THBS1_ITGAV_ITGB1"))



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

# 使用netVisual_bubble绘制特定配体-受体对的气泡图
p <- netVisual_bubble(cellchat, pairLR.use = pairLR.use, sources.use = c(7), targets.use = c(1,2,3,4,5,6),comparison = c(1, 2,3),angle.x = 45)
p <- p + theme(
  legend.title = element_text(size = 18),  # 图例标题大小
  legend.text = element_text(size = 16),   # 图例文本大小
  plot.title = element_text(size = 18),    # 标题大小
  axis.title = element_text(size = 16),    # 坐标轴标题大小
  axis.text = element_text(size = 16)      # 坐标轴刻度文本大小
)

# 保存为PDF文件
ggsave("Compare_LR_bubble_Selected_pairs(PNEC).pdf", p, width = 10, height = 8)



# 网络分析与热图分析
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


# 热图分析
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



## 细胞互作数量不同组对比网络图
# 获取最大权重
weight.max <- getMaxWeight(cellchat_list, attribute = c("org.idents", "count"))
# 循环绘制细胞互作数量不同组对比的网络图并保存为PDF
for (i in 1:length(cellchat_list)) {
  # 定义PDF文件名
  pdf_filename <- paste0("circle_compare_Number_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Number of interactions - ", names(cellchat_list)[i])
  # 设置图形参数
  par(oma = c(0, 0, 0, 0))
  par(cex = 2)
  # 绘制细胞互作网络图
  netVisual_circle(cellchat_list[[i]]@net$count, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 10)
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off() 
}

## 细胞互作强度不同组对比网络图
weight.max.weight <- getMaxWeight(cellchat_list, attribute =c("org.idents", "weight"))
for (i in 1:length(cellchat_list)) {
  # 定义PDF文件名
  pdf_filename <- paste0("circle_compare_Weight_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Weight of interactions - ", names(cellchat_list)[i])
  # 设置图形参数
  par(oma = c(0, 0, 0, 0))
  par(cex = 2)
  # 绘制细胞互作网络图
  netVisual_circle(cellchat_list[[i]]@net$weight, 
                   weight.scale = TRUE, 
                   label.edge = FALSE, 
                   edge.weight.max = weight.max.weight[2], 
                   edge.width.max = 30)
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off() 
}


#########绘制弦图##########
library(circlize)

# 获取最大权重
weight.max <- getMaxWeight(cellchat_list, attribute = c("org.idents", "count"))

# 循环绘制细胞互作数量不同组对比的弦图并保存为PDF
for (i in 1:length(cellchat_list)) {
  # 定义PDF文件名
  pdf_filename <- paste0("chord_compare_Number_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Number of interactions - ", names(cellchat_list)[i])
  
  # 获取互作数量数据
  interaction_data <- cellchat_list[[i]]@net$count
  
  # 设置弦图参数
  chordDiagram(interaction_data, 
               preAllocate = 1, 
               direction.type = c("diffHeight"), 
               link.arr.type = "big.arrow", 
               annotationTrack = c("name", "grid"),  # 只保留一次定义
               link.border = NA,
               transparency = 0.5)
  
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off() 
}



# 获取最大权重
weight.max.weight <- getMaxWeight(cellchat_list, attribute =c("org.idents", "weight"))

# 循环绘制细胞互作强度不同组对比的弦图并保存为PDF
for (i in 1:length(cellchat_list)) {
  # 定义PDF文件名
  pdf_filename <- paste0("chord_compare_Weight_", names(cellchat_list)[i], ".pdf")
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Weight of interactions - ", names(cellchat_list)[i])
  
  # 获取互作强度数据
  interaction_weight_data <- cellchat_list[[i]]@net$weight
  
  # 设置弦图参数
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



## 保守和特异性信号通路的识别与可视化
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strength.pdf", p, width = 10, height = 10)


saveRDS(cellchat, "cellchat_compare.rds")

## 细胞信号模式对比
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


# 特定信号通路对比
df.net <- subsetCommunication(cellchat_obj)
table(df.net$pathway_name)  ###查看具体有哪些信号通路
#levels(df.net$source)


# 设置感兴趣的信号通路列表
pathways_list <- c( "TGFb", "BMP", "EGF", "FGF", "WNT", "GAS", "THBS", "NOTCH", "VEGF", "COLLAGEN", "LAMININ")


# 遍历信号通路列表并生成图
# 遍历信号通路列表并生成图
for (pathways.show in pathways_list) {
  # 检查pathways.show是否在cellchat_list中的pathways里
  available_pathways <- unique(unlist(lapply(cellchat_list, function(x) x@netP$pathways)))
  if (!(pathways.show %in% available_pathways)) {
    message(paste("Skipping pathway:", pathways.show, "- not found in any cellchat_list"))
    next
  }
  
  # 确保weight.max使用正确的参数
  weight.max <- tryCatch({
    getMaxWeight(cellchat_list, slot.name = "netP", attribute = pathways.show)
  }, error = function(e) {
    message(paste("Error in getMaxWeight for pathway:", pathways.show))
    return(NULL)
  })
  
  weight.max <- getMaxWeight(cellchat_list, slot.name = c("netP"), attribute = pathways.show)
  
  # 遍历生成基因表达热图、贡献分析图、热图
  for (i in 1:length(cellchat_list)) {
    SampleOutput <- paste0("compare_", pathways.show, "_", names(cellchat_list)[i])
    
    # 生成并保存单独的基因表达热图（PNG格式）
    png(file = paste0(SampleOutput, "_GeneExpression_heatmap.png"), width = 1200, height = 1400, res = 300)
    p3 <- plotGeneExpression(cellchat_list[[i]], signaling = pathways.show)
    print(p3)
    dev.off()
    
    # 生成并保存单独的贡献分析图（PNG格式）
    png(file = paste0(SampleOutput, "_contribution.png"), width = 1200, height = 800, res = 300)
    p4 <- netAnalysis_contribution(cellchat_list[[i]], signaling = pathways.show)
    print(p4)
    dev.off()
    
    # 生成并保存单独的热图可视化（PNG格式）
    png(file = paste0(SampleOutput, "_netVisual_heatmap.png"), width = 1200, height = 1200, res = 300)
    p5 <- netVisual_heatmap(cellchat_list[[i]], signaling = pathways.show, color.heatmap = "Reds")
    print(p5)
    dev.off()
  }
  
  # 合并基因表达热图
  gene_expr_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_GeneExpression_heatmap.png"))
  gene_expr_imgs <- lapply(gene_expr_filenames, function(x) ggdraw() + draw_image(x))
  combined_gene_expr_plot <- plot_grid(plotlist = gene_expr_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_GeneExpression_heatmap.pdf"), combined_gene_expr_plot, width = 10, height = 8)
  
  # 合并贡献分析图
  contribution_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_contribution.png"))
  contribution_imgs <- lapply(contribution_filenames, function(x) ggdraw() + draw_image(x))
  combined_contribution_plot <- plot_grid(plotlist = contribution_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_contribution.pdf"), combined_contribution_plot, width = 8, height = 8)
  
  # 合并热图可视化
  heatmap_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_netVisual_heatmap.png"))
  heatmap_imgs <- lapply(heatmap_filenames, function(x) ggdraw() + draw_image(x))
  combined_heatmap_plot <- plot_grid(plotlist = heatmap_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_netVisual_heatmap.pdf"), combined_heatmap_plot, width = 8, height = 8)
  
  
  # 循环绘制细胞互作路径网络图并保存为PDF
  for (i in 1:length(cellchat_list)) {
    # 定义PDF文件名
    pdf_filename <- paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_net.pdf")
    pdf(pdf_filename, width = 10, height = 10)  # 设置PDF文件的大小
    
    plot_title <- paste0(pathways.show, " - ", names(cellchat_list)[i])  # 设置图标题
    # 设置图形参数
    par(oma = c(0, 0, 0, 0))
    par(cex = 2)
    
    # 绘制细胞互作网络图
    netVisual_aggregate(cellchat_list[[i]], signaling = pathways.show, layout = "circle", 
                        edge.weight.max = weight.max[1], edge.width.max = 30,
                        vertex.label.cex = 1.2)
    
    title(main = plot_title, line = 1, cex.main = 1.5)  # 设置标题
    dev.off()  # 关闭PDF设备
  }
  
  # 遍历生成chord图
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
  
  # 合并chord图
  chord_filenames <- lapply(1:length(cellchat_list), function(i) paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_chord.png"))
  chord_imgs <- lapply(chord_filenames, function(x) ggdraw() + draw_image(x))
  combined_chord_plot <- plot_grid(plotlist = chord_imgs, ncol = 2)
  ggsave(paste0("Combined_compare_", pathways.show, "_chord.pdf"), combined_chord_plot, width = 8, height = 8)
}






