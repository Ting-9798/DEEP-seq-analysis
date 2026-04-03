
# 清空环境
rm(list = ls())

# 加载所需的包
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(grid)
library(tidydr)
#install.packages("tidydr")


setwd("D:/R/GS/WH/20250403-8C/data1/")
# 定义文件夹路径
data_dir <- "D:/R/GS/WH/20250403-8C/data1/"

folders <- c("0331-0330-WF/","0331-0329-WF-2/")

# 初始化 Seurat 对象列表
seurat.list <- list()

# 循环读取每个子文件夹中的表达矩阵文件并进行 Seurat 分析
for (folder in folders) {
  folder_path <- file.path(data_dir, folder)
  
  # 读取10X数据格式
  sample_data <- Read10X(data.dir = folder_path)
  sample_name <- basename(folder)
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 10, min.features = 1000)
  
  # 添加线粒体基因比例
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # 计算核糖体蛋白比例
  seurat_obj[["percent.rps"]] <- PercentageFeatureSet(seurat_obj, pattern = c("^RPS","^RPL"))
  
  #seurat_obj <- subset(seurat_obj, subset = percent.rps < 10)
  
  # 质控过滤
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 25)
  
  
  # 数据标准化与高变基因寻找
  #seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  
  # 添加 treatment 列
  seurat_obj$treatment <- sample_name 
  
  
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
saveRDS(combined_seurat, file = "D:/R/GS/WH/20250403-8C/out(WF合并)/combined_seurat.rds")


#######################Seurat分析#####################
# 设置输出目录
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


setwd("D:/R/GS/WH/20250403-8C/out(WF合并)/")
outdir <- "D:/R/GS/WH/20250403-8C/out(WF合并)/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# 拼接完整路径
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)
View(ScRNA@meta.data)

# 生成小提琴图，显示质控指标
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 6)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"), ncol = 4, group.by = "treatment", pt.size = 0,cols = col)
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


svg(paste(OUTPUT, "QC-ViolinPlot.svg"), width = 12, height = 6)
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

p3 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = percent.mt, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "percent.mt", x = "", y = "") +
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
(p1 | p2 | p3 ) + plot_layout(ncol = 3)
dev.off()


#QC:基因数与线粒体基因以及RNA数量的相关性
pdf(paste(OUTPUT,"cor-plot.pdf"),width = 15,height = 6)
plot1 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = col)
plot2 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols = col)
CombinePlots(plots = list(plot1, plot2),legend = "right")
dev.off()

# 计算细胞周期得分
#pdf(paste(OUTPUT, "cellcycle.pdf"), width = 9, height = 6)
#s.genes <- cc.genes$s.genes    ##S期
#g2m.genes <- cc.genes$g2m.genes    ##G2/M期
#ScRNA <- CellCycleScoring(ScRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#VlnPlot(ScRNA, features = c("S.Score", "G2M.Score"), group.by = "treatment", pt.size = 1,cols = col)
#dev.off()

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

#pdf(paste(OUTPUT,  "Dimplot-corret.pdf"),width = 12,height = 6)
#DimPlot(object = ScRNA, reduction = "harmony",
#        pt.size = 0.1, group.by = "treatment")
#dev.off()

#pdf(paste(OUTPUT, "vlnplot-corret.pdf"),width = 12,height = 6)
#VlnPlot(object = ScRNA, features = "harmony_1", 
#        group.by = "treatment", pt.size =0)
#dev.off()


save(ScRNA, file = "ScRNA（批次矫正后分群前）.RData")



#### 7.细胞分群与注释 ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#FEDD81", "#FF9A84",  "#094867",'#E59CC4',"#6666CC",
 "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # 深蓝→绿→浅绿 梯度
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq 浅色
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)


load("ScRNA（批次矫正后分群前）.RData")
#细胞分群
#ScRNA <- ScRNA %>% 
#  RunUMAP(dims = 1:20) %>% 
#  RunTSNE(dims = 1:20) %>%
#  FindNeighbors(dims = 1:20)

ScRNA <- ScRNA %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
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
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident##根据聚类树选择你要的resolution
table(Idents(ScRNA))

# 确保 "treatment" 因子水平按照 Non-infected 和 Infected 顺序排列
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, label.size = 5,repel = TRUE, split.by = "treatment", cols = col)+
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

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 40, height = 5)
DimPlot(ScRNA, reduction = "umap", label = TRUE, label.size = 5,repel = TRUE, split.by = "orig.ident", cols = col)+
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

# 单独生成umap图
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

#将肿瘤与正常展示在一起
pdf(paste(OUTPUT, "cluster-diff_umap.pdf"),width=6,height=6)
DimPlot(ScRNA, repel = TRUE,
        reduction = "umap",
        group.by ="treatment")+
  scale_color_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = c(.01, .1))+
  labs(title = "Sample Origin")
dev.off()



# 单独生成umap图
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

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 100,height = 10)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis()：倾斜X轴文本
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # 增大X轴标题文本大小
        axis.title.y = element_text(size = 18))  # 增大Y轴标题文本大小
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=100*dpi,h=10*dpi,units = "px",res = dpi,type='cairo')
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

setwd("D:/R/GS/WH/20250403-8C/out(WF合并)/")
outdir <- "D:/R/GS/WH/20250403-8C/out(WF合并)/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)


# 在UMAP上展示nFeature_RNA并保存为PDF
pdf(file = paste(output, "nFeature_RNA_UMAP.pdf", sep='/'),width = 7, height = 6)
FeaturePlot(scedata, features = "nFeature_RNA", reduction = "umap", cols = c("lightgrey", "#FF3366"))
dev.off()

# 在UMAP上展示nCount_RNA并保存为PDF
pdf(file = paste(output, "nCount_RNA_UMAP.pdf", sep='/'),width = 7, height = 6)
FeaturePlot(scedata, features = "nCount_RNA", reduction = "umap", cols = c("lightgrey", '#6A4C93'))
dev.off()



# 获取TPRX1基因的表达数据
TPRX1_expression <- scedata[["RNA"]]@data["TPRX1", ]

# 创建一个新的列，定义为'celltype'，根据TPRX1基因的表达情况进行分类
scedata$celltype <- ifelse(TPRX1_expression > 0, "8CLC", "non-8CLC")
#View(scedata@meta.data)

# 绘制细胞类型的umap图
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

######查看tdt表达########
# 获取 tdTomato 基因表达数据
tdTomato_expr <- FetchData(scedata, vars = "TPRX1")

# 添加表达量信息到细胞元数据中
scedata$tdTomato_expr <- tdTomato_expr$TPRX1

# 计算 tdTomato 表达比例（非零表达的细胞比例）
expressed_cells <- sum(scedata$tdTomato_expr > 0)
total_cells <- nrow(scedata@meta.data)
expression_ratio <- expressed_cells / total_cells * 100

# 按表达量排序细胞数据
cells_ordered <- scedata@meta.data[order(scedata$tdTomato_expr, decreasing = FALSE), ]

# 提取排序后的细胞名称
cell_names_ordered <- rownames(cells_ordered)

# 设置标题，标注表达比例
plot_title <- paste0("TPRX1 Expression (", round(expression_ratio, 2), "%)")

# 绘制 FeaturePlot 按排序后的顺序绘制
pdf(paste0(output, "/TPRX1_FeaturePlot_umap.pdf"), width = 4, height = 4)
FeaturePlot(
  scedata, 
  features = "TPRX1",
  reduction = "umap", 
  cells = cell_names_ordered,  # 按细胞顺序绘制
  ncol = 1,
  cols = c('#E5D2DD',  "#FF3366")
) +
  ggtitle(plot_title) +  # 添加标题
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + 
  NoAxes()  # 删除坐标
dev.off()

svg(paste0(output, "/TPRX1_FeaturePlot_umap.svg"), width = 4, height = 4)
FeaturePlot(
  scedata, 
  features = "TPRX1",
  reduction = "umap", 
  cells = cell_names_ordered,  # 按细胞顺序绘制
  ncol = 1,
  cols = c('#E5D2DD',  "#FF3366")
) +
  ggtitle(plot_title) +  # 添加标题
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + 
  NoAxes()  # 删除坐标
dev.off()



####### 计算细胞比例 ###########
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

# 计算各组不同细胞群的细胞数
# 按样本分组计算每种细胞类型的数量
cell_counts_group <- as.data.frame(table(scedata$orig.ident,scedata$celltype))
colnames(cell_counts_group) <- c("Sample", "CellType","Counts")

cell_counts_group$CellType <- factor(cell_counts_group$CellType, levels = c( "non-8CLC","8CLC"))

# 按照"CellType"排序
cell_counts_group <- cell_counts_group %>%
  arrange(CellType)

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

file_path <- paste0(output, "/genecount_8c.pdf")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount_8c.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)


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

file_path <- paste0(output, "/geneRatio_8c.pdf")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio_8c.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)





# 定义不同细胞类型的marker基因
cellmarker <- c(
  "GLYATL3", "OR14C36", "OR7E24", "OR2M7", "PRB2",   ##  Zygote to 4 cell
  
  "ZSCAN4", "ZSCAN5A", "ZSCAN5B","ZNF280A","TPRX1", "ARGFX", "DPRX", "DUXA","DUXAB", "MDB3L2", "MBD3L3", 
  "SLC34A2", "TUBA3D",  "KDM4E", "CCNA1", "RFPL2", "MFSD2A" , "TRIM43","ZGA" ,      # 8CLC
  "DPPA3", "ALPPL2","KLF17","KLF5","MAEL","SUSD2","TRIM60","TFAP2C",    # Naive/8CLC
  "LEFTY1", "FGF4", "TCF7",     # Naive
  "CDH1","DNMT3A","FGFR4",     # Primed/naive
  "DUSP6", "SALL2","ZIC2", "FGFR1", "SALL1","TCF7L1"      # Primed
)


cellmarker <- c(
  "TBX6","MESP1","T","MIXL1","MSGN1","SIX1","EOMES","BRA","LMX1A",   ## 中胚层前体细胞（Mesodermal precursor cells）
  "BRA","EOMES","SNAI1","GSC","MESP2","MSGN1","SIX1",     ## 早期中胚层细胞（Early mesoderm cells）
  "SOX17","FOXA2","HNF4A","GATA4","SOX7","CXCR4","FOXA1","SIX4","EOMES","PDX1",   ## 内胚层细胞（Endodermal cells）
  "CK7","GATA3","TEAD4","TP63","CDX2","PLAC8","KRT18","HLA-G","PSG1","HAND1",   ## 滋养层细胞（Trophoblast Cells）
  
  "OCT4","SOX2","NANOG","KLF4","LIN28","STEM101","UTF1","TDGF1","SALL4","REX1","DPPA4",   ## 多能干细胞（Pluripotent stem cells）
  "OCT4","SOX2","NANOG","KLF4","MYC","LIN28","TDGF1","ZFP42","SALL4","DPPA4",   ## 胚胎干细胞（Embryonic Stem Cells, ESCs）
  "CD44","CD73","CD90","CD105","CD146","ENG","PDGFRB","VIM","ALCAM","THY1",   ## 间充质干细胞（Mesenchymal Stem Cells, MSCs）
  
  "PRDM1","PRDM14","NANOS3","DPPA3","DND1","PIWIL2","DAZL","VASA","STELLA","SOX17",    ## 生殖系细胞（Germ Line Cells）
  "CDH1","KRT8","KRT18","EPCAM","OCLN","CLDN1","MUC1","DSG1","JAM1","TP63",     ## 上皮细胞（Epithelial Cells） 
  "PAX6","SOX1","NES","DCX","VIM","HES1","NOTCH1","NESTIN","TLX","CD133",    ## 神经祖细胞（Neural progenitor cells）
  "MAP2","NEUROD1","NEUROD2","SYP","TUBB3","GAD1","GAP43","SNAP25","TH","CHAT",   ## 神经细胞（Neural cells）
  "TUBB3","MAP2","NEUN","SYP","DCX","NEUROD1","GAP43","SNAP25","RBFOX3","CHAT",   ## 神经元细胞（Neuronal Cells）
  "PAX6","SOX2","NES","VIM","HES1","FABP7","NOTCH1","MKI67","CDH2","GLI3"   ## 神经上皮细胞（Neuroepithelial Cells）
  
)

cellmarker <- c(
  "AIF1L", "DAG1", "DLX3", "GATA2", "GATA3", "GRAMD2A", "LMNA", "PSD4", "ST8SIA4", "TET2", "TFEB", "TINAGL1", 
  "VGLL1", "LAD1", "OVOL1", "PIEZO1", "HSPB8", "CTSL", "CDC42SE1", "ELL", "MICAL2", "FADS3", "BZW2", "WNT7A", 
  "PABPC4", "MYO1E", "CTNND1", "ARL8B", "SERPINB9", "RCOR1", "DNM2", "MKNK1", "DOT1L", "IPPK", "DDB1", "SRPRA", 
  "EPS15L1", "CAPN1", "ZC3H12A", "LIMK2", "ANXA7"
)


cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)

cellmarker <- c(
  #"CDC20","UBE2C", "CCNB1","CKS2", "HMGN2", "KPNA2","TOP2A",         ## 神经祖细胞（Neural progenitor cells）
  #"H4C3",   ## 胚胎干细胞（Embryonic Stem Cells, ESCs）
  #"DDX4",  # 精原细胞
  #"GLYATL3","OR14C36","OR7E24","OR7E24","PRB2",  #  Zygote to 4 cell
  
  "FAM32A","SERTAD1", "ATF3","ZSCAN5B", "ZNF280A","LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "DPPA3",  "KHDC1L",    # 8C
  
  "PLPP1","KLF4",    ## 桑椹胚细胞  Morula cell (Blastomere)
  "CCNB1","S100A6", "MYC",  ## 滋养层细胞（Trophoblast Cells）
  "KRT19", "KRT8",     # 滋养外胚层 Trophectoderm cell
  "FN1", "CDH11","FRZB","OTX2",         ## 原始内胚层细胞 Primitive endoderm cell  "CDH2","APOA2","FST",
  #"SOX2","GATA6",       #外胚层 Ectoderm
  #"MSGN1","TBX6",   # 中胚层 Mesoderm
  #"TFAP2A","MIXL1","EOMES","GATA2","GATA3","HAND1","T","CDX2","FOXC1","FOXC2",     ## *胚外中胚层（ExE Mesoderm）
  #"SOX17","FOXA2","HNF4A","GATA4","SOX7","CXCR4","FOXA1","SIX4","EOMES","PDX1",   ## 内胚层细胞 Endodermal cells
  "KLF17","DPPA5", "SERINC5","USP28","POU5F1",     ## 上胚层细胞  Epiblast cell
  "GATA4","GATA6","SOX17","FOXA2"       # 下胚层 Hypoblast
  #"PECAM1",  # 内皮细胞
  #"COL1A1","COL1A2","COL3A1" ,"PDGFRA"      ## 间充质干细胞（Mesenchymal Stem Cells, MSCs）
)


cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)


cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# 使用DotPlot可视化免疫细胞marker基因的表达
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

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep='/'), plot = plot, width = 14, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep='/'), plot = plot, width = 14, height = 5)


# 绘制点图
plot <- DotPlot(scedata, features = unique(cellmarker)) + RotatedAxis() +
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18))
# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_2.pdf", sep='/'), plot = plot, width = 16, height = 5)
ggsave(filename = paste(output, "marker_DotPlot_2.svg", sep='/'), plot = plot, width = 16, height = 5)


#########绘制细胞注释热图############
# change annotation color
library("scales")
library(ggsci)
library(scRNAtoolVis)

# 设定注释颜色
mycol1 <- pal_simpsons()(18)


#file_path <- file.path(outdir, "celltype.rds")
#scedata <- readRDS(file_path)

pdf(file = paste(output, "ann_Heatmap_11.pdf",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # 自定义高值颜色
dev.off()

svg(file = paste(output, "ann_Heatmap_11.svg",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # 自定义高值颜色
dev.off()



#########绘制marker基因表达箱琴图

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

# 筛选存在于数据集中的marker基因
existing_markers <- cellmarker[cellmarker %in% rownames(scedata[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# 提取表达数据并转换为适合绘图的数据格式
vln.df <- as.data.frame(scedata[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

# 继续原有步骤
anno <- scedata@meta.data[, c("CB", "seurat_clusters")]
vln.df <- inner_join(vln.df, anno, by = "CB")
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)


# 绘制Violin Plot，将X轴和Y轴调换
plot <- vln.df %>%
  ggplot(aes(exp, seurat_clusters)) +
  geom_violin(aes(fill = seurat_clusters), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # 调整facet_grid，以基因作为列
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

# 保存Violin Plot
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

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)

# 对Cluster进行细胞类型注释
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

#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom"))


# 绘制细胞类型的umap图
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
#        legend.position = c(0.99, 0.12),  # 将图例移到右下角
#        legend.justification = c("right", "bottom")) +
dev.off()


pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),width=6*length(unique(scedata$treatment)),height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=1,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
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


###########绘制分群注释点图
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

# 生成DotPlot
p <- DotPlot(scedata, features = unique(cellmarker))

# 获取数据
dat <- p$data

# 获取cluster注释并合并
anno <- distinct(data.frame(id = scedata$celltype, celltype = scedata$seurat_clusters))
colnames(anno) <- c("id", "celltype")
df <- left_join(dat, anno, by = "id")

# 定义cluster顺序
cluster.order <- c(0,2,1,3,4,5,6,7,8,9,10,11,12)
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
ggsave(plot = p, filename = paste0(output, "/ann_DotPlot.pdf"), width = 10, height = 5)
ggsave(plot = p, filename = paste0(output, "/ann_DotPlot.svg"), width = 10, height = 5)



#########绘制细胞注释热图############
# change annotation color
library("scales")
library(ggsci)
library(scRNAtoolVis)

# 设定注释颜色
mycol1 <- pal_simpsons()(18)


#file_path <- file.path(outdir, "celltype.rds")
#scedata <- readRDS(file_path)

pdf(file = paste(output, "ann_Heatmap.pdf",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # 自定义高值颜色
dev.off()

svg(file = paste(output, "ann_Heatmap.svg",sep = '/'), width = 6, height = 10)
averageHeatmap(object = scedata,
               markerGene = cellmarker)   # 自定义高值颜色
dev.off()




#########marker基因在细胞中的表达趋势#######

col <- c('#E5D2DD','#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FF3366",'#58A4C3',"#FFFFCC",'#E5D2DD',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 挑选差异细胞展示
output <- paste(outdir,'cluster', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
ScRNA <- readRDS(file_path)

cellmarker <- c(
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)

cellmarker <- c(
  "ZSCAN5B","ZNF280A","TPRX1","DUXA","DNMT3L","KLF17", "DPPA3","KHDC1L","KLF4","CD24"
)

cellmarker <- c(
  "OCT4","SOX2", "NANOG","DPPA3","KLF17","TFAP2C","ZSCAN4","DUXA","DUXB","TPRX1","ZNF280A",
  "MBD3L2","FAM151A","TRIM43","GATA6"
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

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
          legend.position = "none") +
    theme(
      axis.text.x = element_text(size = 18, color = "black", angle = 30, hjust = 1),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  
  # 获取基因表达数据，确保返回为数值向量
  gene_expr <- FetchData(ScRNA, vars = gene)[[gene]]  # 提取为数值向量
  
  # 将基因表达数据加入元数据
  ScRNA@meta.data[[paste0(gene, "_expr")]] <- gene_expr  # 添加表达量信息到元数据
  
  # 使用数值向量进行排序，避免直接索引 ScRNA 对象
  cells_ordered <- ScRNA@meta.data[order(ScRNA@meta.data[[paste0(gene, "_expr")]], 
                                         decreasing = FALSE), ]
  cell_names_ordered <- rownames(cells_ordered)  # 提取排序后的细胞名称
  
  # FeaturePlot 按排序后的顺序绘制
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", 
                                       cells = cell_names_ordered,  # 指定细胞顺序
                                       ncol = 1, cols = c('#E5D2DD', "#FF3366")) +
    theme(legend.position = "none",
          plot.title = element_text(size = 24, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, size = 1)) +   
    NoAxes()  # 删除坐标轴  
}


# 保存 RidgePlot 图
#pdf(paste0(out, "cellmarker_RidgePlot.pdf"), width = 25, height = 12)
#print(cowplot::plot_grid(plotlist = ridge_plots, ncol = 4))
#dev.off()

# 保存 FeaturePlot 图
pdf(paste(output, "cellmarker_FeaturePlot_umap.pdf",sep = '/'), width = 20, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 5))
dev.off()
svg(paste(output, "cellmarker_FeaturePlot_umap.svg",sep = '/'), width = 20, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 5))
dev.off()

pdf(paste(output, "cellmarker_VlnPlot_umap.pdf",sep = '/'), width = 20, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =5))
dev.off()
svg(paste(output, "cellmarker_VlnPlot_umap.svg",sep = '/'), width = 20, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =5))
dev.off()



#########绘制marker基因表达箱琴图

library(ggplot2)
library(reshape2)
library(dplyr)


cellmarker <- c(
  "OCT4","SOX2", "NANOG","DPPA3","KLF17","TFAP2C","ZSCAN4","DUXA","DUXB","TPRX1","ZNF280A",
  "MBD3L2","FAM151A","TRIM43","GATA6"
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# 筛选存在于数据集中的marker基因
existing_markers <- cellmarker[cellmarker %in% rownames(ScRNA[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# 提取表达数据并转换为适合绘图的数据格式
vln.df <- as.data.frame(ScRNA[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

# 继续原有步骤
anno <- ScRNA@meta.data[, c("CB", "celltype")]
vln.df <- inner_join(vln.df, anno, by = "CB")
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)


# 绘制Violin Plot，将X轴和Y轴调换
plot <- vln.df %>%
  ggplot(aes(exp, celltype)) +
  geom_violin(aes(fill = celltype), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # 调整facet_grid，以基因作为列
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

# 保存Violin Plot
ggsave(filename = paste(output, "marker_ViolinPlot_ann.pdf", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)
ggsave(filename = paste(output, "marker_ViolinPlot_ann.svg", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)



####### 计算细胞比例 ###########
col <- c('#E5D2DD','#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FF3366",'#58A4C3',"#FFFFCC",'#E5D2DD',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

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
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)


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
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 3*length(unique(scedata$orig.ident)), height = 8, dpi = 800)



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
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 20),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)

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
ggsave(file_path, plot = p2, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 4*length(unique(scedata$treatment)), height = 8, dpi = 800)




############# 8C #################

# 加载所需的包
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)

#######################Seurat分析#####################
# 设置输出目录
col <- c('#E5D2DD','#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#4F6272',"#A4CDE1",'#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


output <- paste(outdir,'8C', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype.rds")
data <- readRDS(file_path)

# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
Cells.sub <- subset(data@meta.data, celltype == c("8CLC"))
summary(Cells.sub$celltype)
scedata <- subset(data, cells=row.names(Cells.sub))


# 根据Cbr2表达量定义exp列
Cbr2_threshold <- 0  # 假设表达阈值为0，可根据实际需求调整
scedata@meta.data$exp <- ifelse(scedata@assays$RNA@data["TPRX1", ] > Cbr2_threshold, "TPRX1+", "TPRX1-")

# 检查定义结果
#View(scedata@meta.data)

library(ggsci)
# 绘制细胞类型的umap图
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

# 绘制细胞类型的umap图（SVG格式）
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

# 保存更新后的数据对象
saveRDS(scedata,"celltype(8C).rds")




####### 计算细胞比例 ###########

output <- paste(outdir,'8C', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype(8C).rds")
scedata <- readRDS(file_path)

table(scedata$exp)

# 计算各组不同细胞群的细胞数
# 按样本分组计算每种细胞类型的数量
cell_counts_group <- as.data.frame(table(scedata$orig.ident, scedata$exp))
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
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)


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
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 2.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)




############# 比较WF和SNJ的8C ############
setwd("D:/R/GS/WH/20250403-8C/out(WF合并)/8C/WF_SNJ/")
outdir <- "D:/R/GS/WH/20250403-8C/out(WF合并)/8C/WF_SNJ/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

file_path <- file.path("D:/R/GS/WH/20250403-8C/out(WF合并)/celltype(8C).rds")
WF <- readRDS(file_path)

file_path <- file.path("D:/R/GS/WH/20250403-8C/out(SNJ1)/celltype(8C).rds")
SNJ <- readRDS(file_path)


# 给每个样本添加标识
WF$treatment <- "WF"
SNJ$treatment <- "SNJ"

# 合并两个Seurat对象
merged_seurat <- merge(WF, SNJ, add.cell.ids = c("WF", "SNJ"))

# 可选：添加样本信息到metadata
merged_seurat$orig.ident <- paste0(merged_seurat$treatment, "_", merged_seurat$orig.ident)
View(merged_seurat@meta.data)

# 保存合并后的对象
saveRDS(merged_seurat, file = file.path("merged_WF_SNJ_8C.rds"))



file_path <- file.path("D:/R/GS/WH/20250403-8C/out(WF合并)/8C/WF_SNJ/merged_WF_SNJ_8C.rds")
ScRNA <- readRDS(file_path)


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

col <- c('#58A4C3','#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272',
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
genes <- c("TPRX1")
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



col <- c("#CCCCCC",'#58A4C3','#FF6666','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066")

output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)


# 根据Cbr2表达量定义exp列
Cbr2_threshold <- 0  # 假设表达阈值为0，可根据实际需求调整
scedata@meta.data$exp <- ifelse(scedata@assays$RNA@data["TPRX1", ] > Cbr2_threshold, "TPRX1+", "TPRX1-")

# 检查定义结果
#View(scedata@meta.data)

library(ggsci)
# 绘制细胞类型的umap图
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

# 绘制细胞类型的umap图（SVG格式）
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



##### 添加各组的细胞总数 #####
# 统计每个 treatment 的细胞数
cell_counts <- scedata@meta.data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(n = n()) %>%
  mutate(label = paste0(treatment, " (", n, " cells)"))

# 构建命名向量用于替换 facet 标签
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# 绘制 PDF
pdf(paste(output, "ann-diff-umap.pdf",sep = '/'),
    width = 4*length(unique(scedata$treatment)), height = 4)

DimPlot(scedata, reduction = "umap",group.by = "exp",  split.by = "treatment",
        pt.size = 1, label = FALSE, cols = col) +
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


# 保存更新后的数据对象
saveRDS(scedata,"celltype(8C).rds")




####### 计算细胞比例 ###########

output <- paste(outdir,'celltype/', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype(8C).rds")
scedata <- readRDS(file_path)

table(scedata$exp)

# 计算各组不同细胞群的细胞数
# 按样本分组计算每种细胞类型的数量
cell_counts_group <- as.data.frame(table(scedata$orig.ident, scedata$exp))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# 添加分组信息 (假设分组变量为 `treatment`)
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])  # 确保分组信息的唯一性
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# 计算每个样本中每种细胞类型的比例
cell_counts_group <- cell_counts_group %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(Ratio = Counts / sum(Counts))

p1 <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) + 
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
  #  scale_x_discrete(labels = c("WF-1", "WF-2")) +  # 修改X轴标签
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 20),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 22), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 20))  # 修改图例文本大小
# 添加细胞比例文本标签
p2 <- p2 + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 1.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 1.5*length(unique(scedata$orig.ident)), height = 6, dpi = 800)







############分组############
cell_counts_treatment <- as.data.frame(table(scedata$treatment, scedata$exp))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# 计算每个处理组中每种细胞类型的比例
cell_counts_treatment <- cell_counts_treatment %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(Ratio = Counts / sum(Counts))

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
# 添加细胞数文本标签
p1 <- p1 + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)

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
# 添加细胞数文本标签
p2 <- p2 + geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 7)


file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)
file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 3*length(unique(scedata$treatment)), height = 7, dpi = 800)







#########绘制marker基因表达箱琴图

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

# 筛选存在于数据集中的marker基因
existing_markers <- cellmarker[cellmarker %in% rownames(scedata[["RNA"]]@data)]
existing_markers <- unique(existing_markers)
# 提取表达数据并转换为适合绘图的数据格式
vln.df <- as.data.frame(scedata[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "expression")

# 删除CB列中的前缀 "WF-" 和 "SNJ-"
vln.df$CB <- gsub("^(WF_|SNJ_)", "", vln.df$CB)

# 你可以用下列命令确认前缀是否被去掉
head(vln.df$CB)

# 继续原有步骤
anno <- scedata@meta.data[, c("CB", "exp","treatment")]
vln.df <- inner_join(vln.df, anno, by = "CB")
# 设置基因和treatment顺序
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)
vln.df$treatment <- factor(vln.df$treatment, levels = c("WF", "SNJ"))  # 👈 强制设定顺序

# 绘制Violin Plot，将X轴和Y轴调换
plot <- vln.df %>%
  ggplot(aes(expression, treatment)) +
  geom_violin(aes(fill = treatment), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # 调整facet_grid，以基因作为列
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

# 保存Violin Plot
ggsave(filename = paste(output, "marker_ViolinPlot_ann.pdf", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)
ggsave(filename = paste(output, "marker_ViolinPlot_ann.svg", sep = '/'), plot = plot, width = 18,height=3,limitsize = FALSE)







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
output <- paste(outdir,'差异分析(WF vs SNJ)', sep='/')
dir.create(output, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype(8C).rds")
scRNAsub <- readRDS(file_path)
colnames(scRNAsub@meta.data)


# 寻找 Res 和 Sen 组之间的差异基因
logFCfilter <- 0.25        # 定义 log2FC 过滤值
adjPvalFilter <- 0.05   # 定义矫正后 P 值过滤值

# 寻找 Epi_cisplatin_res 和 Epi_other 组之间的差异基因
scRNAsub.cluster.markers <- FindMarkers(object = scRNAsub, 
                                        ident.1 = "WF",
                                        ident.2 =  "SNJ",
                                        group.by = "treatment", 
                                        logfc.threshold = 0, 
                                        min.pct = 0.25, 
                                        test.use = "wilcox")
scRNAsub.cluster.markers$gene <- rownames(scRNAsub.cluster.markers)

# 添加显著性标注
scRNAsub.cluster.markers <- scRNAsub.cluster.markers %>%
  mutate(Significance = ifelse(p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter, 
                               ifelse(avg_log2FC > 0, "Up", "Down"), "Normal"))
write.table(scRNAsub.cluster.markers, file = file.path(output,"sig.markers_ann_WF_vs_SNJ.txt"), sep = "\t",row.names = T, quote = FALSE)

saveRDS(scRNAsub.cluster.markers, file = file.path(output, "ScRNA.sig.markers.rds"))

# 分别保存上调基因和下调基因
upregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")
downregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")
write.csv(upregulated_genes, file = file.path(output, "upregulated_genes_WF_vs_SNJ.csv"), row.names = TRUE, quote = FALSE)
write.csv(downregulated_genes, file = file.path(output, "downregulated_genes_WF_vs_SNJ.csv"), row.names = TRUE, quote = FALSE)

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
  geom_text_repel(data = top_genes_upregulated, aes(label = top_genes_upregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  geom_text_repel(data = top_genes_downregulated, aes(label = top_genes_downregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  #geom_text_repel(data = interested_genes, aes(label = interested_genes$gene), size = 5, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  theme_classic() +
  labs(title = "WF vs SNJ", 
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
ggsave(file.path(output, "WF_vs_SNJ_volcano_plot.svg"), p, width = 8, height = 7, dpi = 300)
ggsave(file.path(output, "WF_vs_SNJ_volcano_plot.pdf"), p, width = 8, height = 7, dpi = 300)





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














