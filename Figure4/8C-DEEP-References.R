


############转化为csv格式###############
# 设置工作目录
setwd("D:/R/GS/WH/20250208-8C参考/data/raw/")
data_dir <- "D:/R/GS/WH/20250208-8C参考/data/raw/"

# 定义一个函数来计算RPKM
calculate_rpkm <- function(expression_data) {
  # 假设数据框中的列为细胞样本，行是基因表达
  # 获取总的转录本数（每列求和）
  total_reads <- colSums(expression_data)
  
  # 计算每千个转录本
  rpkm_data <- sweep(expression_data, 2, total_reads, FUN = "/") * 1e6
  return(rpkm_data)
}

# 遍历每个文件，读取数据并计算RPKM
for (file in file_list) {
  # 读取文件
  gene_data <- read.csv(file, header = TRUE, row.names = 1)
  
  # 计算RPKM
  rpkm_result <- calculate_rpkm(gene_data)
  
  # 可以保存RPKM结果到新的CSV文件
  output_file <- paste0("rpkm_", basename(file))
  write.csv(rpkm_result, file = output_file)
  
  # 打印当前文件的RPKM结果
  print(paste("RPKM calculated for:", file))
}



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
library(tidyr)
library(grDevices)


#install.packages("tidydr")


# 设置工作目录
setwd("D:/R/GS/WH/20250403-8C参考/data/")
# 定义数据文件路径
data_dir <- "D:/R/GS/WH/20250403-8C参考/data/"


###############Mazid et al.###############

# 定义数据文件路径
data_4CL <- "D:/R/GS/WH/20250403-8C参考/data/4CL/"

# 获取所有4CL相关的表达量文件
files <- list.files(path = data_4CL, pattern = ".*\\.counts\\.tsv\\.gz$", full.names = TRUE)

# 初始化Seurat对象列表
seurat_list <- list()

# 遍历文件，构建Seurat对象
for (file in files) {
  message("正在读取: ", basename(file))
  
  # 读取表达矩阵
  counts <- fread(file, data.table = FALSE)
  rownames(counts) <- counts[,1]   # 第一列是基因名
  counts <- counts[,-1]            # 去掉基因名列
  
  # 根据文件名添加treatment信息
  filename <- basename(file)
  treatment <- case_when(
    str_detect(filename, "^H9_4CL_D1_scRNA") ~ "4CL-D1",
    str_detect(filename, "^H9_4CL_D2") ~ "4CL-D2",
    str_detect(filename, "^H9_4CL_D3") ~ "4CL-D3",
    str_detect(filename, "^H9_4CL_D5") ~ "4CL-D5",
    str_detect(filename, "^H9_4CL_D8") ~ "4CL-D8",
    str_detect(filename, "^H9_4CL_D12_scRNA") ~ "4CL-D12",
    str_detect(filename, "^H9_e4CL_direct_D3") ~ "direct e4CL-D3",
    str_detect(filename, "^H9_e4CL_direct_D5") ~ "direct e4CL-D5",
    str_detect(filename, "^H9_e4CL_direct_D7") ~ "direct e4CL-D7",
    str_detect(filename, "^H9_e4CL_stepwise_D1") ~ "e4CL-D1",
    str_detect(filename, "^H9_e4CL_stepwise_D2") ~ "e4CL-D2",
    str_detect(filename, "^H9_e4CL_stepwise_D3") ~ "e4CL-D3",
    str_detect(filename, "^H9_e4CL_stepwise_D5") ~ "e4CL-D5",
    TRUE ~ "Unknown"
  )
  
  # 创建 Seurat 对象
  sce <- CreateSeuratObject(counts = counts,
                            project = "4CL_scRNA",
                            min.cells = 0,
                            min.features = 0)
  
  # 提取真实的 cell barcode（从最后一个下划线开始保留后面部分）
  clean_barcodes <- make.unique(sub(".*_", "", colnames(sce)))
  
  # 添加 treatment 前缀，确保唯一性
  sce <- RenameCells(sce, new.names = paste(treatment, clean_barcodes, sep = "_"))
  
  # 添加treatment信息到meta.data
  sce <- AddMetaData(sce, metadata = treatment, col.name = 'treatment')
  
  # 表达量标准化
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # 筛选表达量变化显著的基因
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  
  # 计算线粒体基因百分比（假设MT-开头）
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  
  # 计算核糖体基因百分比（假设RPS开头）
  sce[["percent.rps"]] <- PercentageFeatureSet(sce, pattern = "^RPS")
  
  # 将细胞条形码作为meta信息保存
  sce@meta.data$CB <- rownames(sce@meta.data)
  
  # 添加到列表
  seurat_list[[filename]] <- sce
}

# 合并所有Seurat对象
combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)
combined_seurat$celltype <-combined_seurat$treatment
combined_seurat@meta.data$orig.ident <- combined_seurat@meta.data$treatment

# 可视化元数据检查
#View(combined_seurat@meta.data)

# 保存合并后的Seurat对象
saveRDS(combined_seurat, file = "D:/R/GS/WH/20250403-8C参考/out/sce_female(e4CL-primed).rds")


###############Yan et al.###############

# 读取数据，并设置第一列为行名
female_data <- fread(file.path(data_dir, "rpkm_merged_expression_matrix.csv"), header = TRUE, data.table = FALSE)
rownames(female_data) <- female_data[, 1]  # 设置第一列为行名
female_data <- female_data[, -1]  # 删除第一列，因为它已经作为行名使用

# 转换为矩阵
female_data <- as.matrix(female_data)

# 查看数据的维度
dim(female_data)  # 查看行和列的数量

# 查看列名
head(colnames(female_data))

# 按照发育时期划分细胞群
treatment <- stringr::str_split(colnames(female_data), '_', simplify = TRUE)[,1]

# 查看划分的发育时期
head(treatment)

# 将列名和阶段信息进行关联
names(treatment) <- colnames(female_data)

# 查看不同发育阶段的细胞数量
table(treatment)

# 创建Seurat对象
sce_female <- CreateSeuratObject(counts = female_data,
                                 project = "sce_female",
                                 min.cells = 0,
                                 min.features = 0)


# 使用AddMetaData添加额外的元数据信息
# 将发育时期信息添加到Seurat对象中
sce_female <- AddMetaData(object = sce_female,
                          metadata = treatment,
                          col.name = 'treatment')

#### 5.表达量标准化 ####
sce_female <- NormalizeData(sce_female, normalization.method = "LogNormalize",
                            scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
sce_female <- FindVariableFeatures(sce_female, selection.method = "vst",
                                   nfeatures = 2000)

# 计算percent.mt（线粒体基因百分比）
# 假设线粒体基因的前缀为 "MT-"
sce_female[["percent.mt"]] <- PercentageFeatureSet(sce_female, pattern = "^MT-")

# 计算percent.rps（假设与核糖体相关的基因前缀为 "RPS"）
# 根据具体数据修改前缀，以下是一个假设的例子
sce_female[["percent.rps"]] <- PercentageFeatureSet(sce_female, pattern = "^RPS")

sce_female@meta.data$CB <- rownames(sce_female@meta.data)
sce_female@meta.data$orig.ident <- sce_female@meta.data$treatment

# 查看Seurat对象的元数据
#View(sce_female@meta.data)

saveRDS(sce_female, file = "D:/R/GS/WH/20250403-8C参考/out/sce_female(E0-E6).rds")

                                                                                                                                                             

############### Petropoulos et al.###############

# 读取数据，并设置第一列为行名
female_data <- fread(file.path(data_dir, "rpkm_E3-E7.csv"), header = TRUE, data.table = FALSE)
rownames(female_data) <- female_data[, 1]  # 设置第一列为行名
female_data <- female_data[, -1]  # 删除第一列，因为它已经作为行名使用

# 转换为矩阵
female_data <- as.matrix(female_data)

# 查看数据的维度
dim(female_data)  # 查看行和列的数量

# 查看列名
head(colnames(female_data))

# 按照发育时期划分细胞群
treatment <- sub("\\..*", "", colnames(female_data))

# 查看划分的发育时期
head(treatment)

# 将列名和阶段信息进行关联
names(treatment) <- colnames(female_data)

# 查看不同发育阶段的细胞数量
table(treatment)

# 为每个细胞名称添加样本前缀，确保唯一性
colnames(female_data) <- paste(treatment, colnames(female_data), sep = "_")

# 更新 treatment 信息，保持与新的细胞名一致
names(treatment) <- colnames(female_data)

# 创建Seurat对象
sce_female <- CreateSeuratObject(counts = female_data,
                                 project = "sce_female",
                                 min.cells = 0,
                                 min.features = 0)


# 使用AddMetaData添加额外的元数据信息
# 将发育时期信息添加到Seurat对象中
sce_female <- AddMetaData(object = sce_female,
                          metadata = treatment,
                          col.name = 'treatment')

#### 5.表达量标准化 ####
sce_female <- NormalizeData(sce_female, normalization.method = "LogNormalize",
                            scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
sce_female <- FindVariableFeatures(sce_female, selection.method = "vst",
                                   nfeatures = 2000)

# 计算percent.mt（线粒体基因百分比）
# 假设线粒体基因的前缀为 "MT-"
sce_female[["percent.mt"]] <- PercentageFeatureSet(sce_female, pattern = "^MT-")

# 计算percent.rps（假设与核糖体相关的基因前缀为 "RPS"）
# 根据具体数据修改前缀，以下是一个假设的例子
sce_female[["percent.rps"]] <- PercentageFeatureSet(sce_female, pattern = "^RPS")

sce_female@meta.data$CB <- rownames(sce_female@meta.data)
sce_female@meta.data$orig.ident <- sce_female@meta.data$treatment

# 查看Seurat对象的元数据
View(sce_female@meta.data)

saveRDS(sce_female, file = "D:/R/GS/WH/20250403-8C参考/out/sce_female(E3-E7).rds")


#############################整合####################################
# 读取8C和4C的Seurat对象
seurat_8C <- readRDS("D:/R/GS/WH/20250403-8C/out(WF合并)/celltype(8C).rds")
seurat_4C <- readRDS("D:/R/GS/WH/20250403-8C参考/out/sce_female(e4CL-primed).rds")
seurat_E0 <- readRDS("D:/R/GS/WH/20250403-8C参考/out/sce_female(E0-E6).rds")
#seurat_E1 <- readRDS("D:/R/GS/WH/20250403-8C参考/out/sce_female(E1-E4).rds")
seurat_E3 <- readRDS("D:/R/GS/WH/20250403-8C参考/out/sce_female(E3-E7).rds")
#seurat_E5 <- readRDS("D:/R/GS/WH/20250403-8C参考/out/sce_female(E5-E7).rds")

# 修改seurat_8C的orig.ident和treatment字段
seurat_8C$orig.ident <- "8CLC"
seurat_8C$treatment <- "8CLC"

# 从seurat_8C中提取celltype列值为8CLC的数据
seurat_8C_8CLC <- subset(seurat_8C, subset = celltype == "8CLC")
table(seurat_8C_8CLC@meta.data$celltype)
# 对Cluster进行细胞类型注释
#seurat_8C_8CLC <- RenameIdents(seurat_8C_8CLC, c("8CLC"="WF-sorted 8CLC"))

# 将细胞类型添加到meta数据中
#seurat_8C_8CLC$celltype <- seurat_8C_8CLC@active.ident
#view(seurat_8C_8CLC@meta.data)

# 赋予seurat_8C_8CLC一个新的source列，命名为"This study"
seurat_8C_8CLC$source <- "This study"
seurat_4C$source <- "Mazid et al."
seurat_E0$source <- "Yan et al."
#seurat_E1$source <- "Xue et al."
seurat_E3$source <- "Petropoulos et al."
#seurat_E5$source <- "Yanagida et al."

# 创建对象列表
seurat.list <- list(seurat_8C_8CLC, seurat_4C, seurat_E3)

# 查找整合锚点
anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)

# 进行数据整合
combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:20)    #,k.weight=20

# 归一化和降维
DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(combined_seurat, npcs = 30, verbose = FALSE)
combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:20)

# 添加celltype列，基于orig.ident来创建celltype
combined_seurat@meta.data$celltype <- combined_seurat@meta.data$treatment
combined_seurat@meta.data$seurat_clusters <- combined_seurat@meta.data$treatment

# 删除包含 NA 值的列
combined_seurat@meta.data <- combined_seurat@meta.data[, colSums(is.na(combined_seurat@meta.data)) == 0]

View(combined_seurat@meta.data)

# 可选：保存整合后的Seurat对象
saveRDS(combined_seurat, file = "D:/R/GS/WH/20250403-8C参考/out/combined_seurat_8c.rds")


#######################Seurat分析#####################
# 设置输出目录
col <- c('#FF6666',"#A4CDE1",'#FF9999',"#66CCCC",'#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         "#66CCCC",'#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#CCFFCC","#00CC66","#99FFFF", 
         "#99CCFF", '#FF6600',"#FF3366","#CC0066","#000066","#990000","#99CC66","#CC9966",
         "#FF3300","#CC99CC","#9999FF","#CCCCFF","#FF6699","#6699CC","#FFFFCC")

col <- c('#FF6666','#6A4C93','#4F6272','#E5D2DD','#58A4C3','#F9BB72', '#57C3F3', '#E59CC4',
         "#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",
         '#437eb8', "#99CCFF",'#F3B1A0', '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")



col<- c(
  # 深蓝→绿→浅绿 梯度
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # UMAP
  "#31CDEE", "#D0F199", "#79BC98", "#3C8487", "#094867",'#E59CC4',"#6666CC",
  "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF","#FEDD81", 

  # Sum-seq 浅色
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)



setwd("D:/R/GS/WH/20250403-8C参考/out/")
outdir <- "D:/R/GS/WH/20250403-8C参考/out/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# 拼接完整路径
file_path <- file.path(outdir, "combined_seurat_8c.rds")
ScRNA <- readRDS(file_path)


# 将celltype重新排序
celltype_order <- c("E3","E4","E5","E6","E7",
                    "direct e4CL-D3","direct e4CL-D5","direct e4CL-D7",
                    "4CL-D1","4CL-D2","4CL-D3","4CL-D5","4CL-D8","4CL-D12",
                    "e4CL-D1","e4CL-D2","e4CL-D3","e4CL-D5",
                    "8CLC")

#"Oocyte", "Zygote", "2 cell", "4 cell", "8 cell", "ES p0", "Morula", "Late blast", "hESC p0", "hESC p10",
#"ICM-TE-Transition","Early-TE","TE","ICM","Epiblast","Transitioning","Hypoblast"

# 更新scedata中的celltype的顺序
ScRNA$celltype <- factor(ScRNA$celltype, levels = celltype_order)

source_order <- c("This study", "Mazid et al.", "Petropoulos et al.")
ScRNA$source <- factor(ScRNA$source, levels = source_order)

View(ScRNA@meta.data)


# 提取基因表达矩阵
#gene_expression_matrix <- as.data.frame(ScRNA@assays$RNA@counts)
#output_file <- file.path(outdir, "gene_expression_matrix.csv")
#write.csv(gene_expression_matrix, file = output_file, row.names = TRUE)

# 生成小提琴图，显示质控指标
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 6)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"), ncol = 4, group.by = "source", pt.size = 0,cols = col)
dev.off()

# 生成小提琴图，显示质控指标
svg(paste(OUTPUT, "QC-BoxPlot.svg"), width = 8, height = 6)
p1 <- ggplot(data = ScRNA@meta.data, aes(x = source, y = nFeature_RNA, color = source)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nFeature_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p2 <- ggplot(data = ScRNA@meta.data, aes(x = source, y = nCount_RNA, color = source)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nCount_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

CombinePlots(plots = list(p1, p2))
dev.off()


pdf(paste(OUTPUT, "QC-ViolinPlot.pdf"), width = 12, height = 6)
# 小提琴图1：nFeature_RNA
p1 <- ggplot(data = ScRNA@meta.data, aes(x = source, y = nFeature_RNA, fill = source)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values =col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 小提琴图2：nCount_RNA
p2 <- ggplot(data = ScRNA@meta.data, aes(x = source, y = nCount_RNA, fill = source)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p3 <- ggplot(data = ScRNA@meta.data, aes(x = source, y = percent.mt, fill = source)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "percent.mt", x = "", y = "") +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p4 <- ggplot(data = ScRNA@meta.data, aes(x = source, y = percent.rps, fill = source)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "percent.rps", x = "", y = "") +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black", angle = 30, hjust = 1), 
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )


# 合并图像
library(patchwork)
(p1 | p2 | p3 |p4) + plot_layout(ncol = 4)
dev.off()


#scale_color_jco()  #20
#scale_color_npg()  #10
#scale_color_aaas()   #10
#scale_color_locuszoom()  #7
#scale_color_futurama()  #12
#scale_color_simpsons()  #16
#scale_color_rickandmorty()  #12

# 展示聚类，按Non-infected和Infected顺序展示
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 20, height = 6)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, cols = col,split.by = "source", group.by = "celltype")
dev.off()

# 单独生成umap图
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 8, height = 6)
DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE,cols = col, pt.size=0.1,group.by = "celltype") + 
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),  # 去除坐标轴刻度
        axis.title = element_blank(), # 去除坐标轴标题
        axis.ticks = element_blank(), # 去除坐标轴刻度线
        axis.line = element_blank(), # 去除坐标轴刻度线
        plot.title = element_blank(), # 去除图表标题
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
dev.off()

# 单独生成umap图
pdf(paste(OUTPUT, "cluster_umap_celltype.pdf"), width = 8, height = 6)
DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, pt.size=0.1, group.by ="celltype") + 
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        plot.title = element_blank(),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
dev.off()

# 单独生成umap图
pdf(paste(OUTPUT, "cluster_umap_source.pdf"), width = 8, height = 6)
DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, pt.size=0.1, group.by ="source") + 
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        plot.title = element_blank(),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
dev.off()

# 单独生成umap图
pdf(paste(OUTPUT, "cluster_umap_celltype_source.pdf"), width = 8, height = 6)
DimPlot(ScRNA, reduction = "umap", 
        label = FALSE, repel = TRUE, 
        pt.size=0.8,
        cols = col,
        group.by = "celltype", 
        shape.by = "source") + 
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        plot.title = element_blank(),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12)) +
  guides(shape = guide_legend(title = "Data Source"), 
         color = guide_legend(title = "Cell Type"))
dev.off()

#saveRDS(ScRNA, "ScRNA（分群后）.rds")




# 为了模仿图像结构，我们将分别绘制每个source并叠加在同一UMAP上
p1 <- DimPlot(ScRNA, reduction = "umap", group.by = "celltype", cells.highlight = list("This study" = WhichCells(ScRNA, expression = source == "This study")), cols.highlight = col, cols = "lightgrey", pt.size = 0.1) + 
  ggtitle("This study") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

p2 <- DimPlot(ScRNA, reduction = "umap", group.by = "celltype", cells.highlight = list("Mazid et al." = WhichCells(ScRNA, expression = source == "Mazid et al.")), cols.highlight = col, cols = "lightgrey", pt.size = 0.1) + 
  ggtitle("Mazid et al.") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

p3 <- DimPlot(ScRNA, reduction = "umap", group.by = "celltype", cells.highlight = list("Petropoulos et al." = WhichCells(ScRNA, expression = source == "Petropoulos et al.")), cols.highlight = col, cols = "lightgrey", pt.size = 0.1) + 
  ggtitle("Petropoulos et al.") + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


# 合并绘图
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

# 输出 PDF 文件
pdf(paste0(OUTPUT, "cluster_umap_by_source.pdf"), width = 12, height = 4.5)
print(combined_plot)
dev.off()



# 加载需要的包
library(Seurat)
library(ggplot2)
library(dplyr)

# 设置工作目录和输出路径
setwd("D:/R/GS/WH/20250403-8C参考/out/")
outdir <- "D:/R/GS/WH/20250403-8C参考/out/"
MergeOUT <- file.path(outdir, "Merge")
dir.create(MergeOUT, showWarnings = FALSE)
OUTPUT <- file.path(MergeOUT, "Multiple_")

# 加载合并后的Seurat对象
file_path <- file.path(outdir, "combined_seurat_8c.rds")
ScRNA <- readRDS(file_path)

########## Mazid ##########

Mazid <- c("direct e4CL-D3","direct e4CL-D5","direct e4CL-D7",
           "4CL-D1","4CL-D2","4CL-D3","4CL-D5","4CL-D8","4CL-D12",
           "e4CL-D1","e4CL-D2","e4CL-D3","e4CL-D5")

# 创建一个新的元数据列，用于标记选中的细胞类型和未选中的细胞
ScRNA$SelectedType <- ifelse(ScRNA$seurat_clusters %in% Mazid, 
                             ScRNA$seurat_clusters, 
                             "Unselected")

# 将 SelectedType 转换为 factor 并设置顺序
ScRNA$SelectedType <- factor(ScRNA$SelectedType, levels = c(Mazid, "Unselected"))

# 设置颜色：给Mazid指定颜色，Unselected为灰色

Mazid_colors <- colorRampPalette(c("#0571b0", "#92c5de", "#f4a582", "#ca0020"))(length(Mazid))


Mazid_colors <- colorRampPalette(c("#0571b0", "#92c5de","#66CCCC","#0099CC","#FFFFCC"))(length(Mazid))
all_colors <- c(Mazid_colors, "grey80")
names(all_colors) <- c(Mazid, "Unselected")

# 绘制UMAP图
p <- DimPlot(ScRNA, group.by = "SelectedType", 
             reduction = "umap", 
             cols = all_colors) + 
     ggtitle(" Mazid et al.") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),  # 去除坐标轴刻度
        axis.title = element_blank(), # 去除坐标轴标题
        axis.ticks = element_blank(), # 去除坐标轴刻度线
        axis.line = element_blank(), # 去除坐标轴刻度线
        plot.title = element_blank(), # 去除图表标题
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18))


# 保存为PDF
pdf(file = file.path(MergeOUT, "UMAP_Selected_Mazid.pdf"), width = 8, height = 5)
print(p)
dev.off()


########## Petropoulos ##########

Petropoulos <- c("E3","E4","E5","E6","E7")

# 创建一个新的元数据列，用于标记选中的细胞类型和未选中的细胞
ScRNA$SelectedType <- ifelse(ScRNA$seurat_clusters %in% Petropoulos, 
                             ScRNA$seurat_clusters, 
                             "Unselected")

# 将 SelectedType 转换为 factor 并设置顺序
ScRNA$SelectedType <- factor(ScRNA$SelectedType, levels = c(Petropoulos, "Unselected"))

# 设置颜色：给Mazid指定颜色，Unselected为灰色

Mazid_colors <- colorRampPalette(c("#0571b0", "#92c5de", "#f4a582", "#ca0020"))(length(Petropoulos))


Mazid_colors <- colorRampPalette(c("#ca0020", "#f4a582","#FF3366","#FFCCCC"))(length(Petropoulos))
all_colors <- c(Mazid_colors, "grey80")
names(all_colors) <- c(Petropoulos, "Unselected")

# 绘制UMAP图
p <- DimPlot(ScRNA, group.by = "SelectedType", 
             reduction = "umap", 
             cols = all_colors) + 
  ggtitle(" Petropoulos et al.") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),  # 去除坐标轴刻度
        axis.title = element_blank(), # 去除坐标轴标题
        axis.ticks = element_blank(), # 去除坐标轴刻度线
        axis.line = element_blank(), # 去除坐标轴刻度线
        plot.title = element_blank(), # 去除图表标题
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18))


# 保存为PDF
pdf(file = file.path(MergeOUT, "UMAP_Selected_Petropoulos.pdf"), width = 7.5, height = 5)
print(p)
dev.off()



########## This_study ##########
col <- c('#FF6666','#6A4C93','#4F6272','#E5D2DD','#58A4C3','#F9BB72', '#57C3F3', '#E59CC4',
         "#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",
         '#437eb8', "#99CCFF",'#F3B1A0', '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


This_study  <- c("8CLC")

# 创建一个新的元数据列，用于标记选中的细胞类型和未选中的细胞
ScRNA$SelectedType <- ifelse(ScRNA$seurat_clusters %in% This_study, 
                             ScRNA$seurat_clusters, 
                             "Unselected")

# 将 SelectedType 转换为 factor 并设置顺序
ScRNA$SelectedType <- factor(ScRNA$SelectedType, levels = c(This_study, "Unselected"))

# 设置颜色：给Mazid指定颜色，Unselected为灰色

Mazid_colors <- colorRampPalette(c("#0571b0", "#92c5de", "#f4a582", "#ca0020"))(length(This_study))


Mazid_colors <- colorRampPalette(c('#6A4C93'))(length(This_study))
all_colors <- c(Mazid_colors, "grey80")
names(all_colors) <- c(This_study, "Unselected")

# 绘制UMAP图
p <- DimPlot(ScRNA, group.by = "SelectedType", 
             reduction = "umap", 
             cols = all_colors) + 
  ggtitle(" This_study et al.") +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),  # 去除坐标轴刻度
        axis.title = element_blank(), # 去除坐标轴标题
        axis.ticks = element_blank(), # 去除坐标轴刻度线
        axis.line = element_blank(), # 去除坐标轴刻度线
        plot.title = element_blank(), # 去除图表标题
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 18))


# 保存为PDF
pdf(file = file.path(MergeOUT, "UMAP_Selected_This_study.pdf"), width = 7.5, height = 5)
print(p)
dev.off()






###########细胞手动注释########

# 设置输出目录
col <- c('#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("D:/R/GS/WH/20250403-8C参考/out/")
outdir <- "D:/R/GS/WH/20250403-8C参考/out/"

output <- paste(outdir,"celltype", sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)


# 获取TPRX1基因的表达数据
TPRX1_expression <- scedata[["RNA"]]@data["TPRX1", ]

# 创建一个新的列，定义为'celltype'，根据TPRX1基因的表达情况进行分类
scedata$celltype <- ifelse(TPRX1_expression > 0, "8CLC", "non-8CLC")
#View(scedata@meta.data)

# 绘制细胞类型的umap图
pdf(paste(output, "ann_umap_8c.pdf", sep='/'), width = 6, height = 5)
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

#saveRDS(scedata,  "celltype_8C.rds")

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



# 定义不同细胞类型的marker基因

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)

cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)


cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]


# 指定要保留的celltype
selected_celltypes <- c("8CLC","e4CL-D1","e4CL-D2","e4CL-D3","e4CL-D5",
                        "4CL-D1","4CL-D2","4CL-D3","4CL-D5","4CL-D8",
                        "4CL-D12","direct e4CL-D3","direct e4CL-D5","direct e4CL-D7")

# 筛选数据，只保留指定celltype的细胞
scedata_sub <- subset(scedata, subset = celltype %in% selected_celltypes)

# 加载绘图包
library(ggplot2)

# 绘制DotPlot，翻转X和Y轴
plot <- DotPlot(scedata_sub, features = unique(cellmarker), group.by = "celltype") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, angle = 30,size = 18),  # 调整为翻转后的坐标
    axis.text.y = element_text(size = 18),
    legend.title  = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#330066', '#336699', '#66CC66', '#FFCC33')
  ) +
  coord_flip()  # 翻转X和Y轴

# 保存图片
ggsave(filename = paste(output, "marker_DotPlot_1_8CLC.pdf", sep = '/'), plot = plot, width = 10, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1_8CLC.svg", sep = '/'), plot = plot, width = 10, height = 8)


# 绘制点图
plot <- DotPlot(scedata_sub, features = unique(cellmarker), group.by = "celltype") + RotatedAxis() +
  scale_color_gradientn(colors = c('#E5D2DD', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 20), 
        axis.title.x = element_text(size = 22), 
        axis.title.y = element_text(size = 22),
        legend.title  = element_text(size = 20),
        legend.text = element_text(size = 18))+
  coord_flip()  # 翻转坐标轴

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_2_8CLC.pdf", sep='/'), plot = plot, width = 10, height = 10)
ggsave(filename = paste(output, "marker_DotPlot_2_8CLC.svg", sep='/'), plot = plot, width = 10, height = 10)



# 指定要保留的celltype
selected_celltypes <- c("E3","E4","E5","E6","E7")

# 筛选数据，只保留指定celltype的细胞
scedata_sub <- subset(scedata, subset = celltype %in% selected_celltypes)

# 加载绘图包
library(ggplot2)

# 绘制DotPlot，翻转X和Y轴
plot <- DotPlot(scedata_sub, features = unique(cellmarker), group.by = "celltype") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, angle = 30,size = 18),  # 调整为翻转后的坐标
    axis.text.y = element_text(size = 18),
    legend.title  = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#330066', '#336699', '#66CC66', '#FFCC33')
  ) +
  coord_flip()  # 翻转X和Y轴

# 保存图片
ggsave(filename = paste(output, "marker_DotPlot_1_E3.pdf", sep = '/'), plot = plot, width = 7, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1_E3.svg", sep = '/'), plot = plot, width = 7, height = 8)


# 绘制点图
plot <- DotPlot(scedata_sub, features = unique(cellmarker), group.by = "celltype") + RotatedAxis() +
  scale_color_gradientn(colors = c('#E5D2DD', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 20), 
        axis.title.x = element_text(size = 22), 
        axis.title.y = element_text(size = 22),
        legend.title  = element_text(size = 20),
        legend.text = element_text(size = 18))+
  coord_flip()  # 翻转坐标轴

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_2_E3.pdf", sep='/'), plot = plot, width = 7, height = 10)
ggsave(filename = paste(output, "marker_DotPlot_2_E3.svg", sep='/'), plot = plot, width = 7, height = 10)

# 绘制点图
plot <- DotPlot(scedata, features = unique(cellmarker), group.by = "celltype") + RotatedAxis() +
  scale_color_gradientn(colors = c('#E5D2DD', "white", "#FF3366")) +
  theme(axis.text = element_text(size = 20), 
        axis.title.x = element_text(size = 22), 
        axis.title.y = element_text(size = 22),
        legend.title  = element_text(size = 20),
        legend.text = element_text(size = 18))+
  coord_flip()  # 翻转坐标轴

# 保存DotPlot图
ggsave(filename = paste(output, "marker_DotPlot_2.pdf", sep='/'), plot = plot, width = 16, height = 10)
ggsave(filename = paste(output, "marker_DotPlot_2.svg", sep='/'), plot = plot, width = 16, height = 10)


#########绘制marker基因表达箱琴图
library(ggplot2)
library(reshape2)
library(dplyr)

# 筛选存在于数据集中的marker基因
existing_markers <- cellmarker[cellmarker %in% rownames(scedata[["RNA"]]@data)]
existing_markers <- unique(existing_markers)

# 提取表达数据并转换为适合绘图的数据格式
vln.df <- as.data.frame(scedata[["RNA"]]@data[existing_markers,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

# 继续原有步骤
anno <- scedata@meta.data[, c("CB", "celltype")]  # 将seurat_clusters替换为celltype
vln.df <- inner_join(vln.df, anno, by = "CB")
vln.df$gene <- factor(vln.df$gene, levels = existing_markers)

# 绘制Violin Plot，将X轴和Y轴调换，并根据celltype分组
plot <- vln.df %>%
  ggplot(aes(exp, celltype)) +  # 使用celltype进行分组
  geom_violin(aes(fill = celltype), scale = "width") +
  facet_grid(. ~ gene, scales = "free_x") +  # 调整facet_grid，以基因作为列
  scale_fill_manual(values = col) +
  scale_x_continuous("") + scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 25),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),  
    strip.text = element_text(size = 18, face = "bold"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# 保存Violin Plot
ggsave(filename = paste(output, "marker_ViolinPlot.pdf", sep = '/'), plot = plot, width = 60, height = 10, limitsize = FALSE)
ggsave(filename = paste(output, "marker_ViolinPlot.svg", sep = '/'), plot = plot, width = 60, height = 10, limitsize = FALSE)



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


#########marker基因在细胞中的表达趋势#######
col <- c('#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 挑选差异细胞展示
output <- paste(outdir,'cluster', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)

cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
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
  vln_plots[[gene]] <- VlnPlot(ScRNA, features = gene, ncol = 1, pt.size = 0, cols = col, group.by = "celltype") +
    theme(axis.title.x = element_blank(),
          legend.position = "none") 
  
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
                                       ncol = 1, cols = c('#CCCCCC', "#FF3366")) +
    theme(legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, size = 1)) + 
    NoAxes()  # 删除坐标轴  
}


# 保存 RidgePlot 图
#pdf(paste0(out, "cellmarker_RidgePlot.pdf"), width = 25, height = 12)
#print(cowplot::plot_grid(plotlist = ridge_plots, ncol = 4))
#dev.off()

# 保存 FeaturePlot 图
pdf(paste(output, "cellmarker_FeaturePlot_umap.pdf",sep = '/'), width = 20, height = 20)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 5))
dev.off()
svg(paste(output, "cellmarker_FeaturePlot_umap.svg",sep = '/'), width = 20, height = 20)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 5))
dev.off()

pdf(paste(output, "cellmarker_VlnPlot_umap.pdf",sep = '/'), width = 25, height = 15)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =5))
dev.off()
svg(paste(output, "cellmarker_VlnPlot_umap.svg",sep = '/'), width = 25, height = 15)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =5))
dev.off()



cellmarker <- c(
  
  "FAM32A", "H2AFZ", "HBEGF", "ZNF23", "ZNF34", "MED26", "CDK5R1", "EPC2", "AFTPH", "TUT1", "DIO3", "GPATCH3",
  "HIST1H2BK", "HIST1H2BG", "SERTAD1", "ATF3", "ZNF266", "ZNF394", "PLAGL1", "PHC2", "ZNF337", "SLC6A16", "ZBTB16",
  "NCALD", "PRTG", "RFX4", "ZEB1", "GADD45A", "GADD45B", "SNAI1", "PRAMEF1", "ZSCAN4B", "ZSCAN5B", "ZNF280A",
  "LEUTX", "TPRX1", "DUXA", "DUXB", "DNMT3L", "KLF17", "DPPA3", "DPPA5", "KHDC1L", "POU5F1", "SOX2", "NANOG","KLF4",
  "EPCAM", "DNMT3B", "CD24", "OTX2", "CER1", "ZIC2"
)

cellmarker <- c(
  "ZNF280A","TPRX1","DUXA","KLF17", "DPPA3","CD24"
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
pdf(paste0(output, "/spacial_FeaturePlot_umap.pdf"), width = 14, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

svg(paste0(output, "/spacial_FeaturePlot_umap.svg"), width = 14, height = 8)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

pdf(paste0(output, "/spacial_VlnPlot_umap.pdf"), width = 20, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()
svg(paste0(output, "/spacial_VlnPlot_umap.svg"), width = 20, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()










#########marker基因在细胞中的表达趋势#######
col <- c('#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

library(Seurat)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(scales)

# 设置输出目录
output <- paste(outdir, '相关性', sep='/')
dir.create(output)

# 读取数据
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)

# 将celltype列添加到meta数据中
colnames(ScRNA@meta.data)
ScRNA@meta.data <- unite(ScRNA@meta.data, 
                         "celltype", 
                         celltype, 
                         remove = FALSE)

# 设置标识符为celltype
Idents(ScRNA) <- ScRNA$celltype

# 需要保留的 celltype
selected_celltypes <- c("8CLC","e4CL-D1", "e4CL-D2", "e4CL-D3", "e4CL-D5", 
                        "4CL-D1", "4CL-D2", "4CL-D3", "4CL-D5", "4CL-D8", "4CL-D12", 
                        "direct e4CL-D3", "direct e4CL-D5", "direct e4CL-D7")

selected_celltypes <- c("8CLC", "E3","E4","E5","E6","E7")

# 筛选特定 celltype
ScRNA_sub <- subset(ScRNA, idents = selected_celltypes)

# 计算平均表达谱（只针对选定的细胞类型）
exp <- AverageExpression(ScRNA_sub)$RNA
#exp <- AverageExpression(ScRNA)$RNA

# 选择表达量变化最大的2000个基因
features <- names(tail(sort(apply(exp, 1, sd)), 2000))
exp <- exp[which(row.names(exp) %in% features),]
#exp_norm <- log2(exp+1)

# 表达量归一化（每行归一化：对每个基因在不同celltype之间归一化）
#exp_norm <- t(apply(exp, 1, function(x) rescale(x, to = c(0, 1))))  # 0-1标准化

# 计算Spearman相关性矩阵
exp_cor <- cor(exp, method = "pearson")
#exp_cor <- cor(exp, method = "spearman")


# 创建热图对象并添加注释
heatmap_plot <- pheatmap(exp_cor, 
                         cluster_rows = TRUE, 
                         cluster_cols = TRUE,
                         show_rownames = TRUE, 
                         show_colnames = TRUE,
                         fontsize = 16,  # 增大字体大小
                         fontsize_row = 14,  # 增大行的文本大小
                         fontsize_col = 14,  # 增大列的文本大小
                         border_color = NA,
                         #border_color = "white",
                         main = "",
                         legend_title = "Pearson correlation",
                         color = colorRampPalette(c("#339999", "white", "#FF6666"))(100))  # 设置颜色

ggsave(file = paste(output,"Gene_Expression_Correlation_Heatmap.pdf",sep = '/'), plot = heatmap_plot$gtable, width = 10, height = 9)
ggsave(file = paste(output,"Gene_Expression_Correlation_Heatmap.svg",sep = '/'), plot = heatmap_plot$gtable, width = 10, height = 9)










######### 所有细胞的基因表达相关性趋势 #########

# 设置颜色
col <- c('#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 加载所需的库
library(Seurat)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# 设置输出目录
output <- paste(outdir, '相关性', sep='/')
dir.create(output)

# 读取数据
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)

# 确保celltype列存在于meta.data中
colnames(ScRNA@meta.data)

# 设置标识符为celltype
Idents(ScRNA) <- ScRNA$celltype

# 计算所有细胞的基因表达数据
exp <- GetAssayData(ScRNA, slot = "data")

# 将稀疏矩阵转换为普通矩阵
exp <- as.matrix(exp)

write.csv(file = paste(output, "Gene_Expression_Correlation_all_cells.csv", sep = '/'),exp)
# 计算细胞之间的相关性矩阵
exp_cor <- cor(exp, method = "pearson")

# 创建热图对象并添加注释
heatmap_plot <- pheatmap(exp_cor, 
                         cluster_rows = TRUE, 
                         cluster_cols = TRUE,
                         show_rownames = FALSE, 
                         show_colnames = FALSE,
                         fontsize = 16,  # 增大字体大小
                         fontsize_row = 14,  # 增大行的文本大小
                         fontsize_col = 14,  # 增大列的文本大小
                         border_color = NA,
                         main = "所有细胞基因表达相关性热图",
                         legend_title = "Pearson 相关性",
                         color = colorRampPalette(c("#339999", "white", "#FF6666"))(100))  # 设置颜色

# 保存热图
ggsave(file = paste(output, "Gene_Expression_Correlation_Heatmap_all_cells.pdf", sep = '/'), plot = heatmap_plot$gtable, width = 10, height = 9)
ggsave(file = paste(output, "Gene_Expression_Correlation_Heatmap_all_cells.svg", sep = '/'), plot = heatmap_plot$gtable, width = 10, height = 9)


######### 所有细胞的基因表达相关性趋势 #########

# 设置颜色
col <- c('#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 加载所需的库
library(Seurat)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# 设置输出目录
output <- paste(outdir, '相关性', sep='/')
dir.create(output)

# 读取数据
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)

# 确保celltype列存在于meta.data中
colnames(ScRNA@meta.data)

# 设置标识符为celltype
Idents(ScRNA) <- ScRNA$celltype

# 计算所有细胞的基因表达数据
exp <- GetAssayData(ScRNA, slot = "data")

# 将稀疏矩阵转换为普通矩阵
exp <- as.matrix(exp)

# 计算细胞之间的相关性矩阵
exp_cor <- cor(exp, method = "pearson")

# 获取行注释（细胞类型）
celltype_annotation <- data.frame(celltype = ScRNA$celltype)
rownames(celltype_annotation) <- colnames(exp)  # 确保注释与行匹配

# 创建热图对象并添加注释
heatmap_plot <- pheatmap(exp_cor, 
                         cluster_rows = TRUE, 
                         cluster_cols = TRUE,
                         show_rownames = FALSE, 
                         show_colnames = FALSE,
                         fontsize = 16,  # 增大字体大小
                         fontsize_row = 14,  # 增大行的文本大小
                         fontsize_col = 14,  # 增大列的文本大小
                         border_color = NA,
                         main = "所有细胞基因表达相关性热图",
                         legend_title = "Pearson 相关性",
                         color = colorRampPalette(c("#339999", "white", "#FF6666"))(100),  # 设置颜色
                         annotation_row = celltype_annotation,  # 为行添加注释（细胞类型）
                         annotation_col = celltype_annotation)  # 为列添加注释（基因类型，假设的示例）

# 保存热图
ggsave(file = paste(output, "Gene_Expression_Correlation_Heatmap_all_cells_with_annotations.pdf", sep = '/'), plot = heatmap_plot$gtable, width = 10, height = 9)
ggsave(file = paste(output, "Gene_Expression_Correlation_Heatmap_all_cells_with_annotations.svg", sep = '/'), plot = heatmap_plot$gtable, width = 10, height = 9)













#####################拟时序分析 monocle2 ###############################
library(Seurat)
library(monocle)
library(igraph)
#devtools::install_version("igraph", version = "2.0.8", repos = "http://cran.us.r-project.org")

#BiocManager::install("monocle")
#install.packages("igraph")
dpi=300

col <- c('#FF6666',"#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")
output <- paste(outdir,'monocle2', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
data <- readRDS(file_path)

#data <- subset(data, subset = treatment == "Tumor")
#View(data@meta.data)

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


col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

######## 轨迹可视化
## Pseudotime表示拟时值，State表示细胞状态，celltype表示细胞类型
## 以不同的类型进行轨迹可视化
types <- c("Pseudotime", "State", "celltype", "treatment")


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
         plot = plot_cell_traj, width = 10, height = 6)
  
  # 保存为 SVG 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 7, height = 5)
}


#saveRDS(cds,  file = file.path(output, "monocle2.rds"))
saveRDS(cds,"monocle2.rds")


#####拟时序后计算细胞比例#####

col <- c('#FF6666','#E5D2DD',"#BC8F8F","#66CCCC",'#FFCC99','#FF9999','#4F6272','#58A4C3',"#CC0066",
         '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#F9BB72',
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

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
#cds <- orderCells(cds,root_state=4)
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
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # 对于离散型数据（celltype 和 treatment），使用自定义颜色 custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 1, show_backbone = TRUE) +
      scale_color_manual(values = custom_colors) +  # 使用 custom_colors 作为颜色
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
      scale_color_manual(values = custom_colors) +  # 应用离散颜色
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





#####################拟时序分析 monocle3###############################
library(ggplot2)
library(Seurat)
library(ggsci) 
library(grid) 
library(monocle3)
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor'))
#BiocManager::install("HDF5Array")
#options(repos = c(CRAN = "https://mirrors.pku.edu.cn/CRAN/"))
#install.packages("sf")
#devtools::install_github("cole-trapnell-lab/monocle3")

dpi=300

col <- c('#E5D2DD','#FF6666',"#BC8F8F","#66CCCC",'#FFCC99','#FF9999','#4F6272','#58A4C3',"#CC0066",
         '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#99CCFF", '#3399CC',"#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",'#F9BB72',
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'monocle3', sep='/')
dir.create(output)

file_path <- file.path(outdir, "combined_seurat_8c.rds")
data <- readRDS(file_path)

#View(data@meta.data)
# 提取表达矩阵
expr_data <- GetAssayData(data, assay = 'RNA', slot = 'counts')

# 提取元数据信息
cell_metadata <- data@meta.data

# 提取基因注释信息
gene_annotation <- data.frame(gene_short_name = rownames(expr_data))
rownames(gene_annotation) <- rownames(expr_data)

# 创建 Monocle3 的 CellDataSet 对象
cds <- new_cell_data_set(
  expr_data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

# 数据预处理（标准化、降维）
cds <- preprocess_cds(cds, num_dim = 50)
pca_plot <- plot_pc_variance_explained(cds)
ggsave(paste0(output, "/PCA_variance.svg"), pca_plot, width = 6, height = 6)
ggsave(paste0(output, "/PCA_variance.pdf"), pca_plot, width = 6, height = 6)

# 降维：UMAP
cds <- reduce_dimension(cds, preprocess_method = "PCA")

# 调整 UMAP 坐标以匹配 Seurat 的降维结果
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(data, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

# 聚类分析
cds <- cluster_cells(cds)

# 绘制 UMAP 图，按 Seurat 聚类结果着色
#umap_plot <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "celltype")
umap_plot <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "celltype",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  group_label_size = 2,
  cell_size = 0.8) +
  scale_color_manual(values = col) + 
  theme_minimal(base_size = 14) + # 设置基础主题
  theme(
    panel.grid = element_blank(), # 移除网格线
    panel.background = element_blank(), # 移除背景
    axis.title = element_text(face = "bold", size = 12), # 加粗轴标题
    axis.text = element_text(size = 14), # 设置轴标签大小
    legend.title = element_blank(), # 图例标题大小
    legend.text = element_text(size = 14), # 图例文字大小
    plot.title = element_blank() # 移除图标题
  ) +
  labs(x = "UMAP_1", y = "UMAP_2")  # 设置坐标轴标题
ggsave(paste0(output, "/UMAP_clusters.svg"), umap_plot, width = 5, height = 4)
ggsave(paste0(output, "/UMAP_clusters.pdf"), umap_plot, width = 5, height = 4)


# 学习轨迹图
cds <- learn_graph(cds)

# 绘制学习轨迹图并美化
trajectory_plot <- plot_cells(cds,
                              color_cells_by = "celltype",
                              label_cell_groups = FALSE,
                              label_leaves = TRUE,
                              label_branch_points = TRUE,
                              group_label_size = 2,
                              cell_size = 0.1) +
  theme(panel.grid = element_blank(),                      # 清除背景网格
        panel.background = element_blank(),                # 设置背景为白色
        axis.line = element_line(color = "black"),         # 添加坐标轴线
        axis.title = element_text(face = "bold", size = 12), # 坐标标题加粗
        axis.text = element_text(size = 12),               # 坐标刻度字体大小
        legend.title = element_blank(),            # 图例标题字体大小
        legend.text = element_text(size = 12),             # 图例内容字体大小
        plot.title = element_blank()) +                    # 移除默认标题
  scale_color_manual(values = col) + 
  labs(x = "UMAP_1", y = "UMAP_2")  # 设置坐标轴标题

ggsave(paste0(output, "/trajectory.svg"), trajectory_plot, width = 8, height = 4)
ggsave(paste0(output, "/trajectory.pdf"), trajectory_plot, width = 8, height = 4)

# 为细胞排序并手动定义伪时间
cds <- order_cells(cds)
pseudotime_plot <- plot_cells(cds, color_cells_by = "pseudotime", 
                              label_cell_groups = FALSE, 
                              label_leaves = TRUE, 
                              label_branch_points = TRUE,
                              group_label_size = 2,
                              cell_size = 0.6) +
  theme(panel.grid = element_blank(),                      # 清除背景网格
        panel.background = element_blank(),                # 设置背景为白色
        axis.line = element_line(color = "black"),         # 添加坐标轴线
        axis.title = element_text(face = "bold", size = 12), # 坐标标题加粗
        axis.text = element_text(size = 12),               # 坐标刻度字体大小
        legend.title = element_blank(),            # 图例标题字体大小
        legend.text = element_text(size = 12),             # 图例内容字体大小
        plot.title = element_blank()) +                    # 移除默认标题
  labs(x = "UMAP_1", y = "UMAP_2")+
  scale_color_gradientn(colors = c('#6A4C93', '#FF9999', "#FFFF99"))
ggsave(paste0(output, "/pseudotime.svg"), pseudotime_plot, width = 5, height = 4)
ggsave(paste0(output, "/pseudotime.pdf"), pseudotime_plot, width = 5, height = 4)

# 手动定义根节点
#embed <- data.frame(Embeddings(data, reduction = "umap"))
#embed <- subset(embed, UMAP_1 > -6.75 & UMAP_1 < -6.5 & UMAP_2 > 0.24 & UMAP_2 < 0.25)
#root.cell <- rownames(embed)
#cds <- order_cells(cds, root_cells = root.cell)
#pseudotime_manual_plot <- plot_cells(cds, color_cells_by = "pseudotime")
#ggsave(paste0(output, "pseudotime_manual.svg"), pseudotime_manual_plot, width = 6, height = 6)
#ggsave(paste0(output, "pseudotime_manual.pdf"), pseudotime_manual_plot, width = 6, height = 6)



# 动态定义根节点函数
get_earliest_principal_node <- function(cds, time_bin = "E3") {
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  root_pr_nodes
}

# 使用函数定义根节点并排序
cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
pseudotime_dynamic_plot <- plot_cells(cds, color_cells_by = "pseudotime", 
                                      label_cell_groups = FALSE, 
                                      label_leaves = TRUE, 
                                      label_branch_points = TRUE,
                                      group_label_size = 2,
                                      cell_size = 0.6) +
  theme(panel.grid = element_blank(),                      # 清除背景网格
        panel.background = element_blank(),                # 设置背景为白色
        axis.line = element_line(color = "black"),         # 添加坐标轴线
        axis.title = element_text(face = "bold", size = 12), # 坐标标题加粗
        axis.text = element_text(size = 12),               # 坐标刻度字体大小
        legend.title = element_blank(),            # 图例标题字体大小
        legend.text = element_text(size = 12),             # 图例内容字体大小
        plot.title = element_blank()) +                    # 移除默认标题
  labs(x = "UMAP_1", y = "UMAP_2")+
  scale_color_gradientn(colors = c('#6A4C93', '#FF9999', "#FFFF99"))
ggsave(paste0(output, "/pseudotime_dynamic.svg"), pseudotime_dynamic_plot, width = 5, height = 4)
ggsave(paste0(output, "/pseudotime_dynamic.pdf"), pseudotime_dynamic_plot, width = 5, height = 4)





















####### 计算细胞比例 ###########
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 挑选差异细胞展示
output <- paste(outdir,'celltype', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)

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
cell_counts_group <- as.data.frame(table(Idents(scedata), scedata$orig.ident))
colnames(cell_counts_group) <- c("CellType", "Sample", "Counts")

# 计算每组不同细胞群的比例
cell_counts_group$Ratio <- ave(cell_counts_group$Counts, cell_counts_group$Sample, FUN = function(x) x / sum(x))


###### 在细胞比例图中添加细胞数 ######
# 更新每组的图例标签，包含细胞计数
cell_counts_group_agg <- aggregate(Counts ~ CellType, cell_counts_group, sum)
cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", cell_counts_group_agg$Counts, ")")

legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)

# 绘制数量图
p <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Counts') +
  scale_fill_manual(values = col,labels = legend_labels_group) +
  #scale_x_discrete(labels = c("WF-1", "WF-2")) +  # 修改X轴标签
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 18, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 16),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 18), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 18))  # 修改图例文本大小 

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 12, height = 8, dpi = 800)
file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 12, height = 8, dpi = 800)



# 更新每组的图例标签，包含细胞计数
cell_counts_group_agg <- aggregate(Ratio ~ CellType, cell_counts_group, sum)
cell_counts_group_agg$LegendLabel <- paste0(cell_counts_group_agg$CellType, " (", round(cell_counts_group_agg$Ratio * 100, 2), "%)")

legend_labels_group <- setNames(cell_counts_group_agg$LegendLabel, cell_counts_group_agg$CellType)

######绘制细胞比例
p <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='', y = 'Ratio') +
  scale_fill_manual(values = col,labels = legend_labels_group) +
  #scale_x_discrete(labels = c("WF-1", "WF-2")) +  # 修改X轴标签
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 18, angle = 30, hjust = 1),  # 修改X轴文本大小并旋转30度
        axis.text.y = element_text(size = 16),  # 修改Y轴文本大小
        axis.title.y = element_text(size = 18), # 修改Y轴标题大小
        legend.title = element_blank(),         # 删除图例标题
        legend.text = element_text(size = 18))  # 修改图例文本大小

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 6, height = 8, dpi = 800)
file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 6, height = 8, dpi = 800)




####8.所有细胞群进行差异分析####
# 加载必要的库
library(Seurat)
library(dplyr)
library(scRNAtoolVis)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggrepel)
output <- paste(outdir,'celltype', sep='/')
dir.create(output)
ScRNA  <- readRDS("ScRNA（分群后）.rds")

logFCfilter = 1        # 定义log2FC过滤值
adjPvalFilter = 0.05   # 定义矫正后P值过滤值
ScRNA.markers <- FindAllMarkers(object = ScRNA,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)  # 计算所有cluster marker基因

ScRNA.sig.markers <- ScRNA.markers[abs(ScRNA.markers$avg_log2FC) > logFCfilter & 
                                     ScRNA.markers$p_val_adj < adjPvalFilter, ]
write.table(ScRNA.sig.markers, file=paste(output, "sig.markers_all_cluster.txt",sep = '/'), sep="\t", row.names=F, quote=F)
head(ScRNA.sig.markers)
# 设置颜色
col <- c('#437eb8','#FF6666',"#FFFFCC",'#FFCC99','#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFCCCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")

# 绘制标准火山图
volcano_plot <- jjVolcano(diffData = ScRNA.sig.markers,
                          log2FC.cutoff = 0.25, 
                          size  = 3.5, 
                          fontface = 'italic', 
                          tile.col = col,
                          topGeneN = 10) + 
  theme(legend.title = element_text(size = 15), 
        legend.text = element_text(size = 13), 
        legend.position = c(0, 1),  
        legend.justification = c(0, 1),  
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(0.5, "lines"))

ggsave(paste(output, "markers火山图_all_cluster.pdf",sep = '/'), plot = volcano_plot, width = 12, height = 6)
ggsave(paste(output, "markers火山图_all_cluster.svg",sep = '/'), plot = volcano_plot, width = 12, height = 6)

saveRDS(ScRNA.sig.markers,  "ScRNA.sig.markers.rds")


# 进行功能注释与富集分析
# 将SYMBOL转成ENTREZID
library("org.Hs.eg.db")
## "hsa"
#library("org.Mm.eg.db")
## "mmu"
library(clusterProfiler)

ScRNA.sig.markers <- readRDS("ScRNA.sig.markers.rds")

ids <- bitr(ScRNA.sig.markers$gene, 'SYMBOL', 'ENTREZID', "org.Hs.eg.db")
ScRNA.sig.markers <- merge(ScRNA.sig.markers, ids, by.x = 'gene', by.y = 'SYMBOL')

# 按照cluster分割marker基因
gcSample <- split(ScRNA.sig.markers$ENTREZID, ScRNA.sig.markers$cluster)

# KEGG富集分析
kegg_results <- compareCluster(gcSample, fun = "enrichKEGG", organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
#kegg_results@compareClusterResult$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", kegg_results@compareClusterResult$Description)

# 保存KEGG结果为CSV文件
write.csv(as.data.frame(kegg_results), file = paste(output, "KEGG_results.csv",sep = '/'))

kegg_dotplot <- dotplot(kegg_results) + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16),  
        axis.title = element_text(size = 16, face = "bold"),  
        legend.title = element_text(size = 14, face = "bold"),  
        legend.text = element_text(size = 12)) +
  scale_fill_gradientn(colors = c("#FF3366", "white", "#00CC99")) +
  guides(color = "none")
ggsave(paste(output,"KEGG_dotplot.pdf",sep = '/'), plot = kegg_dotplot, width = 8, height = 8)
ggsave(paste(output, "KEGG_dotplot.svg",sep = '/'), plot = kegg_dotplot, width = 8, height = 8)

# GO富集分析
go_results <- compareCluster(gcSample, fun = "enrichGO", OrgDb ="org.Hs.eg.db", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.01)
# 保存GO结果为CSV文件
write.csv(as.data.frame(go_results), file = paste(output, "GO_results.csv",sep = '/'))

go_dotplot <- dotplot(go_results) + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 16), 
        axis.text.y = element_text(size = 16),  
        axis.title = element_text(size = 16, face = "bold"),  
        legend.title = element_text(size = 14, face = "bold"),  
        legend.text = element_text(size = 12)) +
  scale_fill_gradientn(colors = c("#FF3366", "white", "#00CC99")) +
  guides(color = "none")

ggsave(paste(output, "GO_dotplot.pdf",sep = '/'), plot = go_dotplot, width = 8, height = 8)
ggsave(paste(output, "GO_dotplot.svg",sep = '/'), plot = go_dotplot, width = 8, height = 8)


###############挑选有用的通路################
##############KEGG############
#selected_kegg <- c(
#  ##炎症相关通路
#  "TNF signaling pathway","NF-kappa B signaling pathway","IL-17 signaling pathway","Toll-like receptor signaling pathway",
#  "NOD-like receptor signaling pathway","Cytokine-cytokine receptor interaction","Chemokine signaling pathway",
#  "Complement and coagulation cascades","AGE-RAGE signaling pathway in diabetic complications","Arachidonic acid metabolism",
#  ##肺癌相关通路
#  "Proteoglycans in cancer","PI3K-Akt signaling pathway","Small cell lung cancer","Transcriptional misregulation in cancer",
#  "p53 signaling pathway","Cell cycle","Hippo signaling pathway","Relaxin signaling pathway","Fluid shear stress and atherosclerosis",
#  ##耐药性相关通路
#  "PI3K-Akt signaling pathway","NF-kappa B signaling pathway","MAPK signaling pathway","mTOR signaling pathway",
#  "Wnt signaling pathway","Notch signaling pathway","TGF-beta signaling pathway","Apoptosis","Autophagy",
#  ##糖酵解代谢相关通路
#  "Glycolysis / Gluconeogenesis","Pentose phosphate pathway","Pyruvate metabolism","Citrate cycle (TCA cycle)",
#  "Fructose and mannose metabolism","Galactose metabolism","Starch and sucrose metabolism",
#  "Amino sugar and nucleotide sugar metabolism","Propanoate metabolism","Butanoate metabolism",
#  ##免疫逃逸相关通路
#  "PD-L1 expression and PD-1 checkpoint pathway in cancer","Antigen processing and presentation",
#  "Natural killer cell mediated cytotoxicity","T cell receptor signaling pathway",
#  "B cell receptor signaling pathway","Fc gamma R-mediated phagocytosis","Fc epsilon RI signaling pathway","Jak-STAT signaling pathway","Th1 and Th2 cell differentiation","Th17 cell differentiation"
#  ) 
#View(kegg_results@compareClusterResult)


###############挑选有用的通路################
##############KEGG############
selected_kegg <- c("IL-17 signaling pathway","TNF signaling pathway","NF-kappa B signaling pathway",
                   "Toll-like receptor signaling pathway","NOD-like receptor signaling pathway",
                   "Chemokine signaling pathway","Cytokine-cytokine receptor interaction",
                   "Proteoglycans in cancer","PI3K-Akt signaling pathway","Focal adhesion",
                   "ECM-receptor interaction","Small cell lung cancer","MAPK signaling pathway",
                   "Transcriptional misregulation in cancer","HIF-1 signaling pathway",
                   "Glycolysis / Gluconeogenesis","PD-L1 expression and PD-1 checkpoint pathway in cancer",
                   "Antigen processing and presentation","T cell receptor signaling pathway",
                   "Natural killer cell mediated cytotoxicity")  # 替换为实际的KEGG通路名称
#View(kegg_results@compareClusterResult)
# 筛选指定KEGG通路
kegg_results_filtered <- subset(as.data.frame(kegg_results@compareClusterResult), Description %in% selected_kegg)

# 绘制筛选后的 KEGG 富集点图
# 绘制筛选后的 KEGG 富集点图
filtered_kegg_dotplot <- ggplot(kegg_results_filtered, aes(x = Cluster, y = Description, size = Count, color =p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(colors = c( "#D73027","#91BFDB")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  labs(title = "",y = "",x = "Cluster",size = "Gene Count")
# 保存结果
ggsave(paste(output, "KEGG_dotplot_filtered.pdf",sep = '/'), plot = filtered_kegg_dotplot, width = 12, height = 8)
ggsave(paste(output, "KEGG_dotplot_filtered.svg",sep = '/'), plot = filtered_kegg_dotplot, width = 12, height = 8)


############GO#############
#selected_go_terms <- c("GO:0008150", "GO:0003674", "GO:0005575") # 替换为实际的GO术语ID或描述
# 筛选指定GO术语
#go_filtered_results <- subset(as.data.frame(go_results), Description %in% selected_go_terms)


