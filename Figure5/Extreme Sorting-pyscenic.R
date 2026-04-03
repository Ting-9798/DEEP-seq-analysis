

##可视化
rm(list=ls())

#BiocManager::install("RcisTarget")
#devtools::install_github("aertslab/SCopeLoomR")
devtools::install_github("aertslab/SCENIC")

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)

setwd("D:/R/GS/WH/20250901-肺癌极限分选/out/cancer cell/pyscenic/")


#######准备pyscenic的文件#########
data = readRDS("./celltype.rds")

#注意矩阵一定要转置，不然会报错
write.csv(t(as.matrix(data@assays$RNA@counts)),file = "for.scenic.data.csv")







#### 1.提取 out_SCENIC.loom 信息
loom <- open_loom('out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)



#读入测序数据
seurat.data = readRDS("./celltype.rds")
seurat.data

seurat.data@meta.data$Barcode = colnames(seurat.data)
seurat.data = subset(seurat.data, Barcode %in% colnames(regulonAUC)) # 提取抽样后的细胞
seurat.data
table(seurat.data@meta.data$Seu_Clusters)

DimPlot(seurat.data,reduction = "umap",label=T ) 

## AUC可视化
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data

#确认是否一致
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$Meta_3class))

## 将每个细胞后接分组名称
seurat.data@meta.data$Celltype_Group <- paste0(seurat.data@meta.data$Meta_3class, "_", seurat.data@meta.data$treatment)
table(seurat.data@meta.data$Celltype_Group)

cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$Meta_3class)

Celltype_Group <- data.frame(row.names = colnames(seurat.data), 
                             celltype = seurat.data$Celltype_Group)

head(cellTypes)
head(Celltype_Group)
sub_regulonAUC[1:4,1:4] 

#保存一下
save(sub_regulonAUC,cellTypes,Celltype_Group,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')



### 4.1. TF活性均值
# 看看不同单细胞亚群的转录因子活性平均值
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
#[1] 209   9
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)


Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)


### 4.2. rss查看特异TF
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)

rss_treatment <- calcRSS(AUC=getAUC(sub_regulonAUC), 
                         cellAnnotation=Celltype_Group[colnames(sub_regulonAUC), "celltype"])

#rss_treatment=na.omit(rss_treatment) 
#rssPlot <- plotRSS(rss_treatment)
#plotly::ggplotly(rssPlot$plot)



## 根据Zscore筛选celltype特异的top_n转录因子
rss_data <- rssPlot$plot$data
# 筛选正调控因子
# rss_data = subset(rss_data, grepl("\\(\\+\\)", Topic))
head(rss_data)
dim(rss_data)
library(reshape2)  
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
head(rss_data)
dim(rss_data)
colnames(rss_data)

# Zscore最大值定所属cluster
rss_summary <- data.frame(
  TF = rownames(rss_data),
  Max_Zscore = apply(rss_data, 1, max),  # 计算每一行的最大值
  MaxCluster = colnames(rss_data)[apply(rss_data, 1, which.max)]  # 找到最大值对应的列
)
# 查看结果
head(rss_summary)
table(rss_summary$MaxCluster)
rss_summary$gene <- gsub("\\(\\+\\)|\\(\\-\\)", "", rss_summary$TF)
rss_summary$regulation <- ifelse(str_detect(rss_summary$TF, "\\(\\+\\)"), "Pos", "Neg")
head(rss_summary)
write.csv(rss_summary, file = "rss_summary.csv")

top10_TFs_df = rss_summary[rss_summary$regulation == "Pos", ] # 筛选正调控转录因子
top10_TFs_df <- top10_TFs_df %>%
  # top10_TFs_df <- rss_summary %>%
  group_by(MaxCluster) %>%
  slice_max(order_by = Max_Zscore, n = 5)
head(top10_TFs_df)
table(top10_TFs_df$MaxCluster)
plot_celltypes = unique(top10_TFs_df$MaxCluster)
top10_TFs_df <- top10_TFs_df %>%
  mutate(MaxCluster = factor(MaxCluster, levels = plot_celltypes, ordered = TRUE)) %>%
  arrange(MaxCluster)
top10_TFs = top10_TFs_df[top10_TFs_df$MaxCluster %in% plot_celltypes, ]$TF
top10_TFs
# 查看RSS值
sub_rssPlot <- plotRSS(rss[top10_TFs, plot_celltypes], varName = "Cell_type", thr=0.01, 
                       col.low = '#330066', col.mid = '#66CC66', col.high = '#FFCC33')
sub_rssPlot<-sub_rssPlot$plot
ggsave("rss_dotplot.pdf",sub_rssPlot,, width = 3, height = 4)



#### 热图
rss = regulonActivity_byGroup_Scaled
head(rss)

df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat = data.frame(
                 celltype  = rownames(rss),
                 cluster = colnames(rss)[i],
                 sd.1 = rss[, i],
                 sd.2 = apply(rss[, -i], 1, median)  
               )
             }))

df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)

rowcn = data.frame(celltype = top5$cluster) 
n = rss[top5$celltype,] 

# 定义热图颜色 (例如：从蓝到红的渐变色)
heatmap_colors <- colorRampPalette(c("#4DBBD5FF", "white", "#C71000FF"))(100)

# 定义行注释的颜色 (不同细胞类型的颜色)
col40 <- c("#B22222", "#F4A460", "#FED439FF", "#91D1C2FF", "#79AF97FF", "#00A087FF", "#4DBBD5FF", "#3B4992CC", "#3C5488FF", 
           "#8968CD", "#CD96CD", "#B09C85FF", "#CD8C95", "#FD8CC1FF", "#eacb85", "#FFF68F", "#A2CD5A", "#6E8B3D", 
           "#20B2AA", "#6CA6CD", "#3A5FCD",  "#925E9FFF","#7D26CD", "#8B475D", "#008B45CC", "#631879CC", "#008280CC", 
           "#BB0021CC", "#5F559BCC", "#A20056CC", "#808180CC", "#1B1919CC", "#FF6F00FF", "#C71000FF", "#008EA0FF", 
           "#8A4198FF", "#5A9599FF", "#FF6348FF", "#84D7E1FF", "#F7B6D2CC")

# 创建celltype注释的颜色映射
# 获取唯一的celltype
unique_celltypes <- unique(rowcn$celltype)

# 创建命名颜色向量，确保颜色数量足够
annotation_colors <- list(
  celltype = setNames(col40[1:length(unique_celltypes)], unique_celltypes)
)



# 打开PDF设备，将图形保存为PDF文件
pdf("heatmap.pdf", width = 4, height = 5)
# 绘制热图，使用自定义的注释颜色
pheatmap(n,
         annotation_row = rowcn,  # 添加行注释
         show_rownames = TRUE,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,
         fontsize_row = 10,  # 增大行字体大小
         fontsize_col = 12,
         # column_title  = "Heatmap of Interaction Strength",  # 添加标题
         heatmap_legend_param = list(title = "Relative value", 
                                     title_gp = gpar(fontsize = 11, fontface = "bold", col = "black")),  # 设置图例标题
         color = heatmap_colors,
         annotation_colors = annotation_colors,  # 添加注释颜色设置
         border_color = NA  # 删除方块边框
)

# 关闭PDF设备，保存图形
dev.off()






### 点图
regulon_AUC <- regulonAUC@NAMES
seurat.data@meta.data = cbind(seurat.data@meta.data, t(assay(sub_regulonAUC[regulon_AUC,])))
# View(seurat.data@meta.data)

top10_TFs <- top10_TFs_df$TF

# 方法1：直接从metadata计算平均活性（推荐）
# 因为TF活性数据已经存储在metadata中
calculate_tf_activity_from_metadata <- function(seurat_obj, tf_list, group_var = "Meta_3class") {
  # 提取metadata中的TF活性数据
  metadata <- seurat_obj@meta.data
  tf_columns <- tf_list[tf_list %in% colnames(metadata)]
  
  if(length(tf_columns) == 0) {
    stop("没有找到指定的TF在metadata中")
  }
  
  # 按分组计算平均活性
  groups <- unique(metadata[[group_var]])
  result <- matrix(NA, nrow = length(tf_columns), ncol = length(groups))
  rownames(result) <- tf_columns
  colnames(result) <- groups
  
  for(tf in tf_columns) {
    for(group in groups) {
      group_cells <- which(metadata[[group_var]] == group)
      result[tf, group] <- mean(metadata[group_cells, tf], na.rm = TRUE)
    }
  }
  
  return(result)
}

# 计算TF活性
tf_activity <- calculate_tf_activity_from_metadata(seurat.data, top10_TFs, "Meta_3class")

# 检查结果
print("TF活性矩阵:")
print(tf_activity)

# 计算Tumor/Normal差异并排序
if("Tumor" %in% colnames(tf_activity) & "Normal" %in% colnames(tf_activity)) {
  tf_differences <- tf_activity[, "Tumor"] - tf_activity[, "Normal"]
  sorted_tfs <- names(sort(tf_differences, decreasing = TRUE))
  
  # 设置treatment组的顺序（Tumor在前）
  treatment_order <- c("Tumor", "Normal")
  seurat.data$treatment <- factor(seurat.data$treatment, levels = treatment_order)
  
  cat("TF排序信息（按Tumor-Normal差异从高到低）:\n")
  for(i in 1:length(sorted_tfs)) {
    tf <- sorted_tfs[i]
    tumor_act <- tf_activity[tf, "Tumor"]
    normal_act <- tf_activity[tf, "Normal"]
    diff <- tumor_act - normal_act
    cat(sprintf("%2d. %s: Tumor=%.3f, Normal=%.3f, Difference=%.3f\n", 
                i, tf, tumor_act, normal_act, diff))
  }
} else {
  # 如果列名不是Tumor/Normal，使用自动检测
  cat("检测到的treatment组别:\n")
  print(colnames(tf_activity))
  sorted_tfs <- top10_TFs  # 使用原始顺序
}

# 绘制点图
p_sorted <- DotPlot(seurat.data, 
                    features = sorted_tfs, 
                    group.by = 'Meta_3class') +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45,size=12),
    axis.text.y = element_text(face = "italic",size=12),  # TF名称斜体
    #legend.direction = "horizontal", 
    #legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL, 
    y = NULL,
    title = NULL,
    subtitle = NULL,
  ) +
  coord_flip() +  # 调换X和Y轴
  scale_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#C71000FF",
    midpoint = 0
  )

# 显示图形
print(p_sorted)

# 保存排序后的图形
ggsave("TF_activity_by_treatment_sorted_dotplot.pdf", 
       plot = p_sorted, 
       width = 5, 
       height = 5)






### 直接筛选treatment中各组的前五个正调控转录因子并绘制点图

# 计算每个treatment组内TF的RSS值
treatment_rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
                         cellAnnotation=Celltype_Group[colnames(sub_regulonAUC), "celltype"])

treatment_rss <- na.omit(treatment_rss)

# 提取正调控转录因子
positive_tfs <- rownames(treatment_rss)[grepl("\\(\\+\\)", rownames(treatment_rss))]
treatment_rss_pos <- treatment_rss[positive_tfs, ]

# 获取所有treatment组
treatment_groups <- unique(seurat.data@meta.data$Meta_3class)

# 筛选每个treatment组的前5个正调控TF
top_tfs_per_treatment <- list()

for(treatment_group in treatment_groups) {
  # 提取该treatment组的RSS值
  treatment_data <- treatment_rss_pos[, grep(treatment_group, colnames(treatment_rss_pos)), drop = FALSE]
  
  if(ncol(treatment_data) > 0) {
    # 计算平均RSS值（如果有多个celltype属于同一个treatment）
    if(ncol(treatment_data) > 1) {
      treatment_means <- rowMeans(treatment_data, na.rm = TRUE)
    } else {
      treatment_means <- treatment_data[, 1]
    }
    
    # 按RSS值排序并取前5个
    sorted_tfs <- names(sort(treatment_means, decreasing = TRUE))
    top_tfs <- head(sorted_tfs, 5)
    
    top_tfs_per_treatment[[treatment_group]] <- top_tfs
    
    cat("Treatment:", treatment_group, "\n")
    cat("Top 5 positive TFs:", paste(top_tfs, collapse = ", "), "\n")
    cat("RSS values:", paste(round(treatment_means[top_tfs], 3), collapse = ", "), "\n\n")
  }
}

# 合并所有treatment组的top TFs（去重）
all_top_tfs <- unique(unlist(top_tfs_per_treatment))

cat("Total unique top TFs across all treatments:", length(all_top_tfs), "\n")
print(all_top_tfs)

# 如果筛选到的TF数量不足，补充一些高表达的TF
if(length(all_top_tfs) < 5) {
  # 计算所有正调控TF的平均活性
  all_positive_tfs <- rownames(treatment_rss_pos)
  overall_means <- rowMeans(treatment_rss_pos, na.rm = TRUE)
  additional_tfs <- names(sort(overall_means, decreasing = TRUE))[1:min(10, length(overall_means))]
  additional_tfs <- setdiff(additional_tfs, all_top_tfs)
  all_top_tfs <- c(all_top_tfs, head(additional_tfs, 5 - length(all_top_tfs)))
}

# 将TF活性数据添加到seurat对象的metadata中
regulon_AUC <- regulonAUC@NAMES
seurat.data@meta.data = cbind(seurat.data@meta.data, 
                              t(assay(sub_regulonAUC[regulon_AUC,])))

# 设置treatment组的顺序（按字母顺序或自定义顺序）
treatment_order <- sort(unique(seurat.data@meta.data$treatment))
seurat.data$treatment <- factor(seurat.data$treatment, levels = treatment_order)

# 绘制点图
p_treatment <- DotPlot(seurat.data, 
                       features = all_top_tfs, 
                       group.by = 'Meta_3class') +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(hjust = 1,angle = 45,  vjust = 1, size = 12),
    axis.text.y = element_text( size = 14),  # TF名称斜体  face = "italic",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    #legend.direction = "horizontal", 
    #legend.position = "bottom",
  ) +
  labs(
    x = NULL, 
    y = NULL,
    title = NULL,
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  scale_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#C71000FF",
    midpoint = 0,
    name = "Average\nExpression"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = "Percent\nExpressed"
  )

# 显示图形
print(p_treatment)

# 保存图形
ggsave("TF_activity_by_treatment_top5_dotplot.pdf", 
       plot = p_treatment, 
       width = 8, 
       height = 2.8)



## -------------------------------
## 针对每个 Meta_3class 画该组 top5 TF 的 UMAP（AUC）
## -------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)  # 用于 plot_annotation，如果没有可以安装：install.packages("patchwork")

# 确保有 UMAP 坐标（通常 reduction = "umap"）
# 如果你的降维名字不是 "umap"，比如 "umap_scRNA"，要对应改下面 reduction 的名字

umap_plots <- list()

for (treatment_group in names(top_tfs_per_treatment)) {
  
  tf_use <- top_tfs_per_treatment[[treatment_group]]
  tf_use <- tf_use[tf_use %in% colnames(seurat.data@meta.data)]
  
  if (length(tf_use) == 0) {
    message("No matching TF columns in metadata for group: ", treatment_group)
    next
  }
  
  p_group <- FeaturePlot(
    seurat.data,
    features  = tf_use,
    reduction = "umap",
    cols      = c("lightgrey", "red"),
    order     = TRUE,
    combine   = TRUE,
    ncol      = length(tf_use)
  ) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 18)
    )
  
  umap_plots[[treatment_group]] <- p_group
  print(p_group)
  
  ggsave(
    filename = paste0("UMAP_AUC_", treatment_group, "_top5_TFs.pdf"),
    plot     = p_group,
    width    = 4 * min(length(tf_use), 5),
    height   = 4
  )
}











# 可选：分别绘制每个treatment组的点图
for(treatment_group in treatment_groups) {
  if(!is.null(top_tfs_per_treatment[[treatment_group]])) {
    treatment_tfs <- top_tfs_per_treatment[[treatment_group]]
    
    # 创建该treatment组的子集
    treatment_subset <- subset(seurat.data, treatment == treatment_group)
    
    # 如果有多个celltype，可以按celltype分组显示
    if(length(unique(treatment_subset$celltype)) > 1) {
      p_single <- DotPlot(treatment_subset, 
                          features = treatment_tfs, 
                          group.by = 'celltype') +
        theme_bw() +
        theme(
          panel.grid = element_blank(), 
          axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, size = 10),
          axis.text.y = element_text(face = "italic", size = 12),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
        ) +
        labs(
          x = NULL, 
          y = NULL,
          title = paste("Top TFs in", treatment_group)
        ) +
        scale_color_gradient2(
          low = "#4DBBD5FF",
          mid = "white",
          high = "#C71000FF",
          midpoint = 0
        )
      
      ggsave(paste0("TF_activity_", treatment_group, "_by_celltype_dotplot.pdf"), 
             plot = p_single, 
             width = 6, 
             height = 5)
    }
  }
}









## 提取TF-靶基因调控关系，调控分数
reg = read.csv("reg.csv")
head(reg)

# 先提取第一行和第二行
row1 <- as.character(reg[1, ])
row2 <- as.character(reg[2, ])

# 对于空字符串或缺失值，可以做适当处理（例如保留非空部分）
new_names <- ifelse(nchar(row1) > 0, row1, row2)

# 如果两行都有值，可以组合它们（例如用下划线连接），否则用非空的那一个
new_names <- mapply(function(x, y) {
  if(nchar(x) > 0 & nchar(y) > 0) {
    paste0(x, "_", y)
  } else if(nchar(x) > 0) {
    x
  } else {
    y
  }
}, row1, row2, USE.NAMES = FALSE)

# 将新的列名赋值给数据框
colnames(reg) <- new_names

# 删除前两行
reg <- reg[-c(1,2), ]
head(reg)
colnames(reg)
head(reg[, c("TF", "TargetGenes")])


library(stringr)

# 定义函数：传入 TF 和 TargetGenes 字符串，返回包含 TF、Gene 和 Score 的数据框
parse_target_info <- function(tf, s) {
  # 正则表达式解释：
  # \\( 匹配左括号 #'([^']+)' 捕获单引号内的基因名称 # ,\\s* 匹配逗号及其后可能的空格
  # ([0-9\\.eE+-]+) 捕获数字（支持小数和科学计数法） # \\) 匹配右括号
  pattern <- "\\('([^']+)',\\s*([0-9\\.eE+-]+)\\)"
  matches <- str_match_all(s, pattern)[[1]]
  if(nrow(matches) == 0) return(data.frame(TF = character(), Gene = character(), Score = numeric()))
  
  # 构造数据框，每个匹配对应一行
  data.frame(
    TF = rep(tf, nrow(matches)),
    Gene = matches[,2],
    Score = as.numeric(matches[,3]),
    stringsAsFactors = FALSE
  )
}
# 对 reg 数据框的每一行进行解析，并合并所有结果
TF_Genes_score <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
  parse_target_info(reg$TF[i], reg$TargetGenes[i])
}))
head(TF_Genes_score)
TF_Genes_score = unique(TF_Genes_score) # 删除重复行
dim((TF_Genes_score))
head(TF_Genes_score)
range(TF_Genes_score$Score)
table(TF_Genes_score$Score >= 1)  # 调控分数大于1的TF-gene数量
hist(TF_Genes_score$Score)

write.csv(TF_Genes_score, "All_TF_targets_Genes_score.csv", row.names = F)


# 统计每个TF的靶基因数量
TF_GeneNum = TF_Genes_score %>%
  ungroup() %>%
  group_by(TF) %>%
  summarise(Freq = n(), .groups = "drop")
TF_GeneNum <- as.data.frame(TF_GeneNum)
head(TF_GeneNum)

TF_GeneNum$TF_genenum = paste0(TF_GeneNum$TF, "(", TF_GeneNum$Freq, ")")
head(TF_GeneNum)
write.csv(TF_GeneNum, "TF_GeneNum_score.csv", row.names = F)



# 保存为新的assay
scenic_res = getAUC(sub_regulonAUC) %>% as.matrix()
scenic_res[1:5,1:5]
rownames(scenic_res)[1:5]
# 只保留正调控的TFs
scenic_res_pos <- scenic_res[grep("\\(\\+\\)", rownames(scenic_res)), ]
scenic_res_pos[1:5, 1:5]
# 提取 TF 名称，去掉 "(+)" 符号
tf_names <- gsub("\\(\\+\\)", "", rownames(scenic_res_pos))
tf_map <- setNames(TF_GeneNum$TF_genenum, TF_GeneNum$TF) # 创建 TF 名称到 TF_genenum 的映射
rownames(scenic_res_pos) <- tf_map[tf_names]
scenic_res_pos <- scenic_res_pos[!is.na(rownames(scenic_res_pos)), ]
dim(scenic_res_pos)
scenic_res_pos[1:5, 1:5]

# 添加assay
seurat.data[["scenic_pos_genenum"]] <- SeuratObject::CreateAssayObject(counts = scenic_res_pos)
seurat.data <- SeuratObject::SetAssayData(seurat.data, slot = "scale.data",
                                          new.data = scenic_res_pos, assay = "scenic_pos_genenum")
seurat.data
#保存一下
save(sub_regulonAUC,Celltype_Group,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')



# 绘制Top转录因子热图 (附带TF调控基因数量)

# 1.根据rss值最大确定TF -----------------
rss_summary = read.csv("./rss_summary.csv", row.names = 1)
rss_summary[1:5,1:5]
table(rss_summary$MaxCluster)
# 选择亚群 
TF_rss_df = rss_summary[rss_summary$MaxCluster %in% paste0("Cluster", c(1,3,6)), ]
TF_rss_df = rss_summary
head(TF_rss_df)
table(TF_rss_df$MaxCluster)

# 正调控
TF_rss_df = TF_rss_df[TF_rss_df$regulation == "Pos", ] 
# 选择前5 Zscore的TF
top5_TFs_data = TF_rss_df %>% group_by(MaxCluster) %>% top_n(5, Max_Zscore)
head(top5_TFs_data)
table(top5_TFs_data$MaxCluster)

top5_TFs_data$subcelltype = top5_TFs_data$MaxCluster
top5_TFs_data <- top5_TFs_data %>%
  arrange(subcelltype)
head(top5_TFs_data)

library(Seurat)
seurat.data
DefaultAssay(seurat.data) = "scenic_pos_genenum"
View(seurat.data@meta.data)

head(TF_GeneNum)
rownames(TF_GeneNum) = TF_GeneNum$TF
markergenes = TF_GeneNum[top5_TFs_data$gene, "TF_genenum"]
markergenes

# 使用所有包含的cluster
plot_order <- unique(top5_TFs_data$MaxCluster)
sub_seu = subset(seurat.data, celltype %in% plot_order)
table(sub_seu$celltype)

p1 <- DotPlot(sub_seu,
              features = markergenes,
              cols = c("white", "#DC050C"),
              group.by = "celltype",
              # idents = plot_order
) +
  RotatedAxis() + # 来自Seurat
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank()
  )+ coord_flip()
p1
ggsave("TFs_Clusters_Markers_DotPlot_gene_num.pdf", plot = p1, width = 6, height = 10) 



# 构建TF-Gene 调控网络
library(tidyverse)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggnetwork)

top5_TFs_data$TF_name = gsub("\\(\\+\\)", "", top5_TFs_data$TF)
head(top5_TFs_data)
TF_Genes_score= read.csv("All_TF_targets_Genes_score.csv")
head(TF_Genes_score)
hist(TF_Genes_score$Score)
dim(TF_Genes_score)
sel_TFs = c("REL",
            "TAF7",
            "MYC",
            "ARID3A",
            "ZNF274")   # "NR2F2", "ETV1", "PHF8", "ETV7", "ATF3"
#sel_TFs = top5_TFs_data$TF_name

tf_col =  c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A","#E3BE00", "#FB9A99", "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00")

# 挑选每个 TF 的前 20 个基因
sel_TF_Genes0 = TF_Genes_score[TF_Genes_score$TF %in% sel_TFs, ]
dim(sel_TF_Genes0)
hist(sel_TF_Genes0$Score)
table(sel_TF_Genes0$Score > 2)

# 排序并选择每个 TF 的前 20 个基因
sel_TF_Genes_top20 <- sel_TF_Genes0 %>%
  group_by(TF) %>%
  top_n(20, Score) %>%
  ungroup()

# 过滤掉分数大于 1 的基因
sel_TF_Genes_top20 = sel_TF_Genes_top20[sel_TF_Genes_top20$Score > 1, ]

# 确保基因是唯一的
sel_TF_Genes_top20 = unique(sel_TF_Genes_top20)

head(sel_TF_Genes_top20)

# 保存 sel_TF_Genes_top20 为 CSV 格式
write.csv(sel_TF_Genes_top20, file = "sel_TF_Genes_top20.csv", row.names = FALSE)



#构建网络
gr <- sel_TF_Genes %>% graph_from_data_frame(directed = T)
#添加一些TF、targrt信息，方便后续数据修饰
V(gr)$type = names(degree(gr))
#定义TF和targrt type
V(gr)$type[V(gr)$type %in% sel_TFs] = 'TF' 
V(gr)$type[V(gr)$type %in% sel_TF_Genes$Gene] = 'Target'
V(gr)$size[V(gr)$type =="TF"] = 2
V(gr)$size[V(gr)$type =="Target"] = 1
V(gr)$color   = 'white'
p = ggraph(gr, layout = 'sugiyama') + 
  geom_edge_link(aes(color = "#D6404E"), show.legend = F) + 
  geom_node_point(color = 'white')#网络图

p

#获取网络图数据，添加一些自己需要的信息
pData = p$data
pData = pData[rev(order(pData$type)),]
pData$color[1:length(tf_col)] <- tf_col
pData$color[pData$color == "white"] <- 'grey'


#ggplot格式作图，添加上其他内容即可
p + geom_point(data=pData,aes(x,y,color=color,size=size,stroke=1), show.legend = F) + 
  scale_color_manual(values=pData$color) + 
  geom_text(data=subset(pData, type=='Target'),aes(x,y,label=name), 
            size=3,fontface="italic", angle=45, hjust=1)+
  geom_text(data=subset(pData, type=='TF'),aes(x,y,label=name), size=4,fontface="bold")+
  theme_graph()+
  scale_y_discrete(expand=expansion(mult=c(0.5,0.05)))
p
p = ggraph(gr, layout = 'fr') + 
  geom_edge_link(aes(color = "#D6404E"), show.legend = F) + 
  geom_node_point(color = 'white')#网络图


#获取网络图数据，添加一些自己需要的信息
pData = p$data
pData = pData[rev(order(pData$type)),]
pData$color[1:length(tf_col)] <- tf_col
pData$color[pData$color == "white"] <- 'grey'


#ggplot格式作图，添加上其他内容即可
p = p + geom_point(data=pData,aes(x,y,color=color,size=size,stroke=1), show.legend = F) + 
  scale_color_manual(values=pData$color) + 
  geom_text(data=subset(pData, type=='Target'),aes(x,y,label=name), 
            size=3,fontface="italic")+
  geom_text(data=subset(pData, type=='TF'),aes(x,y,label=name), size=4,fontface="bold")+
  theme_graph()
p
ggsave("TF_graph_score_up_1.pdf", width = 16, height = 13.5)




