


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

#######################Seurat分析#####################
# 设置输出目录
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/cancer cell/")
outdir <- "D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/cancer cell/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

data <- readRDS("D:/R/GS/WH/20250716-超高通/out(与卵巢癌合并)/celltype.rds")

# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
ScRNA <- subset(data, idents = c("Cancer Cells"))
table(ScRNA@meta.data$celltype)


# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
#Cells.sub <- subset(data@meta.data, celltype == c("Proximal Tubule","TAL","DCT","IC","Collecting duct Principal"))
#summary(Cells.sub$celltype)
#ScRNA <- subset(data, cells=row.names(Cells.sub))
#View(ScRNA@meta.data)


#### 6.归一化与PCA降维 ####

# 数据标准化与高变基因寻找
ScRNA <- NormalizeData(ScRNA)
ScRNA <- FindVariableFeatures(ScRNA, selection.method = "vst", nfeatures = 2000)

#归一化
ScRNA<-ScaleData(ScRNA)

#运行PCA
ScRNA<-RunPCA(ScRNA,npcs = 30)

pdf(paste(OUTPUT,"Dimplot.pdf"),width = 9,height = 6)
p1 <- DimPlot(object = ScRNA, reduction = "pca", pt.size = .1, group.by = "treatment",cols = col)
CombinePlots(plots=list(p1))
dev.off()

pdf(paste(OUTPUT,"vlnplot.pdf"),width = 9,height = 6)
p2 <- VlnPlot(object = ScRNA, features = "PC_1", group.by = "treatment", pt.size = 0,cols = col)
CombinePlots(plots=list(p2))
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
ScRNA<-RunHarmony(ScRNA,group.by.vars = c("treatment"), npcs = 30,
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

#Idents(ScRNA) <- "integrated_snn_res.1"
Idents(ScRNA) <- "RNA_snn_res.1"
ScRNA$seurat_clusters <- ScRNA@active.ident##根据聚类树选择你要的resolution
table(Idents(ScRNA))

# 确保 "treatment" 因子水平按照 Non-infected 和 Infected 顺序排列
#ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# 展示聚类，按Non-infected和Infected顺序展示
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
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))

DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
dev.off()


# 单独生成umap图
pdf(paste(OUTPUT, "cluster_umap1.pdf"), width = 5.5, height = 4)
DimPlot(ScRNA, reduction = "umap", label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03,size = 14),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16))
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

saveRDS(ScRNA, "ScRNA（分群后）.rds")




############绘制基因表达热图###########
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)

col <- c("#99CCFF","#FF3366","#66CCCC","#FF9933","#CC0066")

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
ScRNA <- readRDS("ScRNA（分群后）.rds")

# 设置分组顺序
ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Ovarian cancer","Gastric cancer","Pancreatic cancer","WF-14W"))

# 设置待分析的T细胞激活相关基因
genes <- c(
  #"EPCAM","WFDC2","KRT8","CLDN4","CD24","MKI67","TOP2A",      # Cancer Cells  "TP53",
  "WT1","MUC16","PAX8",  #"ESR1","PGR",   ## 卵巢
  "VIL1","GCC1","EPCAM",  #"SPOCK1","APOC1","HER2",    ## 胃癌
  "MUC1","FUT3","CEACAM5","CEACAM6","SPAN1","TIMP1",        ## 胰腺癌
  "KRT20","CDX2","MUC6"   ## 消化道
  
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
  col = list(Treatment = c("WF-14W" = "#99CCFF","Pancreatic cancer" =  "#FF3366","Gastric cancer" = "#66CCCC", "Ovarian cancer" = "#FF9933", "14d" = "#CC0066"))
)

# 保存热图（PDF格式）
pdf(file.path(output, "cellmarker_Heatmap_byTreatment.pdf"), width = 6, height = 4)
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




############绘制基因表达uamp###########
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)


output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
ScRNA <- readRDS("ScRNA（分群后）.rds")

# 设置分组顺序
ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Ovarian cancer","Gastric cancer","Pancreatic cancer","WF-14W"))


#####12.展示已知的细胞marker基因的表达情况####
# 定义要绘制的基因列表
genes <- c(
  "WT1","MUC16","PAX8",  #"ESR1","PGR",   ## 卵巢
  "VIL1","GCC1","EPCAM",  #"SPOCK1","APOC1","HER2",    ## 胃癌
  "MUC1","FUT3","CEACAM5","CEACAM6","TIMP1",        ## 胰腺癌
  "KRT20","CDX2","MUC6"   ## 消化道
  
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
  
  ggsave(filename = file.path(output, paste0("FeaturePlot_", gene, "(不同组).pdf")), plot = p1, device = "pdf", width = 6, height = 6)
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


sco <- readRDS("ScRNA（分群后）.rds")
data <- sco@meta.data
colnames(data)


#构建Roe计算需要的输入表格
#data <- data[,c(1,4,7,11)]
#colnames(data) <- c("sample","tissue","celltype")


# ---- 计算 Roe 并准备矩阵 ----
Roe <- calTissueDist(
  data,
  byPatient = FALSE,
  colname.cluster = "seurat_clusters",
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
output <- paste(outdir,'差异分析(转移和非转移)', sep='/')
dir.create(output, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)
View(scRNAsub@meta.data)



logFCfilter = 0.25        # 定义log2FC过滤值
adjPvalFilter = 0.05   # 定义矫正后P值过滤值
# 选择要分析的细胞类型
selected_celltypes <- c("Metastatic_high","Metastatic_low")

# 筛选出属于选定细胞类型的数据
ScRNA_selected <- subset(ScRNA, Metastatic_3class %in% selected_celltypes)

# ✅ 将 celltype 设置为有序因子（按照 selected_celltypes 顺序）
ScRNA_selected$Metastatic_3class <- factor(ScRNA_selected$Metastatic_3class, levels = selected_celltypes)

# 定义logFC过滤值和矫正后P值过滤值
logFCfilter = 1
adjPvalFilter = 0.05

#Idents(ScRNA_selected) <- 'Metastatic_3class'



# 寻找 Res 和 Sen 组之间的差异基因
logFCfilter <- 0.25        # 定义 log2FC 过滤值
adjPvalFilter <- 0.05   # 定义矫正后 P 值过滤值

# 寻找 Epi_cisplatin_res 和 Epi_other 组之间的差异基因
scRNAsub.cluster.markers <- FindMarkers(object = ScRNA_selected, 
                                        ident.1 = "Metastatic_high",
                                        ident.2 =  "Metastatic_low",
                                        group.by = "Metastatic_3class", 
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
  
  "S100A4","NEDD9","FN1","CD44","MMP19","LOXL2","LOX","SPP1","SERPINE1",
  "VIM","COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","COL5A2","COL6A3",
  "VCAN","THBS1","SPARC","TIMP2","RHOB","RND3","CXCL8","CXCL1","CXCL2","CXCL3",
  "CCL2","CCL3","CCL4","PDGFRA","PDGFRB","VEGFA","S100A10","S100A11",
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
  labs(title = "WF-14W vs Tumor", 
       x = "log2 Fold Change", y = "-log10 Adjusted P-value", color = "Significance") +
  #annotate("text", x = -1.5, y = 300, label = paste("Up-regulated:", upregulated_genes), 
  #         hjust = 0, size = 5, color = "#FF3300", fontface = "bold") +
  #annotate("text", x = -1.5, y = 285, label = paste("Down-regulated:", downregulated_genes), 
  #         hjust = 0, size = 5, color = "#6699CC", fontface = "bold") +
  #annotate("text", x = -1.5, y = 270, label = paste("DEGs:", total_diff_genes), 
  #         hjust = 0, size = 5, color = "black", fontface = "bold") +
  scale_x_continuous(
    limits = c(-2, 2),                      # 设置 X 轴范围
    breaks = seq(-1.5, 2, by = 1),            # 设置刻度
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
  xlim(-0.5, 0.6) +
  ylim(-3, 3) +
  labs(title = "Tumor vs Normal")  # 设置标题和子标题

# 保存为PDF和SVG
ggsave(file.path(output, "volcano_plot.pdf"), volcano_plot, width = 8, height = 6)
ggsave(file.path(output, "volcano_plot.svg"), volcano_plot, width = 8, height = 6)





###############差异基因在不同组中表达##################
library(Seurat)
library(tidyverse)
library(ggsci)


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
  "Focal adhesion", "ECM-receptor interaction", "TGF-beta signaling pathway", 
  "Rap1 signaling pathway", 
  "PI3K-Akt signaling pathway", "MAPK signaling pathway", "NF-kappa B signaling pathway", "HIF-1 signaling pathway", 
  "FoxO signaling pathway", "Apoptosis", "Autophagy - animal", "Hippo signaling pathway", "p53 signaling pathway", 
  "Ras signaling pathway", "Chemokine signaling pathway", "TNF signaling pathway", "IL-17 signaling pathway", 
  "Toll-like receptor signaling pathway", "NOD-like receptor signaling pathway"
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
ggsave(file.path(output, 'kegg_enrich_up_bar.pdf'), plot = kegg_up_plot, width = 10, height = 6)
#ggsave(file.path(output, 'kegg_enrich_down_bar.pdf'), plot = kegg_down_plot, width = 8, height = 5)

# 组合并保存
#combined_GO_plot <- go_up_plot + go_down_plot + plot_layout(guides = 'collect')
#combined_KEGG_plot <- kegg_up_plot + kegg_down_plot + plot_layout(guides = 'collect')

#ggsave(file.path(output, 'combined_GO_bar.pdf'), plot = combined_GO_plot, width = 13, height = 10)
#ggsave(file.path(output, 'combined_KEGG_bar.pdf'), plot = combined_KEGG_plot, width = 13, height = 10)





#### 划分转移和非转移细胞 ########
#BiocManager::install("UCell")
library(UCell)
library(Seurat)

# 设置输出目录
output <- paste(outdir, 'Ucell', sep='/')
dir.create(output, showWarnings = FALSE)

# 1. 读取 Seurat 数据
sce <- readRDS("ScRNA（分群后）.rds")

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
view(sce@meta.data)
a <- colnames(sce@meta.data) %>% str_subset("UCell")

# 4. 可视化
p <- FeaturePlot(
  sce,
  features = a,
  order = TRUE,
  ncol = 2,
  cols = viridis(256)
) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    title = "胰腺癌转移相关基因UCell评分",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "UCell分数"
  )

# 5. 保存 PDF
ggsave(
  filename = file.path(output, "Meta_FeaturePlot.pdf"),plot = p,width = 4,height = 3
)











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
sce <- readRDS("ScRNA（分群后）.rds")

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
  celltype_col <- "seurat_clusters"
  
  # 从 meta.data 抽出要画图的数据
  df_plot <- data.frame(
    seurat_clusters = sce@meta.data[[celltype_col]],
    score    = sce@meta.data[[score_col]],
    stringsAsFactors = FALSE
  )
  
  # 固定 celltype 的顺序（按字母排序或自己指定顺序）
  df_plot$seurat_clusters <- factor(df_plot$seurat_clusters,
                             levels = sort(unique(df_plot$seurat_clusters)))
  
  p_box <- ggplot(df_plot, aes(x = seurat_clusters, y = score, fill = seurat_clusters)) +
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
  
  
  ggsave(file.path(out_dir, paste0(set_name, "_VlnPlot_seurat_clusters.pdf")),
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
ScRNA  <- readRDS("celltype.rds")
table(ScRNA$Metastatic_3class)



logFCfilter = 0.25        # 定义log2FC过滤值
adjPvalFilter = 0.05   # 定义矫正后P值过滤值
# 选择要分析的细胞类型
selected_celltypes <- c("Metastatic_high","Metastatic_low")

# 筛选出属于选定细胞类型的数据
ScRNA_selected <- subset(ScRNA, Metastatic_3class %in% selected_celltypes)

# ✅ 将 celltype 设置为有序因子（按照 selected_celltypes 顺序）
ScRNA_selected$Metastatic_3class <- factor(ScRNA_selected$Metastatic_3class, levels = selected_celltypes)

# 定义logFC过滤值和矫正后P值过滤值
logFCfilter = 1
adjPvalFilter = 0.05

Idents(ScRNA_selected) <- 'Metastatic_3class'


# 计算所有cluster marker基因
ScRNA.markers <- FindAllMarkers(object = ScRNA_selected,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)  # 计算所有cluster marker基因

ScRNA.sig.markers <- ScRNA.markers[abs(ScRNA.markers$avg_log2FC) > logFCfilter & 
                                     ScRNA.markers$p_val_adj < adjPvalFilter, ]
write.table(ScRNA.sig.markers, file=paste(output, "sig.markers_all_cluster.txt",sep = '/'), sep="\t", row.names=F, quote=F)
head(ScRNA.sig.markers)


saveRDS(ScRNA.sig.markers,  "ScRNA.sig.markers.rds")



# 设置颜色
col <- c('#437eb8','#FF6666',"#FFFFCC",'#FFCC99','#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFCCCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC",
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300","#FFCCCC",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC")

# 确保 cluster 按照 selected_celltypes 排序
ScRNA.sig.markers$cluster <- factor(ScRNA.sig.markers$cluster, levels = selected_celltypes)

# ✅ 设置感兴趣的基因列表
interest_genes <- c(
  ## "FAP+ THY1+ CAFs"
  "INHBA","ID3","NBL1","DCN","LTBP1","BMPR2","SKP1","ID1","SKIL","SMAD1","NCOR1","CUL1",
  "SFPR2","SERPINF1","SFRP4","FRZB","NKD2","APCDD1","JUN","CACYBP","TCF7L2","PRICKLE1","SKP1","CAMK2D","CSNK2B","DAAM1","FZD1","TLE4","ROR2","CUL1","CSNK1A1",
  "ACTB","PARD3","NKD2","STK3","PPP1CA","DLG2","YWHAB","TCF7L2","BMPR2","ID1","SMAD1","FRMD6","WWTR1","YWHAQ","FZD1","TEAD1",
  
  ## "iCAFs"
  "PTGS2","IL1B","CXCL8","CD14","CXCL12","BCL2A1","CCL4L2","CCL4","IL1R1","CXCL3","CXCL2","CXCL1","LYN","TNFAIP3","VCAM1","NFKB2"
  
)



# =========================
# ✅ 自定义 multi-volcano 风格绘图
# =========================

# 让 cluster 按 selected_celltypes 顺序
ScRNA.sig.markers$cluster <- factor(ScRNA.sig.markers$cluster, levels = selected_celltypes)

# 1) 背景柱状图需要每个cluster的 logFC 最小/最大值
dbar <- ScRNA.sig.markers %>%
  group_by(cluster) %>%
  summarize(
    logFC_min = min(avg_log2FC),
    logFC_max = max(avg_log2FC),
    .groups = "drop"
  ) %>%
  rename(group = cluster)

# 2) 所有差异基因点（用于展示全部点）
AllGene <- ScRNA.sig.markers %>%
  mutate(
    group = cluster,
    logFC = avg_log2FC,
    label = ifelse(logFC > 0, "Up", "Down")
  )

# 3) 构建用于打标签的基因集
# ---- 上调：只显示感兴趣基因
TopGene_up <- AllGene %>%
  filter(logFC > 0, gene %in% interest_genes)

# ---- 下调：显示 top10（按 logFC 最小排序）
TopGene_down <- AllGene %>%
  filter(logFC < 0) %>%
  arrange(logFC) %>%        # logFC 越小越靠前（更下调）
  slice_head(n = 15)

# 合并为最终打标签的数据
TopGene <- bind_rows(TopGene_up, TopGene_down)

# （可选）提示：感兴趣上调基因缺失
miss_up_genes <- setdiff(
  interest_genes,
  TopGene_up$gene
)
if(length(miss_up_genes) > 0){
  message("⚠️ 以下感兴趣基因未出现在上调显著集中（不显示标签）：\n",
          paste(miss_up_genes, collapse = ", "))
}

# 4) 绘图（只需要把 repel 层用 TopGene）
p <- ggplot() +
  geom_col(
    data = dbar,
    aes(x = group, y = logFC_min),
    fill = "#dcdcdc", alpha = 0.6, width = 0.7
  ) +
  geom_col(
    data = dbar,
    aes(x = group, y = logFC_max),
    fill = "#dcdcdc", alpha = 0.6, width = 0.7
  ) +
  
  # 所有差异基因点
  geom_jitter(
    data = AllGene,
    aes(x = group, y = logFC, color = label),
    size = 1.5, width = 0.35, alpha = 0.9
  ) +
  
  geom_tile(
    data = dbar,
    aes(x = group, y = 0, fill = group),
    height = 0.6, color = "black", alpha = 0.6, show.legend = FALSE
  ) +
  
  # ✅ 标签：上调兴趣基因 + 下调top10
  geom_text_repel(
    data = TopGene,
    aes(x = group, y = logFC, label = gene),
    size = 5, color = "black",
    force = 1.2, max.overlaps = 50
  ) +
  
  geom_text(
    data = dbar,
    aes(x = group, y = 0, label = group),
    size = 5, color = "black"
  ) +
  
  ggsci::scale_fill_npg() +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5")) +
  labs(x = "", y = "log2 Fold Change", color = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20, color = "black", face = "bold"),
    axis.text = element_text(size = 18, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 1),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 18)
  )


# 保存图像
ggsave(filename = file.path(output, "multi_volcano_plot.pdf"),
       plot = p, width = 6, height = 6, bg = "white")
ggsave(filename = file.path(output, "multi_volcano_plot.svg"),
       plot = p, width = 6, height = 6, bg = "white")




# 绘制标准火山图
volcano_plot <- jjVolcano(diffData = ScRNA.sig.markers,
                          log2FC.cutoff = 0.25, 
                          size  = 4, 
                          fontface = 'italic', 
                          tile.col = col,
                          topGeneN = 10) + 
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.position = c(2, 0),  
        legend.justification = c(2, 0),
        axis.title = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(size = 16),
        panel.spacing = unit(0.5, "lines"))

ggsave(paste(output, "markers火山图_all_cluster.pdf",sep = '/'), plot = volcano_plot, width = 5, height = 5)
ggsave(paste(output, "markers火山图_all_cluster.svg",sep = '/'), plot = volcano_plot, width = 5, height = 5)





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
ggsave(paste(output,"KEGG_dotplot.pdf",sep = '/'), plot = kegg_dotplot, width = 16, height = 25)
ggsave(paste(output, "KEGG_dotplot.svg",sep = '/'), plot = kegg_dotplot, width = 16, height = 25)

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

ggsave(paste(output, "GO_dotplot.pdf",sep = '/'), plot = go_dotplot, width = 14, height = 25)
ggsave(paste(output, "GO_dotplot.svg",sep = '/'), plot = go_dotplot, width = 14, height = 25)


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
                   "ECM-receptor interaction","MAPK signaling pathway",
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
ggsave(paste(output, "KEGG_dotplot_filtered.pdf",sep = '/'), plot = filtered_kegg_dotplot, width = 6, height = 6)
ggsave(paste(output, "KEGG_dotplot_filtered.svg",sep = '/'), plot = filtered_kegg_dotplot, width = 6, height = 6)








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
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
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
p1 <- FeaturePlot(sce, features = "MAPK", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("MAPK") +
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
p4 <- FeaturePlot(sce, features = "TNFa", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis::turbo(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("TNFa") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

# 绘制MAPK通路
p5 <- FeaturePlot(sce, features = "NFkB", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis::turbo(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("NFkB") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

# 绘制MAPK通路
p6 <- FeaturePlot(sce, features = "Hypoxia", 
                  coord.fixed = TRUE, 
                  order = TRUE, 
                  cols = viridis::turbo(10)) +
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  ggtitle("Hypoxia") +
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
print(p6)
dev.off()



# 单独保存Hypoxia通路UMAP图
pdf(file.path(output, "MAPK_pathway_UMAP.pdf"), width = 5, height = 4)
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
pdf(file.path(output, "TNFa_pathway_UMAP.pdf"), width = 5, height = 4)
print(p4)
dev.off()

# 单独保存TNFa通路UMAP图
pdf(file.path(output, "NFkB_pathway_UMAP.pdf"), width = 5, height = 4)
print(p5)
dev.off()

# 单独保存TNFa通路UMAP图
pdf(file.path(output, "Hypoxia_pathway_UMAP.pdf"), width = 5, height = 4)
print(p6)
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
pathways_to_plot <- c("MAPK", "TGFb", "EGFR", "Hypoxia", "NFkB","TNFa")

# 指定 CellType 顺序
#celltype_order <- c("Non-malignant cells", "Malignant cells")

# 生成颜色（与顺序一致）
my_colors <- c(
  "Tumor" = "#31CDEE", 
  "WF-14W"     = "#D0F199"
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
      comparisons = list(c("WF-14W", "Tumor")),
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











####################### 单个基因相关性分析 (treatment 分组) ##############################

# 加载必要的R包
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 设置输出目录
output <- paste(outdir, 'Correlation_treatment', sep='/')
dir.create(output, showWarnings = FALSE)

# 读取数据
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scRNAsub <- readRDS(file_path)

# 定义每个基因的颜色
gene_colors <- c(
  "PAX8" = "#FF3366",        # 卵巢癌
  "VIL1" = "#1E90FF",        # 胃癌
  "TIMP1" = "#FF6666",       # 胰腺癌
  "FUT3" = "#33CCCC",         # 胰腺癌
  "CEACAM5"="#993399"   # 胰腺癌
)

# 定义每个癌症类型相关的基因列表
ovarian_cancer_genes <- "PAX8"        # 卵巢癌
gastric_cancer_genes <- "VIL1"        # 胃癌
pancreatic_cancer_genes <- c("TIMP1", "FUT3","CEACAM5")  # 胰腺癌

# 提取基因表达矩阵和 treatment 信息
expr_data <- as.data.frame(t(scRNAsub@assays$RNA@data[c(ovarian_cancer_genes, gastric_cancer_genes, pancreatic_cancer_genes), ]))
expr_data$cell <- rownames(expr_data)
# 删除 cell 名称中的第二个 "_" 及其后缀
expr_data$cell <- sub("^([^_]+_[^_]+)_.*", "\\1", expr_data$cell)

# 提取treatment分组信息
treatment_info <- scRNAsub@meta.data$treatment

# 计算每个基因与WF-14W组之间的相关性并绘制相关性图
for (gene in c(ovarian_cancer_genes, gastric_cancer_genes, pancreatic_cancer_genes)) {
  
  # 提取WF-14W组的表达数据
  expr_WF_14W <- expr_data[treatment_info == "WF-14W", c(gene, "cell"), drop = FALSE]
  
  # 提取每个癌症类型对应的基因的其他组的表达数据，并保存组名
  if (gene == "PAX8") {
    group <- "Ovarian cancer"
    expr_other_group <- expr_data[treatment_info == group, c(gene, "cell"), drop = FALSE]
  } else if (gene == "VIL1") {
    group <- "Gastric cancer"
    expr_other_group <- expr_data[treatment_info == group, c(gene, "cell"), drop = FALSE]
  } else if (gene %in% c("TIMP1", "FUT3","CEACAM5")) {
    group <- "Pancreatic cancer"
    expr_other_group <- expr_data[treatment_info == group, c(gene, "cell"), drop = FALSE]
  }
  
  
  # 挑选前200个表达量最大的细胞
  expr_WF_14W_top200 <- expr_WF_14W %>%
    arrange(desc(!!sym(gene))) %>%
    head(200)
  
  expr_other_group_top200 <- expr_other_group %>%
    arrange(desc(!!sym(gene))) %>%
    head(200)
  
  # 合并WF-14W组和当前组的表达数据
  merged_data <- data.frame(
    expr_WF_14W = expr_WF_14W_top200[[gene]],
    expr_other_group = expr_other_group_top200[[gene]], 
    by = "cell"
  )
  
  # 计算相关性
  cor_value <- cor(merged_data$expr_WF_14W, merged_data$expr_other_group, method = "pearson")
  
  # 线性回归计算R²和y函数
  lm_model <- lm(expr_WF_14W ~ expr_other_group, data = merged_data)
  r2_value <- summary(lm_model)$r.squared
  coef_values <- coef(lm_model)
  intercept <- coef_values[1]
  slope <- coef_values[2]
  
  # 创建y函数公式
  y_function <- paste0("y = ", round(slope, 2), "x ", ifelse(intercept >= 0, "+ ", "- "), round(abs(intercept), 2))
  
  # 打印结果
  print(paste("Correlation for gene", gene, "between WF-14W and", group, ":", cor_value))
  print(paste("R² value for gene", gene, ":", r2_value))
  print(paste("y function for gene", gene, ":", y_function))
  
  # 绘制相关性散点图
  plot <- ggscatter(
    merged_data,
    x = "expr_WF_14W",
    y = "expr_other_group",
    color = I(gene_colors[gene]), 
    add = "reg.line",
    conf.int = TRUE,
    add.params = list(color = "black", fill = "lightgray"),
    cor.coef = FALSE,
    cor.coeff.args = list(size = 7, color = "blue"),
    cor.method = "pearson"
  ) +
    labs(
      title = paste0(gene, "\nR²: ", round(r2_value, 3), " | ", y_function),
      x = "Expression in WF-14W", 
      y = paste0("Expression in ", group)  # 根据组名动态设置
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold"),  
      axis.title.y = element_text(size = 20, face = "bold"),  
      axis.text.x = element_text(size = 18),  
      axis.text.y = element_text(size = 18)  
    )
  
  
  # 保存图表
  ggsave(file.path(output, paste0(gene, "_correlation_", gsub(" ", "_", group), ".pdf")), plot = plot, device = "pdf", width = 5, height = 5.5)
  ggsave(file.path(output, paste0(gene, "_correlation_", gsub(" ", "_", group), ".svg")), plot = plot, device = "svg", width = 5, height = 5.5)
}








####################### 单个基因相关性分析 (treatment 分组) ##############################

# 加载必要的R包
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 设置输出目录
output <- paste(outdir, 'Correlation_treatment', sep='/')
dir.create(output, showWarnings = FALSE)

# 读取数据
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scRNAsub <- readRDS(file_path)

# 为每个基因定义不同的颜色
colors <- c("#FF3366","#1E90FF","#FF6666", "#33CCCC")  # 颜色列表

# 定义每个癌症类型相关的基因列表
ovarian_cancer_genes <- "PAX8"        # 卵巢癌
gastric_cancer_genes <- "VIL1"        # 胃癌
pancreatic_cancer_genes <- c( "TIMP1", "FUT3")  # 胰腺癌

# 提取基因表达矩阵和 treatment 信息
expr_data <- as.data.frame(t(scRNAsub@assays$RNA@data[c(ovarian_cancer_genes, gastric_cancer_genes, pancreatic_cancer_genes), ]))
expr_data$cell <- rownames(expr_data)
# 删除 cell 名称中的第二个 "_" 及其后缀
expr_data$cell <- sub("^([^_]+_[^_]+)_.*", "\\1", expr_data$cell)

# 提取treatment分组信息
treatment_info <- scRNAsub@meta.data$treatment

# 计算每个基因与WF-14W组之间的相关性并绘制相关性图
for (gene in c(ovarian_cancer_genes, gastric_cancer_genes, pancreatic_cancer_genes)) {
  
  # 提取WF-14W组的表达数据
  expr_WF_14W <- expr_data[treatment_info == "WF-14W", c(gene, "cell"), drop = FALSE]
  
  # 提取每个癌症类型对应的基因的其他组的表达数据
  if (gene == "PAX8") {
    expr_other_group <- expr_data[treatment_info == "Ovarian cancer", c(gene, "cell"), drop = FALSE]
  } else if (gene == "VIL1") {
    expr_other_group <- expr_data[treatment_info == "Gastric cancer", c(gene, "cell"), drop = FALSE]
  } else if (gene %in% c( "TIMP1", "FUT3")) {
    expr_other_group <- expr_data[treatment_info == "Pancreatic cancer", c(gene, "cell"), drop = FALSE]
  }
  
  # 随机选择每组的100个细胞（确保每组至少有100个细胞）
  set.seed(123)  # 设置随机种子以确保结果可重复
  expr_WF_14W_sample <- expr_WF_14W %>% sample_n(100)
  expr_other_group_sample <- expr_other_group %>% sample_n(100)
  
  # 合并WF-14W组和当前组的表达数据
  merged_data <- data.frame(
    expr_WF_14W = expr_WF_14W_sample[[gene]],
    expr_other_group = expr_other_group_sample[[gene]]
  )
  
  # 计算该基因在WF-14W组与其他组之间的相关性
  cor_value <- cor(merged_data$expr_WF_14W, merged_data$expr_other_group, method = "pearson")
  print(paste("Correlation for gene", gene, "between WF-14W and other group:", cor_value))
  
  # 绘制相关性散点图
  plot <- ggscatter(
    merged_data,
    x = "expr_WF_14W",
    y = "expr_other_group",
    color = I(colors[i]),
    add = "reg.line",
    conf.int = TRUE,
    add.params = list(color = "black", fill = "lightgray"),
    cor.coef = TRUE,
    cor.coeff.args = list(size = 7, color = "blue"),
    cor.method = "pearson"
  ) +
    labs(title = paste0(gene),
         x = "Expression in WF-14W", y = paste0("Expression in ", group)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold"),  
      axis.title.y = element_text(size = 20, face = "bold"),  
      axis.text.x = element_text(size = 18),  
      axis.text.y = element_text(size = 18)  
    )
  
  # 保存图表为PDF和SVG格式
  ggsave(file.path(output, paste0(gene, "_correlation.pdf")), plot = plot, device = "pdf", width = 5, height = 4.5)
  ggsave(file.path(output, paste0(gene, "_correlation.svg")), plot = plot, device = "svg", width = 5, height = 4.5)
}










####################### 单个基因相关性分析 (treatment 分组) ##############################

# 加载必要的R包
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 设置输出目录
output <- paste(outdir, 'Correlation_treatment', sep='/')
dir.create(output, showWarnings = FALSE)

# 读取数据
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scRNAsub <- readRDS(file_path)

# 定义感兴趣基因列表
interest_genes <- c(
  "PAX8",  #"ESR1","PGR",   ## 卵巢
  "VIL1",    ## 胃癌
  "TIMP1","FUT3"       ## 胰腺癌
)

# 为每个基因定义不同的颜色
colors <- c("#FF3366","#1E90FF","#FF6666", "#33CCCC")  # 颜色列表

# 提取基因表达矩阵和 treatment 信息
expr_data <- as.data.frame(t(scRNAsub@assays$RNA@data[interest_genes, ]))
expr_data$cell <- rownames(expr_data)
# 删除 cell 名称中的第二个 "_" 及其后缀
expr_data$cell <- sub("^([^_]+_[^_]+)_.*", "\\1", expr_data$cell)

# 提取treatment分组信息
treatment_info <- scRNAsub@meta.data$treatment

# 获取WF-14W组之外的所有组名
unique_groups <- unique(treatment_info)
unique_groups <- unique_groups[unique_groups != "WF-14W"]  # 排除"WF-14W"组

# 计算每个基因与"WF-14W"组之间的相关性
for (i in seq_along(interest_genes)) {
  gene <- interest_genes[i]
  
  # 提取WF-14W组的表达数据
  expr_WF_14W <- expr_data[treatment_info == "WF-14W", gene, drop = FALSE]
  
  for (group in unique_groups) {
    
    # 提取当前组的表达数据
    expr_group <- expr_data[treatment_info == group, gene, drop = FALSE]
    
    # 删除表达量为0的细胞
    #expr_WF_14W <- expr_WF_14W[expr_WF_14W[[gene]] != 0, , drop = FALSE]
    #expr_group <- expr_group[expr_group[[gene]] != 0, , drop = FALSE]
    
    # 挑选前50个表达量最大的细胞
    expr_WF_14W_top50 <- expr_WF_14W %>%
      arrange(desc(!!sym(gene))) %>%
      head(50)
    
    expr_group_top50 <- expr_group %>%
      arrange(desc(!!sym(gene))) %>%
      head(50)
    
    # 合并"WF-14W"组和当前组的表达数据
    merged_data <- data.frame(
      expr_WF_14W = expr_WF_14W_top50[[gene]],
      expr_group = expr_group_top50[[gene]]
    )
    
    # 计算该基因在"WF-14W"组与当前组之间的相关性
    cor_value <- cor(merged_data$expr_WF_14W, merged_data$expr_group, method = "pearson")
    print(paste("Correlation for gene", gene, "between WF-14W and", group, ":", cor_value))
    
    # 绘制相关性散点图，并根据基因修改散点颜色
    plot <- ggscatter(
      merged_data,
      x = "expr_WF_14W",
      y = "expr_group",
      color = I(colors[i]),
      add = "reg.line",
      conf.int = TRUE,
      add.params = list(color = "black", fill = "lightgray"),
      cor.coef = TRUE,
      cor.coeff.args = list(size = 7, color = "blue"),
      cor.method = "pearson"
    ) +
      labs(title = paste0(gene),
           x = "Expression in WF-14W", y = paste0("Expression in ", group)) +
      theme_minimal(base_size = 14) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),  
        axis.title.y = element_text(size = 20, face = "bold"),  
        axis.text.x = element_text(size = 18),  
        axis.text.y = element_text(size = 18)  
      )
    
    # 保存图表为PDF和SVG格式
    ggsave(file.path(output, paste0(gene, "_correlation_", group, ".pdf")), plot = plot, device = "pdf", width = 5, height = 5)
    ggsave(file.path(output, paste0(gene, "_correlation_", group, ".svg")), plot = plot, device = "svg", width = 5, height = 5)
  }
}










