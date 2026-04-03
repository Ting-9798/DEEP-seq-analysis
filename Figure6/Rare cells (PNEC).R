
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
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
         "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/PNEC/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/PNEC/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

data <- readRDS("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/celltype.rds")

# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
ScRNA <- subset(data, idents = c("PNEC"))
table(ScRNA@meta.data$celltype)


# 挑选上皮细胞并且属于 "Tumor"、"Res" 或 "Sen" 组的细胞
#Cells.sub <- subset(data@meta.data, celltype == "PNEC")
#summary(Cells.sub$celltype)
#ScRNA <- subset(data, cells=row.names(Cells.sub))
#View(ScRNA@meta.data)

#### 6.归一化与PCA降维 ####
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

save(ScRNA, file = "ScRNA（批次矫正后分群前）.RData")



#### 7.细胞分群与注释 ####

col <- c('#FF6666','#E5D2DD','#6A4C93',"#BC8F8F",'#FFCC99','#FF9999',"#FFCCCC",'#4F6272','#58A4C3',
         "#66CCCC",'#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#CCFFCC","#00CC66","#99FFFF", 
         "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#FF3300","#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

load("ScRNA（批次矫正后分群前）.RData")
#细胞分群
ScRNA <- ScRNA %>% 
  RunUMAP(dims = 1:20) %>% 
  RunTSNE(dims = 1:20) %>%
  FindNeighbors(dims = 1:20)

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
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 8, height = 4)
DimPlot(ScRNA, reduction = "umap", pt.size=2,label = TRUE, repel = TRUE, split.by = "treatment", label.size = 5,cols = col)+
  theme(
    strip.text = element_text(size = 22, face = "bold"),  # 增大子图标题字体
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
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 5.5, height = 5)
DimPlot(ScRNA, reduction = "umap", pt.size=2,label = TRUE, repel = TRUE, label.size = 6,cols = col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),  # 增大X轴标题大小
        axis.title.y = element_text(size = 18, face = "bold"),  # 增大Y轴标题大小
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # 增大标题大小
        legend.title = element_text(size = 22),  # 增大图例标题大小
        legend.text = element_text(size = 22))
dev.off()

pdf(paste(OUTPUT, "cluster_umap1.pdf"), width = 6, height = 4)
DimPlot(ScRNA, reduction = "umap", pt.size=2,label = FALSE, repel = TRUE, cols = col, group.by ="treatment")+ ggtitle(NULL)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),  # 增大X轴标题大小
        axis.title.y = element_text(size = 18, face = "bold"),  # 增大Y轴标题大小
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  # 增大标题大小
        legend.title = element_text(size = 22),  # 增大图例标题大小
        legend.text = element_text(size = 22))
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




##########查看"tdTomato"和"Epcam"的表达比例###########
#####合并绘制######
genes <- c("tdTomato", "Epcam","Uchl1","Resp18")
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


# 需要绘图的分组和基因
treatment_groups <- c("0d", "3d","7d", "14d")
genes <- c("tdTomato", "Epcam","Uchl1","Resp18")

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
DoHeatmap(ScRNA, features = top5$gene,size = 5)+
  scale_fill_gradientn(colors = c("#437eb8", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # 增大X轴标题文本大小
        axis.title.y = element_text(size = 18),  # 增大Y轴标题文本大小
        legend.text = element_text(size = 14),  # 修改图例文本大小
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 30))  # 修改分组标签的颜色和大小
dev.off()

pdf(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 70,height = 10)
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis()：倾斜X轴文本
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 18),  # 增大X轴标题文本大小
        axis.title.y = element_text(size = 18))  # 增大Y轴标题文本大小
dev.off()
dpi=300
png(paste0(output,"/DotPlot_all_cluster_tsne_",max(dim.use),"PC.png"),w=70*dpi,h=10*dpi,units = "px",res = dpi,type='cairo')
DotPlot(ScRNA, features = unique(top5$gene))+RotatedAxis()+  #RotatedAxis()：倾斜X轴文本
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366"))+
  theme(axis.text = element_text(size = 16), 
        axis.title.x = element_text(size = 20),  # 增大X轴标题文本大小
        axis.title.y = element_text(size = 20))  # 增大Y轴标题文本大小
dev.off()






#########marker基因在细胞中的表达趋势#######
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 挑选差异细胞展示
output <- paste(outdir,'cluster(Club+AT2+Cilliated)', sep='/')
dir.create(output)

#file_path <- file.path(outdir, "celltype.rds")
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)


# 设置T细胞激活相关基因
cellmarker <- c(
  "Ccn1","Ccn2","Resp18","Cbr2","Hes1","Krt8" ,"Foxj1","Scgb3a2","Lamp3","Ager"    #PNEC (肺神经内分泌细胞)
  # "Cbr2","Notch2","Hes1","Myb","Calca","Mki67",
  
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
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", pt.size=3,
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
pdf(paste0(output, "/spacial_FeaturePlot_umap.pdf"), width = 15, height = 12)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

svg(paste0(output, "/spacial_FeaturePlot_umap.svg"), width = 15, height = 12)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()

pdf(paste0(output, "/spacial_VlnPlot_umap.pdf"), width = 12, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()
svg(paste0(output, "/spacial_VlnPlot_umap.svg"), width = 12, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()





# 设置颜色（如果样本较多请酌情调整颜色列表）
col_sample <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
                '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
                "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
                "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")


ScRNA$`treatment` <- factor(ScRNA$`treatment`, levels = c("0d", "3d","7d", "14d"))
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
pdf(paste0(output, "/cellmarker_VlnPlot_bySample.pdf"), width = 24, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 6))
dev.off()

svg(paste0(output, "/cellmarker_VlnPlot_bySample.svg"), width = 24, height = 6)
print(cowplot::plot_grid(plotlist = vln_plots_sample, ncol = 6))
dev.off()




############绘制基因表达热图###########
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)

# 设置分组顺序
ScRNA$treatment <- factor(ScRNA$treatment, levels = c("0d", "3d", "7d", "14d"))

# 提取表达矩阵并计算每组 treatment 的平均表达
avg_expr <- AverageExpression(ScRNA, features = cellmarker, group.by = "treatment", assays = "RNA")$RNA

# 转换为矩阵
avg_expr_mat <- as.matrix(avg_expr)

# 确保行名和 gene 对应
rownames(avg_expr_mat) <- rownames(avg_expr)

# 按指定基因顺序排序
gene_order <- cellmarker
avg_expr_mat <- avg_expr_mat[gene_order, , drop = FALSE]

# Z-score 转换
mat_scaled <- t(scale(t(avg_expr_mat)))

# 设置分组信息
group <- colnames(avg_expr_mat)  # 即 treatment 分组
group_anno <- HeatmapAnnotation(
  Treatment = factor(group, levels = c("0d", "3d", "7d", "14d")),
  col = list(Treatment = c("0d" = "#99CCFF", "3d" = "#66CCCC", "7d" = "#FF9933", "14d" = "#CC0066"))
)

# 保存热图（PDF格式）
pdf(file.path(output, "cellmarker_Heatmap_byTreatment.pdf"), width = 6, height = 3.5)
Heatmap(mat_scaled,
        name = "Z-score",
        top_annotation = group_anno,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 12),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = colorRampPalette(c('#3399CC', "white", "#FF3366"))(100),
        heatmap_legend_param = list(
          title = "Z-score",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 10))
) %>% draw()
dev.off()





#########marker基因在细胞中的表达趋势#######
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")



# 加载必要的包
library(ggpubr)
library(ggplot2)
library(dplyr)

#setwd("D:/R/GS/YY/20241218-fei-A/out(rna)/")
#outdir <- "D:/R/GS/YY/20241218-fei-A/out(rna)/"

# 挑选差异细胞展示
output <- paste(outdir,'cluster(Club+AT2+Cilliated)', sep='/')
dir.create(output)

#file_path <- file.path(outdir, "celltype.rds")
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)


# 设置T细胞激活相关基因
cellmarker <- c(
  "Ccn1","Ccn2"    #PNEC (肺神经内分泌细胞)
  # "Cbr2","Notch2","Hes1","Myb","Calca","Mki67",
  
)



# 创建存储图形的列表
boxviolin_plots_sample <- list()

for (gene in cellmarker) {
  # 获取该基因的表达数据
  gene_exp <- ScRNA@assays$RNA@data[gene, ]
  
  # 计算1%和99%分位数
  q01 <- quantile(gene_exp, 0.01, na.rm = TRUE)
  q99 <- quantile(gene_exp, 0.99, na.rm = TRUE)
  
  # 去除低表达的5%细胞和高表达的99%细胞
  filtered_cells <- which(gene_exp > q01)
  
  # 创建一个新的Seurat对象，包含过滤后的细胞
  ScRNA_filtered <- subset(ScRNA, cells = names(filtered_cells))
  ScRNA_filtered <- subset(ScRNA)
  
  # 提取表达数据和元数据用于统计检验
  exp_data <- as.numeric(ScRNA_filtered@assays$RNA@data[gene, ])
  
  treatment_groups <- as.character(ScRNA_filtered@meta.data$treatment)
  treatment_groups <- factor(ScRNA_filtered$treatment, levels = c("0d", "3d", "7d", "14d"))
  
  # 创建数据框用于绘图和统计检验
  plot_data <- data.frame(
    expression = exp_data,
    treatment = factor(treatment_groups)
  )
  
  # 移除NA值
  plot_data <- plot_data[!is.na(plot_data$expression), ]
  
  # 检查是否有足够的组和数据进行统计检验
  unique_treatments <- levels(plot_data$treatment)
  
  # 计算y轴范围，用于调整P值标签高度
  y_max <- max(plot_data$expression, na.rm = TRUE)
  y_range <- y_max - min(plot_data$expression, na.rm = TRUE)
  
  # 动态调整标签高度参数
  step_increase_factor <- 0.15  # 基础增加因子
  bracket_adjustment <- 0.08 * y_range  # 括号调整量
  label_y_adjustment <- 0.1 * y_range  # 标签y位置调整量
  
  # 创建箱琴图（箱线图+小提琴图）
  p <- ggplot(plot_data, aes(x = treatment, y = expression, fill = treatment)) +
    # 小提琴图（设置透明度）
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    # 箱线图
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, 
                 fill = "white", color = "black") +
    # 设置颜色
    scale_fill_manual(values = col) +
    # 主题设置
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    labs(
      title = gene,
      y = "Expression Level"
    )
  
  # 添加统计检验（根据组数选择合适的方法）
  if (length(unique_treatments) >= 2 && nrow(plot_data) > 0) {
    # 检查每组是否有足够的数据点
    group_counts <- table(plot_data$treatment)
    valid_groups <- names(group_counts)[group_counts >= 3]
    
    if (length(valid_groups) >= 2) {
      # 过滤数据，只保留有效组
      plot_data_filtered <- plot_data %>% 
        filter(treatment %in% valid_groups)
      
      # 计算过滤后的y轴范围
      y_max_filtered <- max(plot_data_filtered$expression, na.rm = TRUE)
      y_min_filtered <- min(plot_data_filtered$expression, na.rm = TRUE)
      y_range_filtered <- y_max_filtered - y_min_filtered
      
      # 根据组数选择统计方法
      if (length(valid_groups) == 2) {
        # 两组比较：使用Wilcoxon秩和检验
        test_result <- wilcox.test(expression ~ treatment, data = plot_data_filtered)
        p_val <- test_result$p.value
        
        # 计算显著性标记的y位置
        # 第一个比较的标签位置
        label_y_base <- y_max_filtered + 0.1 * y_range_filtered
        
        # 添加显著性标记
        p <- p + 
          stat_compare_means(
            method = "wilcox.test",
            comparisons = list(valid_groups),
            label = "p.signif",
            size = 6,
            label.y = label_y_base,  # 设置标签的y位置
            tip.length = 0.01,
            bracket.size = 0.6,
            vjust = 0.5  # 垂直对齐
          )
        
        # 添加精确p值作为副标题（可选）
        p_label <- ifelse(p_val < 0.001, "p < 0.001", 
                          ifelse(p_val < 0.01, "p < 0.01",
                                 sprintf("p = %.3f", p_val)))
        
      } else {
        # 多组比较：使用Kruskal-Wallis检验
        kruskal_result <- kruskal.test(expression ~ treatment, data = plot_data_filtered)
        p_val_kruskal <- kruskal_result$p.value
        
        # 添加整体检验结果作为副标题
        p_label_kruskal <- ifelse(p_val_kruskal < 0.001, "p < 0.001", 
                                  ifelse(p_val_kruskal < 0.01, "p < 0.01",
                                         sprintf("p = %.3f", p_val_kruskal)))
        p <- p + labs(subtitle = paste("Kruskal-Wallis test:", p_label_kruskal))
        
        # 如果整体显著，进行两两比较
        if (p_val_kruskal < 0.05) {
          # 生成所有两两组合
          comparisons_list <- combn(valid_groups, 2, simplify = FALSE)
          
          # 计算比较数量
          n_comparisons <- length(comparisons_list)
          
          # 动态调整step.increase基于比较数量
          step_increase <- min(0.15, 0.6 / n_comparisons)
          
          # 计算第一个比较的基础y位置
          base_y <- y_max_filtered + 0.15 * y_range_filtered
          
          # 添加两两比较显著性标记
          p <- p + 
            stat_compare_means(
              method = "wilcox.test",
              comparisons = comparisons_list,
              label = "p.signif",
              size = 6,
              tip.length = 0.01,
              bracket.size = 0.5,
              step.increase = step_increase,  # 动态调整步进
              label.y.npc = "top",  # 从顶部开始
              label.y = base_y,     # 设置基础y位置
              hide.ns = TRUE
            )
        }
      }
    }
  }
  
  # 可选：根据数据范围调整y轴上限，为显著性标记留出空间
  # 计算当前图形数据范围
  current_ymax <- layer_scales(p)$y$range$range[2]
  new_ymax <- current_ymax * 1.2  # 增加20%空间给标签
  
  # 设置y轴限制
  p <- p + coord_cartesian(ylim = c(min(plot_data$expression, na.rm = TRUE), new_ymax))
  
  # 添加到列表
  boxviolin_plots_sample[[gene]] <- p
}

# 保存图形
pdf(paste0(output, "/cellmarker_BoxViolinPlot_bySample_filtered_with_stats.pdf"), 
    width = 16, height = 4)
print(cowplot::plot_grid(plotlist = boxviolin_plots_sample, ncol = 5))
dev.off()

svg(paste0(output, "/cellmarkin_BoxViolinPlot_bySample_filtered_with_stats.svg"), 
    width = 16, height = 4)
print(cowplot::plot_grid(plotlist = boxviolin_plots_sample, ncol = 5))
dev.off()





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
setwd("D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/PNEC/")
outdir <- "D:/R/GS/WH/20250903-稀有细胞/out(0+3+7+14)tdt+(1)/去杂/PNEC/"

# 创建输出目录
output <- file.path(outdir, "差异分析")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# 读取数据
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scRNAsub <- readRDS(file_path)

logFCfilter <- 0.25
adjPvalFilter <- 0.05

genes <- c("Gpnmb","Spp1","Ctsd","Nfkb1","Ifngr1","Camk4","Zeb1","Rora","Cd14","Il1b",
           "Il12rb2","Cxcl2","Fth1","Icos","Stat4","Cd63","Thbs1","Tyrobp")

timepoints <- c("7d","14d")

all_markers <- list()

for (tp in timepoints) {
  
  
  comp_name <- paste0(tp, " vs 3d")
  
  markers <- FindMarkers(object = scRNAsub,
                         ident.1 = tp,
                         ident.2 = "3d",
                         group.by = "treatment",
                         logfc.threshold = 0,
                         min.pct = 0.25,
                         test.use = "wilcox")
  
  markers$gene <- rownames(markers)
  markers <- markers %>%
    mutate(Significance = ifelse(p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter, 
                                 ifelse(avg_log2FC > 0, "Up", "Down"), "Normal")) %>%
    mutate(Group = comp_name) 
  
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
    filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(p_val_adj) 
  top_genes_downregulated <- downregulated_genes_df %>%
    filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
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
    scale_x_continuous(limits = c(-1.5, 2.5), breaks = seq(-1.5, 2.5, by = 1)) +
    scale_y_continuous(limits = c(0, max(-log10(markers$p_val_adj), na.rm = TRUE) + 1))+
    #scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 5)) +
    theme(plot.title = element_text(size = 22, face = "bold", hjust = 0), 
          legend.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 20, hjust = 0.5),
          axis.text = element_text(size = 18))
  
  
  ggsave(file.path(output, paste0(comp_name, "_volcano_plot.svg")), p, width = 8, height = 7)
  ggsave(file.path(output, paste0(comp_name, "_volcano_plot.pdf")), p, width = 8, height = 7)
  
  
  all_markers[[tp]] <- markers
}


# 合并三个组
merged_markers <- bind_rows(all_markers)

# 保存多组输入表格
write.csv(merged_markers, file = file.path(output, "merged_markers_for_mutiVolcano.csv"), row.names = FALSE)


#########绘制多组差异火山图###############
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)
library(ggsci)
library(patchwork)
library(tidyr)

# 加载数据
merged_markers <- read.csv(file.path(output, "merged_markers_for_mutiVolcano.csv"))

# 数据预处理
data <- merged_markers %>%
  rename(logFC = avg_log2FC, adj.P.Val = p_val_adj) %>%
  mutate(change = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.25,
                         ifelse(logFC > 0, "Up", "Down"),
                         "No change")) %>%
  filter(change != "No change") %>%
  mutate(group = Group) %>%
  mutate(label = change)

# 设置分组顺序
group_order <- c("7d vs 3d", "14d vs 3d")
data <- data %>%
  mutate(group = factor(group, levels = group_order))

# 分别提取每组Up和Down中的Top20
TopGene_up <- data %>%
  filter(change == "Up") %>%
  group_by(group) %>%
  distinct(gene, .keep_all = TRUE) %>%
  top_n(20, wt = logFC)

TopGene_down <- data %>%
  filter(change == "Down") %>%
  group_by(group) %>%
  distinct(gene, .keep_all = TRUE) %>%
  top_n(20, wt = -logFC)

TopGene <- bind_rows(TopGene_up, TopGene_down) %>%
  ungroup()

# 背景柱状图数据
dbar <- data %>%
  group_by(group) %>%
  summarise(logFC_min = min(logFC),
            logFC_max = max(logFC))

# 绘图
p <- ggplot() +
  # 背景柱状图
  geom_col(data = dbar, aes(x = group, y = logFC_min), fill = "#dcdcdc", alpha = 0.6, width = 0.7) +
  geom_col(data = dbar, aes(x = group, y = logFC_max), fill = "#dcdcdc", alpha = 0.6, width = 0.7) +
  
  # 所有点图
  #geom_jitter(data = TopGene, aes(x = group, y = logFC, color = label), size = 0.85, width = 0.3) +
  
  # Top基因点图
  geom_jitter(data = TopGene, aes(x = group, y = logFC, color = label), size = 1.5, width = 0.35) +
  
  # 中间tile
  geom_tile(data = TopGene, aes(x = group, y = 0, fill = group), height = 0.6, color = "black", alpha = 0.6, show.legend = FALSE) +
  
  # Top基因标签，标签位置向 logFC 方向偏移 0.1，防止重叠
  geom_text_repel(data = TopGene,
                  aes(x = group, y = logFC - 0.1 * sign(logFC), label = gene),
                  size = 5, color = 'black',
                  force = 1.2,
                  max.overlaps = 50,
                  arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
  
  # 中心分组标签
  geom_text(data = TopGene, aes(x = group, y = 0, label = group), size = 5, color = "white") +
  
  ggsci::scale_fill_npg() +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5"))+
  labs(x="",y = "log2 Fold Change", color = "") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 18)
  )

# 保存图像
ggsave(filename = file.path(output, "multi_volcano_plot.pdf"), plot = p, width = 6, height = 5, bg = "white")
ggsave(filename = file.path(output, "multi_volcano_plot.svg"), plot = p, width = 6, height = 5, bg = "white")




# ======== 功能分析：批量进行 GSEA、GO、KEGG 分析 ========
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)

# 遍历每一组差异分析结果
for (tp in timepoints) {
  comp_name <- paste0(tp, " vs 3d")
  marker_file <- file.path(output, paste0("ScRNA.sig.markers_", comp_name, ".rds"))
  markers <- readRDS(marker_file)
  
  message(">>> 正在处理: ", comp_name)
  
  
  # ---- GO & KEGG 富集分析 ----
  gene_up <- markers$gene[markers$Significance == "Up"]
  gene_down <- markers$gene[markers$Significance == "Down"]
  
  gene_up_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                               keys = gene_up, 
                                                               columns = 'ENTREZID', 
                                                               keytype = 'SYMBOL')[,2]))
  gene_down_entrez <- as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db, 
                                                                 keys = gene_down, 
                                                                 columns = 'ENTREZID', 
                                                                 keytype = 'SYMBOL')[,2]))
  
  # GO 富集分析
  go_up <- enrichGO(gene = gene_up_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.1)
  go_down <- enrichGO(gene = gene_down_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.1)
  
  # 删除 GO 分析中 Description 中的物种后缀
  go_up@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", go_up@result$Description)
  go_down@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", go_down@result$Description)
  
  # 将geneID从ENTREZID转为SYMBOL
  go_up@result$geneID <- sapply(strsplit(go_up@result$geneID, "/"), function(ids) {
    symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
    paste(symbols, collapse = "/")
  })
  
  go_down@result$geneID <- sapply(strsplit(go_down@result$geneID, "/"), function(ids) {
    symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
    paste(symbols, collapse = "/")
  })
  
  write.csv(as.data.frame(go_up), file = file.path(output, paste0(comp_name, "_GO_UP.csv")))
  write.csv(as.data.frame(go_down), file = file.path(output, paste0(comp_name, "_GO_DOWN.csv")))
  
  if (nrow(go_up) > 0) {
    gop_up <- dotplot(go_up) + ggtitle(paste(comp_name, "GO UP"))
    ggsave(file.path(output, paste0(comp_name, "_GO_UP.pdf")), gop_up, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_GO_UP.svg")), gop_up, width = 6, height = 6)
  }
  if (nrow(go_down) > 0) {
    gop_down <- dotplot(go_down) + ggtitle(paste(comp_name, "GO DOWN"))
    ggsave(file.path(output, paste0(comp_name, "_GO_DOWN.pdf")), gop_down, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_GO_DOWN.svg")), gop_down, width = 6, height = 6)
  }
  
  # KEGG 富集分析
  kegg_up <- enrichKEGG(gene = gene_up_entrez, organism = "mmu", pvalueCutoff = 0.1)
  kegg_down <- enrichKEGG(gene = gene_down_entrez, organism = "mmu", pvalueCutoff = 0.1)
  
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
  
  write.csv(as.data.frame(kegg_up), file = file.path(output, paste0(comp_name, "_KEGG_UP.csv")))
  write.csv(as.data.frame(kegg_down), file = file.path(output, paste0(comp_name, "_KEGG_DOWN.csv")))
  
  if (nrow(kegg_up) > 0) {
    kp_up <- dotplot(kegg_up) + ggtitle(paste(comp_name, "KEGG UP"))
    ggsave(file.path(output, paste0(comp_name, "_KEGG_UP.pdf")), kp_up, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_KEGG_UP.svg")), kp_up, width = 6, height = 6)
  }
  if (nrow(kegg_down) > 0) {
    kp_down <- dotplot(kegg_down) + ggtitle(paste(comp_name, "KEGG DOWN"))
    ggsave(file.path(output, paste0(comp_name, "_KEGG_DOWN.pdf")), kp_down, width = 6, height = 6)
    ggsave(file.path(output, paste0(comp_name, "_KEGG_DOWN.svg")), kp_down, width = 6, height = 6)
  }
}






#########绘制感兴趣的通路###############
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(scales)
library(ggsci)
library(patchwork)
library(tidyr)

# 加载数据
go_up <- read.csv(file.path(output, "14d vs 3d_GO_UP.csv"))


# 准备可视化数据
# 提取GO分析结果
go_up_dt <- as.data.frame(go_up)
#go_down_dt <- as.data.frame(go_down)


# 提取感兴趣的GO通路
interested_go <- c(
  "tight junction assembly","tight junction organization","bicellular tight junction assembly","transepithelial transport",
  "epithelial fluid transport","cell-cell junction assembly","lung alveolus development","epithelial cell proliferation",
  "collagen biosynthetic process","collagen metabolic process","regulation of tissue remodeling",
  "regulation of bone remodeling","prostaglandin biosynthetic process","prostaglandin metabolic process",
  "icosanoid biosynthetic process","icosanoid metabolic process","endothelial cell proliferation",
  "vascular transport","chemokine production",
  "interleukin-6 production","interleukin-8 production","tumor necrosis factor-mediated signaling pathway",
  "response to interferon-alpha","response to type II interferon",
  "JAK-STAT cascade","NF-kappaB signaling","ERK1 and ERK2 cascade","p38MAPK cascade",

  "DNA damage response, signal transduction by p53 class mediator","signal transduction by p53 class mediator",
  "signal transduction in response to DNA damage"
)

# 从GO富集结果中筛选感兴趣的通路
go_up_dt <- go_up_dt[go_up_dt$Description %in% interested_go, ]


# 设置颜色分类
classification_colors <- c('#437eb8','#FF6666','#FFCC99','#FF9999', '#80c5d8',"#9999FF",
                           "#FFCCCC","#99CCFF","#FF3366","#CCCCFF","#CC0066","#FFFFCC",
                           "#66CCCC","#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                           "#6699CC","#CC99CC","#FF6699","#FF0000","#6666CC","#FF9966",
                           "#669999","#CC99FF","#FFCCFF",
                           '#437eb8','#FF6666','#FFCC99','#FF9999', '#80c5d8',"#9999FF",
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
  
  #dt <- head(dt,20)  # 选取前20个通路
  dt$Description <- factor(wrap_text(dt$Description), levels = wrap_text(dt$Description))
  
  # 左图：富集p值
  p1 <- ggplot(dt, aes(x = Description, y = log10(p.adjust), fill = Description)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = classification_colors) +
    coord_flip() +
    ylab('-log10(P-value)') +
    xlab('') +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
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
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(size = 14, face = "bold"),
          plot.margin = margin(10, 10, 10, 10),
          panel.border = element_rect(color = "black", fill = NA, size = 1))  # 添加XY轴边框
  
  # 组合两个图
  p_combined <- p1 + p2 + plot_layout(widths = c(2, 1.5))
  return(p_combined)
}

# 绘制GO和KEGG富集分析条形图
go_up_plot <- plot_GO_bar(go_up_dt, "Upregulated Genes GO Enrichment")
#go_down_plot <- plot_GO_bar(go_down_dt, "Downregulated Genes GO Enrichment")

#kegg_up_plot <- plot_GO_bar(kegg_up_dt, "Upregulated Genes KEGG Enrichment")
#kegg_down_plot <- plot_GO_bar(kegg_down_dt, "Downregulated Genes KEGG Enrichment")

# 保存图像
ggsave(file.path(output, 'go_enrich_up_bar.pdf'), plot = go_up_plot, width = 8, height = 6)
#ggsave(file.path(output, 'go_enrich_down_bar.pdf'), plot = go_down_plot, width = 8, height = 6)
#ggsave(file.path(output, 'kegg_enrich_up_bar.pdf'), plot = kegg_up_plot, width = 8, height = 6)
#ggsave(file.path(output, 'kegg_enrich_down_bar.pdf'), plot = kegg_down_plot, width = 8, height = 6)

# 组合并保存
#combined_GO_plot <- go_up_plot + go_down_plot + plot_layout(guides = 'collect')
#combined_KEGG_plot <- kegg_up_plot + kegg_down_plot + plot_layout(guides = 'collect')

#ggsave(file.path(output, 'combined_GO_bar.pdf'), plot = combined_GO_plot, width = 13, height = 10)
#ggsave(file.path(output, 'combined_KEGG_bar.pdf'), plot = combined_KEGG_plot, width = 13, height = 10)







############ 批量绘制通路基因表达热图 ############
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(dplyr)

output <- file.path(outdir, "marker")
dir.create(output, showWarnings = FALSE)

# 加载数据
ScRNA <- readRDS("ScRNA（分群后）.rds")

# 设置分组顺序
ScRNA$treatment <- factor(ScRNA$treatment, levels = c("3d", "7d", "14d"))

# 设置多个通路的基因列表
pathway_genes <- list(
  ## p53 signaling pathway
  p53 = c("Trp53","Mdm2","Ccnd1","Ccnd2","Ccnd3","Ccne1","Ccne2",
          "Ccng1","Ccng2","Cdk1","Cdk2","Cdk4","Cdk6","Cdkn1a","Gadd45a","Gadd45g",
          "Ei24","Ppm1d","Bax","Bbc3","Bid","Apaf1","Ddit4","Pig3",
          "Casp3","Casp8"),
  
  ## Wnt signaling pathway
  Wnt = c("Axin1","Axin2","Ctnnb1","Gsk3β","Wnt3","Wnt5a","Wnt7b",
          "Lrp5","Lrp6","WISP1","Tnc","C-MYC","Cyclin D","VEGF"),
  
  ## Hippo signaling pathway
  Hippo = c("Mst1","Mst2","Stk3","Stk4","Lats1","Lats2","Yap1","Taz","Ajuba",
            "Wwtr1","Ctgf","Serpine1","Tead1","Tead2","Tead3","Tead4","Mob1","Sav1"),
  
  ## JAK-STAT signaling pathway
  JAK_STAT = c("Jak2","Jak3","Tyk2","Stat1","Stat3","Stat5","Stat6","Il11","Socs1","Socs3"),
  
  ## PI3K-Akt signaling pathway
  PI3K_Akt = c("Pik3ca","Pik3cb","Pik3r1","Akt1","Akt2","Akt3","Mtor","Rictor","Rheb",
               "Tsc1","Tsc2","Pten","Sgk1","Fgf2","Igf1","Igf1r","Cd274","Src","Dab2",
               "Hmox1","Nqo1","Nfe2l2"),
  
  spacial = c("Piezo1","Piezo2","Osr2","Egr3"),
  
  Hedgehog = c("Shh","Ptch1","Ptch2","Smo","Sufu","Kif7","Gli3","Gas1","Boc","Cdon","Hhip"),
  
  Tgfb = c("Tgfb1","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2","Tgfbr3","Acvr1","Acvr2a","Acvr2b","Bmpr1a","Bmpr1b","Bmpr2","Smad2","Smad3","Smad4","Smad1","Smad5","Smad9","Smad6","Smad7","Map3k7"),
  
  EGF = c("Tgfa","Hbegf","Areg","Btc","Nrg1","Nrg2","Nrg3","Nrg4","Egfr","Erbb2","Erbb3","Erbb4","Shc1","Grb2","Sos1","Sos2","Gab1","Kras","Hras","Nras","Raf1","Braf","Map2k1","Map2k2","Mapk1","Mapk3","Pik3ca","Pik3cb","Pik3r1","Akt1","Akt2","Akt3","Mtor","Plcg1","Jak2","Stat3","Stat5","Pten","Cbl"),
  
  FGF = c("Fgf1","Fgf2","Fgf7","Fgf9","Fgf18","Fgf21","Fgfr1","Fgfr2","Fgfr3","Frs2","Grb2","Sos1","Gab1","Pik3ca","Pik3cb","Pik3r1","Akt1","Akt2","Akt3","Plcg1","Mapk1","Mapk3","Klb","Kl","Sdc1","Sdc2","Sdc3","Sdc4"),
  
  Rb = c("Rb1","Rbl1","Rbl2","E2f1","E2f2","E2f3","E2f4","E2f5","E2f6","E2f7","Ccnd1","Ccnd2","Ccnd3","Cdk4","Cdk6","Ccne1","Ccne2","Cdk2","Cdkn2b","Cdkn1a","Cdkn1b")
  
 
  
)

# 循环绘制每个通路的热图
for (pathway_name in names(pathway_genes)) {
  
  genes <- pathway_genes[[pathway_name]]
  
  # 提取表达矩阵并计算每组 treatment 的平均表达
  avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "treatment", assays = "RNA")$RNA
  
  avg_expr_mat <- as.matrix(avg_expr)
  rownames(avg_expr_mat) <- rownames(avg_expr)
  
  # Z-score 转换
  mat_scaled <- t(scale(t(avg_expr_mat)))
  
  # 分组信息
  group <- colnames(avg_expr_mat)  # treatment
  group_anno <- HeatmapAnnotation(
    Treatment = factor(group, levels = c("0d", "3d", "7d", "14d")),
    col = list(Treatment = c( "3d" = "#66CCCC", "7d" = "#FF9933", "14d" = "#CC0066"))   # "0d" = "#99CCFF",
  )
  
  # 热图颜色渐变
  heatmap_col <- colorRampPalette(c('#3399CC', "white", "#FF3366"))(100)
  
  # 保存热图
  pdf(file.path(output, paste0("Heatmap_", pathway_name, "_byTreatment.pdf")), width = 6, height = 5)
  Heatmap(mat_scaled,
          name = "Z-score",
          top_annotation = group_anno,
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 12),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          col = heatmap_col,
          heatmap_legend_param = list(
            title = "Z-score",
            title_gp = gpar(fontsize = 12),
            labels_gp = gpar(fontsize = 10)
          )) %>% draw()
  dev.off()
}




## 推断拟时序起点-CytoTRACE ##

#install
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") 
library(CytoTRACE2)
library(tidyverse)
library(Seurat)


col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC",'#E5D2DD','#4F6272',"#CC99CC",
         '#F9BB72', '#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',
         "#FF9933","#9999FF","#00CC66","#99FFFF","#FF3300",
         "#CCCCFF","#FF6699","#6699CC","#FFFFCC")


output <- paste(outdir,'推断拟时序起点CytoTRACE', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
data1 <- readRDS(file_path)
summary(data1$treatment)


data1@meta.data$CB <- rownames(data1@meta.data)
sample_CB <- data1@meta.data %>% 
  group_by(celltype) %>% 
  sample_frac(0.3)
sce3 <- subset(data1,CB %in% sample_CB$CB) 
sce3


#######输入seurat 对象###########
cytotrace2_result_sce <- cytotrace2(sce3, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'human',
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

col<- c(
  # UMAP
  "#31CDEE", "#D0F199", "#6666CC","#79BC98", "#3C8487", "#094867",'#E59CC4',
  "#FEDD81", "#FF9A84", "#9B6194", "#43457B","#1965B0","#CCFFCC","#CCCCFF",
  # 深蓝→绿→浅绿 梯度
  "#390A4A", "#443880", "#39558B", "#31668D", "#28818E",
  "#20958B", "#20A487", "#48C06E", "#6ECC5A", "#A6DC36",
  "#F5E24B",
  # Sum-seq 浅色
  "#82E1F6", "#E2F8C3", "#ADD8C0", "#89B5B2", "#6C92A0",
  "#32CBF1", "#FEDA84", "#FF9B84", "#966392", "#094869"
  
)

output <- paste(outdir,'monocle2(Club+AT2+AT1)', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
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

# 设置目标基因列表（Notch 通路相关）
target_genes <- c(
  "Foxj1","Tppp3","Tubb4b","Tubb1",   #Cilliated纤毛细胞,
  "Scgb1a1","Scgb3a2","Chad",     #club cell
  "Cd14","Cd74","H2-K1",    # Activated Club
  "H2-Ab1","Cst3",     # MHC-II+ Club
  'Ager', 'Hopx', 'Rtkn2', 'Aqp5',"Cav1 ","Spock2",   #AT1
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  "Lrg1","Lcn2","Retnla","Il33","Car8","Ank3","Cftr",         # Activated AT2
  "Birc5","Top2a",       # Proliferating AT2s
  "Cldn4","Sfn","Clu","Krt19","Krt8"    # PATs [pre-AT1 过渡状态]  Epithelial transitional states 
)

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
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # 对于离散型数据（State），使用颜色向量 col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # 使用 col 作为颜色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # 对于离散型数据（celltype 和 treatment），使用自定义颜色 custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
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





##根节点的确认
cds <- orderCells(cds,root_state=6)
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
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # 对于离散型数据（State），使用颜色向量 col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # 使用 col 作为颜色
      theme(legend.text = element_text(size = 14),  # 调整图例文本大小
            legend.title = element_text(size = 16),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # 对于离散型数据（celltype 和 treatment），使用自定义颜色 custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
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
         plot = plot_cell_traj, width = 4, height = 3)
  
  # 保存为 SVG 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 4, height = 3)
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
#topgene <- ordergene[1:50]
# 设置目标基因列表（Notch 通路相关）
target_genes <- c(
  "Foxj1","Tppp3","Tubb4b",   #Cilliated纤毛细胞,
  "Scgb1a1","Scgb3a2","Chad",     #club cell
  "Cd14","Cd74","H2-K1",    # Activated Club
  "H2-Ab1","Cst3",     # MHC-II+ Club
  'Ager', 'Hopx', 'Aqp5',"Cav1 ","Spock2",   #AT1
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  "Lrg1","Lcn2","Il33","Car8","Ank3","Cftr",         # Activated AT2
  "Birc5","Top2a",       # Proliferating AT2s
  "Cldn4","Sfn","Clu","Krt19","Krt8"    # PATs [pre-AT1 过渡状态]  Epithelial transitional states 
)

# 筛选目标基因中实际存在于表达矩阵中的部分
ordering_genes <- intersect(target_genes, rownames(expr_matrix))


plot_pseu_heatmap <- plot_pseudotime_heatmap(cds[ordering_genes,],num_clusters = 5,cores = 1,
                                             show_rownames = T,return_heatmap=T,hmcols = colorRampPalette(c("#1f77b4", "#ffffff", "#FF3366"))(100))
pdf(paste0(output,"/monocle_pheatmap.pdf"),width = 5, height = 4)
print(plot_pseu_heatmap)
dev.off()
ggsave(paste0(output,"/monocle_pheatmap.svg"),plot_pseu_heatmap,width = 5, height = 4)



## 关键驱动基因的表达变化图
#选择前4个基因
#keygenes <- ordergene[1:8] 
#cds_subset <- cds[keygenes,]

#print(ordergene)

# 感兴趣的基因（用户指定顺序）
interested_genes <- c("Cldn4", "Krt8", "Foxj1", "Scgb1a1", "Lamp3", "Ager")

# 确保基因存在于cds中，并按照interested_genes的顺序筛选
keygenes <- intersect(interested_genes, rownames(cds))
#keygenes <- interested_genes[interested_genes %in% keygenes]  # 保持用户指定顺序

# 检查选中的基因
if (length(keygenes) == 0) {
  stop("未找到感兴趣的基因，请检查基因名称是否正确。")
} else {
  message("找到以下感兴趣的基因（按顺序）：", paste(keygenes, collapse = ", "))
}

# 子集数据（保持基因顺序）
cds_subset <- cds[keygenes,]

# 自定义颜色
custom_colors <- c(
  "Tumor" = '#FF6666',
  "Normal" = '#E5D2DD',
  "BC-cancer" = '#FF6666',
  "BC-no_cancer" = '#E5D2DD'
)

# 遍历类型进行绘图
for (type in types) {
  if (type == "Pseudotime") {
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +
      theme(
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "top",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14)
      ) +
      facet_wrap(~ factor(gene_short_name, levels = keygenes), ncol = 2)
    
  } else if (type %in% c("State", "celltype", "treatment")) {
    plot_cell_pseu <- plot_genes_in_pseudotime(cds_subset, color_by = type) +
      xlab("Pseudotime") +
      scale_color_manual(values = col) +
      theme(
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "top",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14)
      ) +
      facet_wrap(~ factor(gene_short_name, levels = keygenes), ncol = 2) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  }
  
  # 输出 PDF 和 SVG 图像
  pdf(file = file.path(output, paste0("keygene_", type, ".pdf")), width = 6, height = 5)
  print(plot_cell_pseu)
  dev.off()
  
  ggsave(
    filename = file.path(output, paste0("keygene_", type, ".svg")),
    plot = plot_cell_pseu,
    width = 6,
    height = 5
  )
}







## 挖掘对细胞分化有重要影响的开关基因-geneswitches
## load libraries
#install.packages("fastglm")
#devtools::install_github("SGDDNB/GeneSwitches")
library(GeneSwitches)
library(SingleCellExperiment)


cardiac_monocle2 <- cds

monocle:::plot_cell_trajectory(cardiac_monocle2, color_by = "State")
plot_monocle_State(cardiac_monocle2)

## 2) 构建 log-normalized 表达矩阵：logexpdata
## 若尚未估计 size factors/离散度，先估计
if (is.null(sizeFactors(cardiac_monocle2)) || any(is.na(sizeFactors(cardiac_monocle2)))) {
  cardiac_monocle2 <- estimateSizeFactors(cardiac_monocle2)
}
if (is.null(dispersionTable(cardiac_monocle2))) {
  cardiac_monocle2 <- estimateDispersions(cardiac_monocle2)
}


## 取原始表达（exprs 是经过 Monocle 处理后的表达矩阵；依数据而定通常是 UMI 计数）
raw_expr <- exprs(cardiac_monocle2)            # dgCMatrix/矩阵，基因为行、细胞为列

## 用 size factor 做列归一，再 log1p
sf <- sizeFactors(cardiac_monocle2)             # 每个细胞的 size factor
norm_expr <- t(t(as.matrix(raw_expr)) / sf)     # 列归一
logexpdata <- log1p(norm_expr)                  # 等价于 log(norm_expr + 1)
## 如果你偏好 log2：
# logexpdata <- log2(norm_expr + 1)

## Input log-normalized gene expression, Monocle2 pseudo-time and dimensionality reduction
## Path1 containing cells in states 3,2,1
sce_p1 <- convert_monocle2(monocle2_obj = cardiac_monocle2, 
                           states = c(3,2,1), expdata = logexpdata)
## Path2 containing cells in states 3,2,5
sce_p2 <- convert_monocle2(monocle2_obj = cardiac_monocle2, 
                           states = c(3,2,5), expdata = logexpdata)


## ---- 仍保留你前面的代码 ----
## 二元化基因表达（Windows 串行，Linux/Mac 并行）
is_windows <- .Platform$OS.type == "windows"
safe_cores <- if (is_windows) 1 else max(1, parallel::detectCores() - 1)

message(sprintf("Running binarize_exp with ncores = %d (OS: %s)", 
                safe_cores, ifelse(is_windows, "Windows", "Unix-like")))

sce_p1 <- binarize_exp(sce_p1, ncores = safe_cores)
sce_p2 <- binarize_exp(sce_p2, ncores = safe_cores)


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

output <- paste(outdir,'monocle2', sep='/')
dir.create(output)

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
data1 <- readRDS(file_path)
summary(data1$treatment)

# 去除 0d
data <- subset(data1, subset = treatment != "0d")
summary(data$treatment)

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
deg <- subset(diff)  #qval < 0.01
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
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # 对于离散型数据（State），使用颜色向量 col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # 使用 col 作为颜色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # 对于离散型数据（celltype 和 treatment），使用自定义颜色 custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 3, show_backbone = TRUE) +
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
cds <- orderCells(cds,root_state=6)
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
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
      scale_color_gradient(low = "#1f77b4", high = "#FF3366") +  # 设置渐变色
      theme(legend.text = element_text(size = 16),  # 调整图例文本大小
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14))
  } else if (type %in%c("State","celltype")) {
    # 对于离散型数据（State），使用颜色向量 col
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
      scale_color_manual(values = col) +  # 使用 col 作为颜色
      theme(legend.text = element_text(size = 14),  # 调整图例文本大小
            legend.title = element_text(size = 16),
            axis.title.x = element_text(size = 16),  
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      guides(color = guide_legend(override.aes = list(size = 3, shape = 16)))
  } else if (type %in% c( "treatment")) {
    # 对于离散型数据（celltype 和 treatment），使用自定义颜色 custom_colors
    plot_cell_traj <- plot_cell_trajectory(cds, color_by = type, cell_size = 2, show_backbone = TRUE) +
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
         plot = plot_cell_traj, width = 6, height = 4)
  
  # 保存为 SVG 格式
  ggsave(filename = paste(output, paste0("monocle_", type, ".svg", sep = ""), sep = "/"), 
         plot = plot_cell_traj, width = 6, height = 4)
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
topgene <- c(
  "Foxj1","Tppp3","Tubb4b","Tubb1",   #Cilliated纤毛细胞,
  "Scgb1a1","Scgb3a2","Chad",     #club cell
  "Cd14","Cd74","H2-K1",    # Activated Club
  "H2-Ab1","Cst3",     # MHC-II+ Club
  'Ager', 'Hopx', 'Rtkn2', 'Aqp5',"Cav1 ","Spock2",   #AT1
  'Lamp3',  'Slc34a2', 'Lpcat1',"Sftpc", "Etv5",    #AT2
  "Lrg1","Lcn2","Retnla","Il33","Car8","Ank3","Cftr",         # Activated AT2
  "Birc5","Top2a",       # Proliferating AT2s
  "Cldn4","Sfn","Clu","Krt19","Krt8"    # PATs [pre-AT1 过渡状态]  Epithelial transitional states 
  
  )
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
sce <- readRDS("ScRNA（分群后）.rds")

# 2. 定义基因集
geneSet <- c(
  "Notch2","Notch1"
  #"Resp18","Uchl1","Syp","Pcsk1","Scg5","Chga","Chgb","Calca","Calcb","Sez6l2",
  #"Hes1","Uchl1","Ascl1","Krt8",
  #"Cbr2","Foxj1","Tppp3","Scgb3a2","Lamp3",'Lpcat1',"Sftpc", "Etv5","Ager"

)  


# 去重并提示哪些基因不在对象里（区分区分大小写）
genes_of_interest <- unique(geneSet)
genes_present     <- genes_of_interest[genes_of_interest %in% rownames(sce)]
genes_missing     <- setdiff(genes_of_interest, genes_present)

if (length(genes_missing) > 0) {
  message("以下基因未在对象中检测到（将跳过）：", paste(genes_missing, collapse = ", "))
}

stopifnot(length(genes_present) > 0)

#====================#
# 3) 确保降维可视化可用
#====================#
if (!"umap" %in% Reductions(sce)) {
  message("未检测到 UMAP，自动计算 PCA + UMAP ...")
  DefaultAssay(sce) <- DefaultAssay(sce) # 保持当前默认
  sce <- suppressMessages(RunPCA(sce, npcs = 50, verbose = FALSE))
  sce <- suppressMessages(RunUMAP(sce, dims = 1:30, verbose = FALSE))
}

#====================#
# 4) 实用函数：打分、分箱/分级、绘图
#====================#
# 连续 → 16 等分位箱（带极端并列值兜底）
make_bin16 <- function(x) {
  probs    <- seq(0, 1, length.out = 17)
  brks_raw <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
  brks     <- unique(brks_raw)
  if (length(brks) < 17) {
    brks <- pretty(range(x, na.rm = TRUE), n = 16)
    brks <- unique(brks)
  }
  if (length(brks) < 2) {
    brks <- range(x, na.rm = TRUE)
    brks[1] <- brks[1] - 1e-8
    brks[2] <- brks[2] + 1e-8
  }
  n_bins     <- length(brks) - 1
  bin_labels <- paste0("B", seq_len(n_bins))
  cut(x, breaks = brks, include.lowest = TRUE, labels = bin_labels, right = TRUE, ordered_result = TRUE)
}

# 连续 → 20/60/20 分级
make_hml <- function(x, prefix = "Met") {
  q20 <- as.numeric(quantile(x, 0.20, na.rm = TRUE, type = 7))
  q80 <- as.numeric(quantile(x, 0.80, na.rm = TRUE, type = 7))
  cls <- dplyr::case_when(
    x >= q80 ~ paste0(prefix, "_high"),
    x <= q20 ~ paste0(prefix, "_low"),
    TRUE     ~ paste0(prefix, "_medium")
  )
  list(class = factor(cls, levels = paste0(prefix, c("_high","_medium","_low"))),
       q20 = q20, q80 = q80)
}

# 选择降维
pick_reduc <- function(obj) if ("umap" %in% Reductions(obj)) "umap" else DefaultDimReduc(obj)

# 可选 split.by：若存在 treatment 列则自动分组
split_by <- if ("treatment" %in% colnames(sce@meta.data)) "treatment" else NULL

#====================#
# 5) 主循环：逐基因 UCell + 分箱/分级 + 绘图
#====================#
thresholds_tbl <- list()

for (g in genes_present) {
  message("Processing gene: ", g)
  g_safe <- g  # 列名里保留原始基因名
  
  # 5.1 UCell 单基因打分
  # 每次传入一个 named list，name = 基因名，value = 基因向量（长度=1）
  sce <- AddModuleScore_UCell(
    sce,
    features = setNames(list(g), g),
    name = "_UCell"  # 生成列名：<list_name>_UCell => 即 <基因>_UCell
  )
  
  score_col <- paste0(g_safe, "_UCell")
  if (!score_col %in% colnames(sce@meta.data)) {
    warning("未找到生成的列：", score_col, "；跳过该基因。")
    next
  }
  scores <- sce@meta.data[[score_col]]
  
  # 5.2 16-bin
  bin_vec <- make_bin16(scores)
  bin_col <- paste0(g_safe, "_bin16")
  sce[[bin_col]] <- bin_vec
  
  # 5.3 20/60/20
  hml_out <- make_hml(scores, prefix = g_safe)
  hml_col <- paste0(g_safe, "_3class")
  sce[[hml_col]] <- hml_out$class
  
  # 保存阈值记录
  thresholds_tbl[[g_safe]] <- data.frame(
    gene = g_safe,
    q20  = hml_out$q20,
    q80  = hml_out$q80,
    stringsAsFactors = FALSE
  )
  
  # 5.4 可视化输出目录
  out_g <- file.path(output, g_safe)
  dir.create(out_g, showWarnings = FALSE, recursive = TRUE)
  
  # 5.5 绘图：16-bin（离散）
  p_bin16 <- DimPlot(
    sce,
    reduction = pick_reduc(sce),
    group.by  = bin_col,
    pt.size   = 2,
    label     = FALSE,
    cols      = viridis(n = max(length(levels(bin_vec)), 3), option = "C")
  ) +
    theme_void(base_size = 14) +
    theme(
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(title = paste0(g_safe, " score (16-quantile bins)"), color = "Bin")
  
  ggsave(file.path(out_g, paste0(g_safe, "_bin16_DimPlot.pdf")),
         plot = p_bin16, width = 5, height = 4, device = cairo_pdf)
  
  # 5.6 绘图：H/M/L（离散）
  p_hml <- DimPlot(
    sce,
    reduction = pick_reduc(sce),
    group.by  = hml_col,
    pt.size   = 2,
    label     = FALSE,
    cols      = c("#FF0066","#CCFFFF","#0066CC")
  ) +
    theme_void(base_size = 14) +
    theme(
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.position = "right"
    ) +
    labs(title = paste0(g_safe, " score (class)"), color = "Class")
  
  ggsave(file.path(out_g, paste0(g_safe, "_HML_DimPlot.pdf")),
         plot = p_hml, width = 5, height = 4, device = cairo_pdf)
  
  # 5.7 绘图：连续 FeaturePlot（整体）
  p_cont <- FeaturePlot(
    sce,
    features = score_col,
    order    = TRUE,
    ncol     = 1,
    pt.size   = 2,
    cols     = c("#CCCCCC","#0066CC")
  ) +
    theme_void(base_size = 14) +
    theme(
      plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12),
      legend.position = "right"
    ) +
    labs(title = paste0(g_safe, " UCell score"), color = "Score")
  
  ggsave(file.path(out_g, paste0(g_safe, "_FeaturePlot.pdf")),
         plot = p_cont, width = 4, height = 3)
  
  # 5.8 绘图：连续 FeaturePlot（按 treatment 分组，如果存在）
  if (!is.null(split_by)) {
    p_cont_split <- FeaturePlot(
      sce,
      features = score_col,
      split.by = split_by,
      pt.size   = 2,
      order    = TRUE,
      ncol     = length(unique(sce[[split_by]][,1])),
      cols     = c("lightgrey", "#0066CC")
    ) +
      theme(
        strip.text   = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        legend.position = "right"
      ) +
      labs(color = paste0(g_safe, "_score"))
    ggsave(file.path(out_g, paste0(g_safe, "_FeaturePlot_", split_by, ".pdf")),
           plot = p_cont_split, width = 15, height = 3)
  }
}



#====================#
# 6) 汇总导出 + 保存对象
#====================#
# 合并阈值表
thresholds_df <- bind_rows(thresholds_tbl)
write.csv(thresholds_df, file = file.path(output, "per_gene_q20_q80_thresholds.csv"), row.names = FALSE)

# 保存带有所有新列的 Seurat 对象
saveRDS(sce, file = "celltype.rds")

# 也可导出 meta.data 方便外部审阅
write.csv(sce@meta.data, file = file.path(output, "meta_with_all_ucell_bins_classes.csv"))











############ 基于days划分 #################
# ========= STARTRAC：按基因批量计算 Ro/e，合并矩阵并统一出图 =========
suppressPackageStartupMessages({
  library(Startrac)
  library(dplyr)
  library(ggplot2)
  library(circlize)
  library(ComplexHeatmap)
})


# 设置输出目录
output <- paste(outdir, 'STARTRAC/days', sep='/')
dir.create(output, showWarnings = FALSE)


# 读取包含 <GENE>_3class 的对象（延续上面保存的）
sco  <- readRDS("celltype.rds")
data <- sco@meta.data

# 感兴趣基因（与上面一致）
genes_of_interest <- unique(c(
  "Resp18","Uchl1","Syp","Pcsk1","Scg5","Chga","Chgb","Hes1","Uchl1","Krt8"

  ))

# 关键列名（如有不同，改这里）
patient_col <- "orig.ident"
tissue_col  <- "treatment"
class_suffix <- "_3class"   # <GENE>_3class

#—— 工具函数 ——#
# 行名颜色（给“Gene|Class”中的 Class 上色，或按基因上色）
mk_gene_colors <- function(genes) {
  cols <- viridis::viridis(max(length(genes), 3), option = "D")
  stats::setNames(cols[seq_along(genes)], genes)
}
mk_class_colors <- function(classes) {
  cols <- viridis::viridis(max(length(classes), 3), option = "C")
  stats::setNames(cols[seq_along(classes)], classes)
}

# 逐基因计算 Ro/e 并收集矩阵
roe_list <- list()
skipped  <- character(0)

for (g in genes_of_interest) {
  cls_col <- paste0(g, class_suffix)
  if (!cls_col %in% colnames(data)) {
    message("[跳过] 未找到列：", cls_col)
    skipped <- c(skipped, g); next
  }
  if (!all(c(patient_col, tissue_col) %in% colnames(data))) {
    stop("找不到列：", paste(setdiff(c(patient_col, tissue_col), colnames(data)), collapse = ", "))
  }
  
  Roe <- tryCatch({
    calTissueDist(
      data,
      byPatient       = FALSE,
      colname.cluster = cls_col,
      colname.patient = patient_col,
      colname.tissue  = tissue_col,
      method          = "chisq",
      min.rowSum      = 0
    )
  }, error = function(e) {
    message("[跳过] calTissueDist 失败：", g, " | ", e$message); return(NULL)
  })
  if (is.null(Roe) || nrow(Roe) == 0 || ncol(Roe) == 0) {
    message("[跳过] Ro/e 为空：", g); skipped <- c(skipped, g); next
  }
  
  # ★★★ 只保留该基因的 high / low（去掉 medium）★★★
  keep_rows <- c(paste0(g, "_high"), paste0(g, "_low"))
  Roe <- Roe[rownames(Roe) %in% keep_rows, , drop = FALSE]
  if (nrow(Roe) == 0) {
    message("[跳过] 该基因未检测到 high/low 层级：", g)
    skipped <- c(skipped, g); next
  }
  
  # 行名加前缀：Gene|Class，列维持 tissue
  rownames(Roe) <- paste(g, rownames(Roe), sep = "|")
  roe_list[[g]] <- Roe
}

if (length(roe_list) == 0) {
  stop("没有可用的 Ro/e 矩阵（可能所有基因都被跳过）。")
}

#—— 合并所有基因的 Ro/e 矩阵（按行堆叠，列（tissue）求并集）——#
# 对齐列：把缺失的 tissue 列补 NA
all_cols <- Reduce(union, lapply(roe_list, colnames))
roe_list_aligned <- lapply(roe_list, function(M) {
  miss <- setdiff(all_cols, colnames(M))
  if (length(miss) > 0) {
    M <- cbind(M, matrix(NA_real_, nrow = nrow(M), ncol = length(miss),
                         dimnames = list(rownames(M), miss)))
  }
  M[, all_cols, drop = FALSE]
})
Roe_combined <- do.call(rbind, roe_list_aligned)  # 大矩阵：行=Gene|Class，列=tissue

# 保存合并矩阵
write.csv(Roe_combined, file = file.path(output, "COMBINED_STARTRAC_Roe_matrix_NE.csv"))

#—— 统一配色范围（跨基因）——#
mat <- as.matrix(Roe_combined)
rng <- range(mat, na.rm = TRUE)
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng, na.rm = TRUE)
col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]), c("#f6f8e6", "#f9a33e", "red"))
leg_at  <- pretty(rng, n = 5)
leg_lab <- formatC(leg_at, format = "f", digits = 2)

#—— 行注释：按“基因”上色；同时保留原 Class 信息——#
row_gene  <- sub("\\|.*$", "", rownames(mat))
row_class <- sub("^.*?\\|", "", rownames(mat))
gene_levels  <- unique(row_gene)
class_levels <- unique(row_class)
gene_cols    <- mk_gene_colors(gene_levels)
class_cols   <- mk_class_colors(class_levels)

row_anno <- ComplexHeatmap::rowAnnotation(
  Gene  = row_gene,
  Class = row_class,
  col = list(
    Gene  = gene_cols,
    Class = class_cols
  ),
  gp = grid::gpar(col = NA)
)

#========================#
# 图1：带数值的“合并”热图
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_value_NE.pdf"), width = 7, height = max(4, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e value",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),   # 行很多时适当调小
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = leg_at, labels = leg_lab,
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 12)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    if (!is.na(v)) grid::grid.text(sprintf("%.2f", v), x, y, gp = grid::gpar(fontsize = 8, col = "black"))
  }
)
dev.off()

#========================#
# 图2：自定义“+++”符号热图（合并版）
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_symbols_NE.pdf"), width = 7, height = max(4, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = c(0, max(mat, na.rm = TRUE)), labels = c("0", "Max."),
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 12)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    value <- mat[i, j]
    symbol <- if (is.na(value) || value == 0) {
      "−"
    } else if (value > 0 & value < 0.2) {
      "+/−"
    } else if (value >= 0.2 & value <= 0.8) {
      "+"
    } else if (value > 0.8 & value <= 1) {
      "++"
    } else {
      "+++"
    }
    txt_col <- ifelse(mean(grDevices::col2rgb(fill)) > 127, "black", "white")
    grid::grid.text(symbol, x, y, gp = grid::gpar(fontsize = 9, col = txt_col))
  }
)
dev.off()

#========================#
# 图3：合并气泡图（按基因分面）
#========================#
bubble_df <- as.data.frame(as.table(mat)) %>%
  dplyr::rename(Tissue = Var2, Row = Var1, Roe = Freq) %>%
  dplyr::mutate(
    Gene     = sub("\\|.*$", "", Row),
    CellType = sub("^.*?\\|", "", Row),
    Enrichment = dplyr::case_when(
      is.na(Roe) ~ NA_character_,
      Roe >= 1   ~ "Enrichment",
      TRUE       ~ "Depletion"
    )
  ) %>%
  dplyr::filter(!is.na(Roe))

pdf(file.path(output, "COMBINED_STARTRAC_Roe_bubble_NE.pdf"),
    width = 10, height = max(4, ceiling(length(unique(bubble_df$Gene))/3) * 3))
ggplot(bubble_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(aes(size = Roe, color = Enrichment)) +
  scale_size_continuous(name = "Ro/e", breaks = c(0.5, 1.0, 1.5), range = c(1, 7)) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),
    labels = c("Enrichment", "Depletion"),
    guide  = guide_legend(override.aes = list(size = 4))
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  scale_y_discrete(limits = function(x) rev(x)) +
  labs(
    title = "STARTRAC - Combined Ro/e (by gene facets)",
    x = "Tissue Group", y = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.box      = "horizontal",
    strip.text      = element_text(face = "bold")
  )
dev.off()

if (length(skipped) > 0) {
  message("以下基因因缺列或空结果被跳过：", paste(skipped, collapse = ", "))
} else {
  message("所有基因均已纳入合并分析。")
}






############ 基于days划分 #################
# ========= STARTRAC：按基因批量计算 Ro/e，合并矩阵并统一出图 =========
suppressPackageStartupMessages({
  library(Startrac)
  library(dplyr)
  library(ggplot2)
  library(circlize)
  library(ComplexHeatmap)
})


# 设置输出目录
output <- paste(outdir, 'STARTRAC/days', sep='/')
dir.create(output, showWarnings = FALSE)


# 读取包含 <GENE>_3class 的对象（延续上面保存的）
sco  <- readRDS("celltype.rds")
data <- sco@meta.data

# 感兴趣基因（与上面一致）
genes_of_interest <- unique(c(
  "Cbr2","Tppp3","Scgb3a2","Sftpc","Ager"
))

# 关键列名（如有不同，改这里）
patient_col <- "orig.ident"
tissue_col  <- "treatment"
class_suffix <- "_3class"   # <GENE>_3class

#—— 工具函数 ——#
# 行名颜色（给“Gene|Class”中的 Class 上色，或按基因上色）
mk_gene_colors <- function(genes) {
  cols <- viridis::viridis(max(length(genes), 3), option = "D")
  stats::setNames(cols[seq_along(genes)], genes)
}
mk_class_colors <- function(classes) {
  cols <- viridis::viridis(max(length(classes), 3), option = "C")
  stats::setNames(cols[seq_along(classes)], classes)
}

# 逐基因计算 Ro/e 并收集矩阵
roe_list <- list()
skipped  <- character(0)

for (g in genes_of_interest) {
  cls_col <- paste0(g, class_suffix)
  if (!cls_col %in% colnames(data)) {
    message("[跳过] 未找到列：", cls_col)
    skipped <- c(skipped, g); next
  }
  if (!all(c(patient_col, tissue_col) %in% colnames(data))) {
    stop("找不到列：", paste(setdiff(c(patient_col, tissue_col), colnames(data)), collapse = ", "))
  }
  
  Roe <- tryCatch({
    calTissueDist(
      data,
      byPatient       = FALSE,
      colname.cluster = cls_col,
      colname.patient = patient_col,
      colname.tissue  = tissue_col,
      method          = "chisq",
      min.rowSum      = 0
    )
  }, error = function(e) {
    message("[跳过] calTissueDist 失败：", g, " | ", e$message); return(NULL)
  })
  if (is.null(Roe) || nrow(Roe) == 0 || ncol(Roe) == 0) {
    message("[跳过] Ro/e 为空：", g); skipped <- c(skipped, g); next
  }
  
  # ★★★ 只保留该基因的 high / low（去掉 medium）★★★
  keep_rows <- c(paste0(g, "_high"), paste0(g, "_low"))
  Roe <- Roe[rownames(Roe) %in% keep_rows, , drop = FALSE]
  if (nrow(Roe) == 0) {
    message("[跳过] 该基因未检测到 high/low 层级：", g)
    skipped <- c(skipped, g); next
  }
  
  # 行名加前缀：Gene|Class，列维持 tissue
  rownames(Roe) <- paste(g, rownames(Roe), sep = "|")
  roe_list[[g]] <- Roe
}

if (length(roe_list) == 0) {
  stop("没有可用的 Ro/e 矩阵（可能所有基因都被跳过）。")
}

#—— 合并所有基因的 Ro/e 矩阵（按行堆叠，列（tissue）求并集）——#
# 对齐列：把缺失的 tissue 列补 NA
all_cols <- Reduce(union, lapply(roe_list, colnames))
roe_list_aligned <- lapply(roe_list, function(M) {
  miss <- setdiff(all_cols, colnames(M))
  if (length(miss) > 0) {
    M <- cbind(M, matrix(NA_real_, nrow = nrow(M), ncol = length(miss),
                         dimnames = list(rownames(M), miss)))
  }
  M[, all_cols, drop = FALSE]
})
Roe_combined <- do.call(rbind, roe_list_aligned)  # 大矩阵：行=Gene|Class，列=tissue

# 保存合并矩阵
write.csv(Roe_combined, file = file.path(output, "COMBINED_STARTRAC_Roe_matrix_Non-NE.csv"))

#—— 统一配色范围（跨基因）——#
mat <- as.matrix(Roe_combined)
rng <- range(mat, na.rm = TRUE)
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng, na.rm = TRUE)
col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]), c("#f6f8e6", "#f9a33e", "red"))
leg_at  <- pretty(rng, n = 5)
leg_lab <- formatC(leg_at, format = "f", digits = 2)

#—— 行注释：按“基因”上色；同时保留原 Class 信息——#
row_gene  <- sub("\\|.*$", "", rownames(mat))
row_class <- sub("^.*?\\|", "", rownames(mat))
gene_levels  <- unique(row_gene)
class_levels <- unique(row_class)
gene_cols    <- mk_gene_colors(gene_levels)
class_cols   <- mk_class_colors(class_levels)

row_anno <- ComplexHeatmap::rowAnnotation(
  Gene  = row_gene,
  Class = row_class,
  col = list(
    Gene  = gene_cols,
    Class = class_cols
  ),
  gp = grid::gpar(col = NA)
)

#========================#
# 图1：带数值的“合并”热图
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_value_Non-NE.pdf"), width = 7, height = max(3, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e value",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),   # 行很多时适当调小
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = leg_at, labels = leg_lab,
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 12)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    if (!is.na(v)) grid::grid.text(sprintf("%.2f", v), x, y, gp = grid::gpar(fontsize = 8, col = "black"))
  }
)
dev.off()

#========================#
# 图2：自定义“+++”符号热图（合并版）
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_symbols_Non-NE.pdf"), width = 7, height = max(3, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = c(0, max(mat, na.rm = TRUE)), labels = c("0", "Max."),
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 12)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    value <- mat[i, j]
    symbol <- if (is.na(value) || value == 0) {
      "−"
    } else if (value > 0 & value < 0.2) {
      "+/−"
    } else if (value >= 0.2 & value <= 0.8) {
      "+"
    } else if (value > 0.8 & value <= 1) {
      "++"
    } else {
      "+++"
    }
    txt_col <- ifelse(mean(grDevices::col2rgb(fill)) > 127, "black", "white")
    grid::grid.text(symbol, x, y, gp = grid::gpar(fontsize = 9, col = txt_col))
  }
)
dev.off()

#========================#
# 图3：合并气泡图（按基因分面）
#========================#
bubble_df <- as.data.frame(as.table(mat)) %>%
  dplyr::rename(Tissue = Var2, Row = Var1, Roe = Freq) %>%
  dplyr::mutate(
    Gene     = sub("\\|.*$", "", Row),
    CellType = sub("^.*?\\|", "", Row),
    Enrichment = dplyr::case_when(
      is.na(Roe) ~ NA_character_,
      Roe >= 1   ~ "Enrichment",
      TRUE       ~ "Depletion"
    )
  ) %>%
  dplyr::filter(!is.na(Roe))

pdf(file.path(output, "COMBINED_STARTRAC_Roe_bubble_Non-NE.pdf"),
    width = 10, height = max(4, ceiling(length(unique(bubble_df$Gene))/3) * 3))
ggplot(bubble_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(aes(size = Roe, color = Enrichment)) +
  scale_size_continuous(name = "Ro/e", breaks = c(0.5, 1.0, 1.5), range = c(1, 7)) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),
    labels = c("Enrichment", "Depletion"),
    guide  = guide_legend(override.aes = list(size = 4))
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  scale_y_discrete(limits = function(x) rev(x)) +
  labs(
    title = "STARTRAC - Combined Ro/e (by gene facets)",
    x = "Tissue Group", y = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.box      = "horizontal",
    strip.text      = element_text(face = "bold")
  )
dev.off()

if (length(skipped) > 0) {
  message("以下基因因缺列或空结果被跳过：", paste(skipped, collapse = ", "))
} else {
  message("所有基因均已纳入合并分析。")
}









########### 基于Cbr2划分 ###################
# ========= STARTRAC：按基因批量计算 Ro/e，合并矩阵并统一出图 =========

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


# 设置输出目录
output <- paste(outdir, 'STARTRAC/NE vs Non-NE(Cbr2)', sep='/')
dir.create(output, showWarnings = FALSE)


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
  colname.cluster = "Cbr2_3class",
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
pdf(file.path(output, "STARTRAC_Roe_value_Crb2.pdf"), width = 6, height = 3)
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
pdf(file.path(output, "STARTRAC_Roe_colored_Crb2.pdf"), width = 6, height = 3)
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

pdf(file.path(output, "STARTRAC_Roe_bubble_Crb2.pdf"), width = 5, height = 3)
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










########### 基于Cbr2划分 ###################
# ========= STARTRAC：按基因批量计算 Ro/e，合并矩阵并统一出图 =========

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


# 设置输出目录
output <- paste(outdir, 'STARTRAC/NE vs Non-NE(Cbr2)', sep='/')
dir.create(output, showWarnings = FALSE)


# 读取包含 <GENE>_3class 的对象（延续上面保存的）
sco  <- readRDS("celltype.rds")
data <- sco@meta.data

# 感兴趣基因（与上面一致）
genes_of_interest <- unique(c(
    "Resp18","Uchl1","Syp","Pcsk1","Scg5","Chga","Chgb","Hes1","Uchl1","Krt8"
))

# 关键列名（如有不同，改这里）
patient_col <- "orig.ident"
tissue_col  <- "Cbr2_3class"
class_suffix <- "_3class"   # <GENE>_3class


#—— 工具函数 ——#
# 行名颜色（给“Gene|Class”中的 Class 上色，或按基因上色）
mk_gene_colors <- function(genes) {
  cols <- viridis::viridis(max(length(genes), 3), option = "D")
  stats::setNames(cols[seq_along(genes)], genes)
}
mk_class_colors <- function(classes) {
  cols <- viridis::viridis(max(length(classes), 3), option = "C")
  stats::setNames(cols[seq_along(classes)], classes)
}

# 逐基因计算 Ro/e 并收集矩阵
roe_list <- list()
skipped  <- character(0)

for (g in genes_of_interest) {
  cls_col <- paste0(g, class_suffix)
  if (!cls_col %in% colnames(data)) {
    message("[跳过] 未找到列：", cls_col)
    skipped <- c(skipped, g); next
  }
  if (!all(c(patient_col, tissue_col) %in% colnames(data))) {
    stop("找不到列：", paste(setdiff(c(patient_col, tissue_col), colnames(data)), collapse = ", "))
  }
  
  Roe <- tryCatch({
    calTissueDist(
      data,
      byPatient       = FALSE,
      colname.cluster = cls_col,
      colname.patient = patient_col,
      colname.tissue  = tissue_col,
      method          = "chisq",
      min.rowSum      = 0
    )
  }, error = function(e) {
    message("[跳过] calTissueDist 失败：", g, " | ", e$message); return(NULL)
  })
  if (is.null(Roe) || nrow(Roe) == 0 || ncol(Roe) == 0) {
    message("[跳过] Ro/e 为空：", g); skipped <- c(skipped, g); next
  }
  
  # ★★★ 只保留该基因的 high / low（去掉 medium）★★★
  keep_rows <- c(paste0(g, "_high"), paste0(g, "_low"))
  Roe <- Roe[rownames(Roe) %in% keep_rows, , drop = FALSE]
  if (nrow(Roe) == 0) {
    message("[跳过] 该基因未检测到 high/low 层级：", g)
    skipped <- c(skipped, g); next
  }
  
  # 行名加前缀：Gene|Class，列维持 tissue
  rownames(Roe) <- paste(g, rownames(Roe), sep = "|")
  roe_list[[g]] <- Roe
}

if (length(roe_list) == 0) {
  stop("没有可用的 Ro/e 矩阵（可能所有基因都被跳过）。")
}

#—— 合并所有基因的 Ro/e 矩阵（按行堆叠，列（tissue）求并集）——#
# 对齐列：把缺失的 tissue 列补 NA
all_cols <- Reduce(union, lapply(roe_list, colnames))
roe_list_aligned <- lapply(roe_list, function(M) {
  miss <- setdiff(all_cols, colnames(M))
  if (length(miss) > 0) {
    M <- cbind(M, matrix(NA_real_, nrow = nrow(M), ncol = length(miss),
                         dimnames = list(rownames(M), miss)))
  }
  M[, all_cols, drop = FALSE]
})
Roe_combined <- do.call(rbind, roe_list_aligned)  # 大矩阵：行=Gene|Class，列=tissue

Roe_combined <- Roe_combined[, setdiff(colnames(Roe_combined), "Cbr2_medium"), drop = FALSE]

# 保存合并矩阵
write.csv(Roe_combined, file = file.path(output, "COMBINED_STARTRAC_Roe_matrix_NE.csv"))

#—— 统一配色范围（跨基因）——#
mat <- as.matrix(Roe_combined)
rng <- range(mat, na.rm = TRUE)
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng, na.rm = TRUE)
col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]), c("#f6f8e6", "#f9a33e", "red"))
leg_at  <- pretty(rng, n = 5)
leg_lab <- formatC(leg_at, format = "f", digits = 2)

#—— 行注释：按“基因”上色；同时保留原 Class 信息——#
row_gene  <- sub("\\|.*$", "", rownames(mat))
row_class <- sub("^.*?\\|", "", rownames(mat))
gene_levels  <- unique(row_gene)
class_levels <- unique(row_class)
gene_cols    <- mk_gene_colors(gene_levels)
class_cols   <- mk_class_colors(class_levels)

row_anno <- ComplexHeatmap::rowAnnotation(
  Gene  = row_gene,
  Class = row_class,
  col = list(
    Gene  = gene_cols,
    Class = class_cols
  ),
  gp = grid::gpar(col = NA)
)

#========================#
# 图1：带数值的“合并”热图
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_value_NE.pdf"), width = 7, height = max(4, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e value",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),   # 行很多时适当调小
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = leg_at, labels = leg_lab,
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 14)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    if (!is.na(v)) grid::grid.text(sprintf("%.2f", v), x, y, gp = grid::gpar(fontsize = 8, col = "black"))
  }
)
dev.off()

#========================#
# 图2：自定义“+++”符号热图（合并版）
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_symbols_NE.pdf"), width = 7, height = max(4, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = c(0, max(mat, na.rm = TRUE)), labels = c("0", "Max."),
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 14)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    value <- mat[i, j]
    symbol <- if (is.na(value) || value == 0) {
      "−"
    } else if (value > 0 & value < 0.2) {
      "+/−"
    } else if (value >= 0.2 & value <= 0.8) {
      "+"
    } else if (value > 0.8 & value <= 1) {
      "++"
    } else {
      "+++"
    }
    txt_col <- ifelse(mean(grDevices::col2rgb(fill)) > 127, "black", "white")
    grid::grid.text(symbol, x, y, gp = grid::gpar(fontsize = 9, col = txt_col))
  }
)
dev.off()

#========================#
# 图3：合并气泡图（按基因分面）
#========================#
bubble_df <- as.data.frame(as.table(mat)) %>%
  dplyr::rename(Tissue = Var2, Row = Var1, Roe = Freq) %>%
  dplyr::mutate(
    Gene     = sub("\\|.*$", "", Row),
    CellType = sub("^.*?\\|", "", Row),
    Enrichment = dplyr::case_when(
      is.na(Roe) ~ NA_character_,
      Roe >= 1   ~ "Enrichment",
      TRUE       ~ "Depletion"
    )
  ) %>%
  dplyr::filter(!is.na(Roe))

pdf(file.path(output, "COMBINED_STARTRAC_Roe_bubble_NE.pdf"),
    width = 10, height = max(4, ceiling(length(unique(bubble_df$Gene))/3) * 3))
ggplot(bubble_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(aes(size = Roe, color = Enrichment)) +
  scale_size_continuous(name = "Ro/e", breaks = c(0.5, 1.0, 1.5), range = c(1, 7)) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),
    labels = c("Enrichment", "Depletion"),
    guide  = guide_legend(override.aes = list(size = 4))
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  scale_y_discrete(limits = function(x) rev(x)) +
  labs(
    title = "STARTRAC - Combined Ro/e (by gene facets)",
    x = "Tissue Group", y = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.box      = "horizontal",
    strip.text      = element_text(face = "bold")
  )
dev.off()

if (length(skipped) > 0) {
  message("以下基因因缺列或空结果被跳过：", paste(skipped, collapse = ", "))
} else {
  message("所有基因均已纳入合并分析。")
}















########### 基于Cbr2划分 ###################
# ========= STARTRAC：按基因批量计算 Ro/e，合并矩阵并统一出图 =========

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


# 设置输出目录
output <- paste(outdir, 'STARTRAC/NE vs Non-NE(Cbr2)', sep='/')
dir.create(output, showWarnings = FALSE)


# 读取包含 <GENE>_3class 的对象（延续上面保存的）
sco  <- readRDS("celltype.rds")
data <- sco@meta.data

# 感兴趣基因（与上面一致）
genes_of_interest <- unique(c(
  "Tppp3","Scgb3a2","Sftpc","Ager"
))

# 关键列名（如有不同，改这里）
patient_col <- "orig.ident"
tissue_col  <- "Cbr2_3class"
class_suffix <- "_3class"   # <GENE>_3class

#✅【新增】仅保留 Cbr2_high / Cbr2_low 样本
data <- data %>% dplyr::filter(!!sym(tissue_col) %in% c("Cbr2_high", "Cbr2_low"))

#—— 工具函数 ——#
# 行名颜色（给“Gene|Class”中的 Class 上色，或按基因上色）
mk_gene_colors <- function(genes) {
  cols <- viridis::viridis(max(length(genes), 3), option = "D")
  stats::setNames(cols[seq_along(genes)], genes)
}
mk_class_colors <- function(classes) {
  cols <- viridis::viridis(max(length(classes), 3), option = "C")
  stats::setNames(cols[seq_along(classes)], classes)
}

# 逐基因计算 Ro/e 并收集矩阵
roe_list <- list()
skipped  <- character(0)

for (g in genes_of_interest) {
  cls_col <- paste0(g, class_suffix)
  if (!cls_col %in% colnames(data)) {
    message("[跳过] 未找到列：", cls_col)
    skipped <- c(skipped, g); next
  }
  if (!all(c(patient_col, tissue_col) %in% colnames(data))) {
    stop("找不到列：", paste(setdiff(c(patient_col, tissue_col), colnames(data)), collapse = ", "))
  }
  
  Roe <- tryCatch({
    calTissueDist(
      data,
      byPatient       = FALSE,
      colname.cluster = cls_col,
      colname.patient = patient_col,
      colname.tissue  = tissue_col,
      method          = "chisq",
      min.rowSum      = 0
    )
  }, error = function(e) {
    message("[跳过] calTissueDist 失败：", g, " | ", e$message); return(NULL)
  })
  if (is.null(Roe) || nrow(Roe) == 0 || ncol(Roe) == 0) {
    message("[跳过] Ro/e 为空：", g); skipped <- c(skipped, g); next
  }
  
  # ★★★ 只保留该基因的 high / low（去掉 medium）★★★
  keep_rows <- c(paste0(g, "_high"), paste0(g, "_low"))
  Roe <- Roe[rownames(Roe) %in% keep_rows, , drop = FALSE]
  if (nrow(Roe) == 0) {
    message("[跳过] 该基因未检测到 high/low 层级：", g)
    skipped <- c(skipped, g); next
  }
  
  # 行名加前缀：Gene|Class，列维持 tissue
  rownames(Roe) <- paste(g, rownames(Roe), sep = "|")
  roe_list[[g]] <- Roe
}

if (length(roe_list) == 0) {
  stop("没有可用的 Ro/e 矩阵（可能所有基因都被跳过）。")
}

#—— 合并所有基因的 Ro/e 矩阵（按行堆叠，列（tissue）求并集）——#
# 对齐列：把缺失的 tissue 列补 NA
all_cols <- Reduce(union, lapply(roe_list, colnames))
roe_list_aligned <- lapply(roe_list, function(M) {
  miss <- setdiff(all_cols, colnames(M))
  if (length(miss) > 0) {
    M <- cbind(M, matrix(NA_real_, nrow = nrow(M), ncol = length(miss),
                         dimnames = list(rownames(M), miss)))
  }
  M[, all_cols, drop = FALSE]
})
Roe_combined <- do.call(rbind, roe_list_aligned)  # 大矩阵：行=Gene|Class，列=tissue

## 删除"Cbr2_medium"
Roe_combined <- Roe_combined[, setdiff(colnames(Roe_combined), "Cbr2_medium"), drop = FALSE]

# 保存合并矩阵
write.csv(Roe_combined, file = file.path(output, "COMBINED_STARTRAC_Roe_matrix_No-NE.csv"))

#—— 统一配色范围（跨基因）——#
mat <- as.matrix(Roe_combined)
rng <- range(mat, na.rm = TRUE)
mid <- if (1 >= rng[1] && 1 <= rng[2]) 1 else mean(rng, na.rm = TRUE)
col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]), c("#f6f8e6", "#f9a33e", "red"))
leg_at  <- pretty(rng, n = 5)
leg_lab <- formatC(leg_at, format = "f", digits = 2)

#—— 行注释：按“基因”上色；同时保留原 Class 信息——#
row_gene  <- sub("\\|.*$", "", rownames(mat))
row_class <- sub("^.*?\\|", "", rownames(mat))
gene_levels  <- unique(row_gene)
class_levels <- unique(row_class)
gene_cols    <- mk_gene_colors(gene_levels)
class_cols   <- mk_class_colors(class_levels)

row_anno <- ComplexHeatmap::rowAnnotation(
  Gene  = row_gene,
  Class = row_class,
  col = list(
    Gene  = gene_cols,
    Class = class_cols
  ),
  gp = grid::gpar(col = NA)
)

#========================#
# 图1：带数值的“合并”热图
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_value_No-NE.pdf"), width = 6.5, height = max(3, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e value",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),   # 行很多时适当调小
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = leg_at, labels = leg_lab,
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 12)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- mat[i, j]
    if (!is.na(v)) grid::grid.text(sprintf("%.2f", v), x, y, gp = grid::gpar(fontsize = 8, col = "black"))
  }
)
dev.off()

#========================#
# 图2：自定义“+++”符号热图（合并版）
#========================#
pdf(file.path(output, "COMBINED_STARTRAC_Roe_symbols_No-NE.pdf"), width = 6.5, height = max(3, nrow(mat) * 0.18))
ComplexHeatmap::Heatmap(
  mat,
  name = "Ro/e",
  col  = col_fun,
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "right", column_names_side = "bottom",
  show_row_names = TRUE,  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 12),
  heatmap_legend_param = list(
    at = c(0, max(mat, na.rm = TRUE)), labels = c("0", "Max."),
    title_gp = grid::gpar(fontsize = 14),
    labels_gp = grid::gpar(fontsize = 12)
  ),
  left_annotation = row_anno,
  cell_fun = function(j, i, x, y, w, h, fill) {
    value <- mat[i, j]
    symbol <- if (is.na(value) || value == 0) {
      "−"
    } else if (value > 0 & value < 0.2) {
      "+/−"
    } else if (value >= 0.2 & value <= 0.8) {
      "+"
    } else if (value > 0.8 & value <= 1) {
      "++"
    } else {
      "+++"
    }
    txt_col <- ifelse(mean(grDevices::col2rgb(fill)) > 127, "black", "white")
    grid::grid.text(symbol, x, y, gp = grid::gpar(fontsize = 9, col = txt_col))
  }
)
dev.off()

#========================#
# 图3：合并气泡图（按基因分面）
#========================#
bubble_df <- as.data.frame(as.table(mat)) %>%
  dplyr::rename(Tissue = Var2, Row = Var1, Roe = Freq) %>%
  dplyr::mutate(
    Gene     = sub("\\|.*$", "", Row),
    CellType = sub("^.*?\\|", "", Row),
    Enrichment = dplyr::case_when(
      is.na(Roe) ~ NA_character_,
      Roe >= 1   ~ "Enrichment",
      TRUE       ~ "Depletion"
    )
  ) %>%
  dplyr::filter(!is.na(Roe))

pdf(file.path(output, "COMBINED_STARTRAC_Roe_bubble_No-NE.pdf"),
    width = 10, height = max(4, ceiling(length(unique(bubble_df$Gene))/3) * 3))
ggplot(bubble_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(aes(size = Roe, color = Enrichment)) +
  scale_size_continuous(name = "Ro/e", breaks = c(0.5, 1.0, 1.5), range = c(1, 7)) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#2E75B6", "Depletion" = "#E36C8C"),
    labels = c("Enrichment", "Depletion"),
    guide  = guide_legend(override.aes = list(size = 4))
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  scale_y_discrete(limits = function(x) rev(x)) +
  labs(
    title = "STARTRAC - Combined Ro/e (by gene facets)",
    x = "Tissue Group", y = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.box      = "horizontal",
    strip.text      = element_text(face = "bold")
  )
dev.off()

if (length(skipped) > 0) {
  message("以下基因因缺列或空结果被跳过：", paste(skipped, collapse = ", "))
} else {
  message("所有基因均已纳入合并分析。")
}








###############差异基因在不同组中表达##################
library(Seurat)
library(tidyverse)
library(ggsci)

# 设置输出目录
output <- paste(outdir, 'STARTRAC/NE vs Non-NE(Cbr2)', sep='/')
dir.create(output, showWarnings = FALSE)

# 读取包含 <GENE>_3class 的对象（延续上面保存的）
ScRNA <- readRDS("celltype.rds")

genes <- c(
  "Resp18","Uchl1","Syp","Pcsk1","Scg5","Chga","Chgb","Calca","Calcb","Sez6l2",
  "Uchl1","Ascl1","Krt8",
  "Cbr2","Foxj1","Scgb3a2","Lamp3",'Lpcat1',"Sftpc", "Ager"
)  

# 过滤在表达矩阵中实际存在的基因
genes <- genes[genes %in% rownames(ScRNA)]

# 删除 Cbr2_medium 类别
# 方法1: 直接过滤掉 Cbr2_medium 的细胞
cells_to_keep <- which(ScRNA$Cbr2_3class != "Cbr2_medium")
ScRNA_filtered <- subset(ScRNA, cells = colnames(ScRNA)[cells_to_keep])

# 重新设置因子水平，移除空的 Cbr2_medium 水平
ScRNA_filtered$Cbr2_3class <- factor(ScRNA_filtered$Cbr2_3class)

# 使用过滤后的对象继续分析
ScRNA <- ScRNA_filtered

# 提取treatment信息
treatment_groups <- levels(ScRNA$Cbr2_3class)
print("处理组别:")
print(treatment_groups)

# 获取表达矩阵
expr_matrix <- GetAssayData(ScRNA, slot = "data")[genes, ]

# 计算平均表达量
avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "Cbr2_3class")$RNA
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
                                     group.by = "Cbr2_3class",
                                     assays = "RNA")$RNA

print("处理组别（平均表达量计算后）:")
print(colnames(avg_exp))

# 计算Tumor组相对于其他组的标准化表达量
# 这里我们计算Tumor组表达量减去其他组平均表达量
other_groups <- setdiff(colnames(avg_exp), "Cbr2_high")
tumor_normalized <- avg_exp[,"Cbr2_high"] - rowMeans(avg_exp[,other_groups, drop = FALSE])

# 根据标准化表达量排序基因
# 先按是否为正值分组，再按绝对值大小排序
genes_positive <- names(sort(tumor_normalized[tumor_normalized > 0], decreasing = TRUE))
genes_negative <- names(sort(tumor_normalized[tumor_normalized <= 0], decreasing = FALSE))

# 合并基因顺序
genes_ordered <- c(genes_positive, genes_negative)

# 确保所有要绘制的基因都在排序列表中
genes_to_plot <- intersect(genes_ordered, genes_to_plot)

# 绘制基因表达量点图（根据treatment分组）
plot <- DotPlot(ScRNA, features = genes_to_plot, group.by = "Cbr2_3class") + 
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






############### 差异基因在不同组中表达（NE子集比例图） ##################
library(Seurat)
library(tidyverse)
library(ggsci)
library(Matrix)
library(scales)

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
         "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
         "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
         "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
         "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
         "#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#4F6272','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# ----------------------
# 基础设置与读入对象
# ----------------------
output <- file.path(outdir, "STARTRAC/NE vs Non-NE(Cbr2)")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

ScRNA <- readRDS("celltype.rds")

genes <- c(
  "Calca","Scg5","Chgb","Scg2","Pcsk1n","Chga","Scg3","Cartpt",
  "Dbh","Vgf","Inha","Agt","Calcb","Iapp","Inhbb","Adcyap1",
  "Pomc","Adm","Agrp","Crh","Edn3","Gal","Ghrh","Gnrh1",
  "Nmb","Nppb","Npw","Nppa","Uts2","Nts","Ucn2","Ucn3",
  "Inhba","Prok2","Edn1","Oxt","Ghrl","Lhb","Nppc"
)

# 过滤仅保留表达矩阵中存在的基因
genes <- genes[genes %in% rownames(ScRNA)]

# ----------------------
# 选择 Cbr2_low 并重命名为 NE
# ----------------------
cells_ne <- colnames(ScRNA)[ScRNA$Cbr2_3class == "Cbr2_low"]

# 新增分组列（仅用于标注）
ScRNA$Group <- ifelse(ScRNA$Cbr2_3class == "Cbr2_low", "NE", ScRNA$Cbr2_3class)
# View(ScRNA@meta.data)

# ----------------------
# 计算 NE 内各基因的表达比例与“表达细胞中的 LN(CPM+1) 平均值”
# ----------------------
cts <- GetAssayData(ScRNA, slot = "counts")[genes, cells_ne, drop = FALSE]

expr_pos <- cts > 0
prop_ne  <- Matrix::rowMeans(expr_pos)

lib_size <- Matrix::colSums(cts)
lib_size[lib_size == 0] <- 1
cpm <- Matrix::t(Matrix::t(cts) / lib_size * 1e6)
logcpm <- log1p(cpm)

n_expr <- Matrix::rowSums(expr_pos)
sum_log_expr <- Matrix::rowSums(logcpm * (expr_pos))
mean_log_expr <- sum_log_expr / pmax(n_expr, 1)
mean_log_expr[n_expr == 0] <- NA_real_

df_plot <- tibble(
  gene = rownames(cts),
  proportion = as.numeric(prop_ne),
  mean_log_cpm1_expressing = as.numeric(mean_log_expr),
  n_expressing_cells = as.numeric(n_expr)
) %>%
  arrange(desc(proportion))

# ----------------------
# 绘图（将均值标签放在 y = 1.0；增大坐标字体；修改颜色）
# ----------------------
p <- ggplot(
  df_plot,
  aes(x = reorder(gene, -proportion), y = proportion)
) +
  geom_col(aes(fill = gene), width = 0.75, show.legend = FALSE) +

  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.9),                 # 顶部留白，确保 y=1.0 的标签不被裁切
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "",
    y = "% NE cells",
    title = ""
  ) +
  scale_fill_manual(values = col) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 16, angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  )

# 根据基因数量自适应宽度
pdf_file <- file.path(output, "NE_genes_expression_proportion_with_logCPM.pdf")
ggsave(pdf_file, p, width = max(6, nrow(df_plot) * 0.30), height = 5.0, units = "in")

# 导出统计表
write_csv(df_plot, file.path(output, "NE_genes_expression_stats.csv"))












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

file_path <- file.path(outdir, "ScRNA（分群后）.rds")
data1 <- readRDS(file_path)
summary(data1$treatment)


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
annotation <- data.frame(phenotype = sce3@meta.data$treatment) %>% 
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
p1 <- ggboxplot(cytotrace2_result_sce@meta.data, x="treatment", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",#轮廓颜色
                fill="treatment",#填充
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

###指定组比较
my_comparisons <- list(c( "PNEC","AT2"), c("T", "un"),c("Myeloid", "un"))

# 添加统计比较
p1 <- p1 + stat_compare_means(comparisons = my_comparisons, 
                              method = "wilcox.test", 
                              label = "p.signif") # 添加显著性标签
# 使用ggsave将图形保存为PDF
ggsave(file.path(output, "CytoTRACE2_Boxplot_byPheno.pdf"), plot = p1,w=4,h=3)











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
output <- paste(outdir,'差异分析(NE vs Non-NE)', sep='/')
dir.create(output, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)
colnames(scRNAsub@meta.data)


# 寻找 Res 和 Sen 组之间的差异基因
logFCfilter <- 0.25        # 定义 log2FC 过滤值
adjPvalFilter <- 0.05   # 定义矫正后 P 值过滤值

# 寻找 Epi_cisplatin_res 和 Epi_other 组之间的差异基因
scRNAsub.cluster.markers <- FindMarkers(object = scRNAsub, 
                                        ident.1 = "Cbr2_low",
                                        ident.2 =  "Cbr2_high",
                                        group.by = "Cbr2_3class", 
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
  geom_text_repel(data = top_genes_upregulated, aes(label = top_genes_upregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  geom_text_repel(data = top_genes_downregulated, aes(label = top_genes_downregulated$gene), size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  #geom_text_repel(data = interested_genes, aes(label = interested_genes$gene), size = 5, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  theme_classic() +
  labs(title = "Cbr2_low vs Cbr2_high", 
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
library(org.Mm.eg.db) # 小鼠数据库
library(org.Hs.eg.db) # 小鼠数据库
library(clusterProfiler)
library(enrichplot)
library(DOSE)

scRNAsub.cluster.markers <- readRDS(file.path(output, "ScRNA.sig.markers.rds"))

# 从差异分析中获取差异基因列表
deg <- scRNAsub.cluster.markers[, c('avg_log2FC', 'p_val_adj')]
colnames(deg) <- c('log2FoldChange', 'pvalue')  # 更改列名

# SYMBOL转换为ENTREZID
gene <- bitr(rownames(deg), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# 匹配logFC信息
gene$logFC <- deg$log2FoldChange[match(gene$SYMBOL, rownames(deg))]

# 构建genelist
geneList <- gene$logFC
names(geneList) <- gene$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)  # 根据logFC降序排序

# GSEA分析（KEGG通路）
kk_gse <- gseKEGG(geneList = geneList,
                  organism =  "mmu",  
                  ## 'mmu',  # 小鼠
                  nPerm = 1000,       ##置换次数：表示随机打乱基因集1000次进行模拟
                  minGSSize = 10,     ##用于富集分析的基因集至少需要包含10个基因。
                  pvalueCutoff = 0.25,
                  verbose = FALSE)

# 将ENTREZID转换为可读的SYMBOL名称
kk_gse <- DOSE::setReadable(kk_gse, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID')
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

# 将geneID从ENTREZID转为SYMBOL
go_up@result$geneID <- sapply(strsplit(go_up@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
  paste(symbols, collapse = "/")
})

go_down@result$geneID <- sapply(strsplit(go_down@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = ids, columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL
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
kegg_up <- enrichKEGG(gene = gene_up_entrez, organism = 'mmu', pAdjustMethod = "BH", pvalueCutoff = 0.1)
kegg_down <- enrichKEGG(gene = gene_down_entrez, organism = 'mmu', pAdjustMethod = "BH", pvalueCutoff = 0.1)

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

#file_path <- file.path(outdir, "celltype.rds")
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
ScRNA <- readRDS(file_path)


# 设置T细胞激活相关基因
cellmarker <- c(
  "Resp18","Cbr2","Notch2","Hes1",
  
  "Krt8" ,"Foxj1","Scgb3a2","Lamp3","Ager"    #PNEC (肺神经内分泌细胞)
  
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
  feature_plots[[gene]] <- FeaturePlot(ScRNA, features = gene, reduction = "umap", pt.size=3,
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
pdf(paste0(output, "/spacial_FeaturePlot_umap.pdf"), width = 15, height = 12)
print(cowplot::plot_grid(plotlist = feature_plots, ncol = 3))
dev.off()


pdf(paste0(output, "/spacial_VlnPlot_umap.pdf"), width = 12, height = 8)
print(cowplot::plot_grid(plotlist = vln_plots, ncol =3))
dev.off()







###############根据Hes1基因表达定义细胞群###########
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

# 挑选差异细胞展示
output <- paste(outdir,'celltype(Notch2)', sep='/')
dir.create(output)

#file_path <- file.path(outdir, "celltype.rds")
file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

# 根据Cbr2表达量定义exp列
Cbr2_threshold <- 0  # 假设表达阈值为0，可根据实际需求调整
scedata@meta.data$exp <- ifelse(scedata@assays$RNA@data["Notch2", ] > Cbr2_threshold, "Notch2+", "Notch2-")

# 检查定义结果
View(scedata@meta.data)

library(ggsci)
# 绘制细胞类型的umap图
pdf(file = paste0(output, "/ann_umap.pdf"), width = 4, height = 3)
DimPlot(object=scedata,group.by = "exp",reduction='umap',pt.size=1,label=TRUE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_blank())
#scale_color_npg()

dev.off() 

# 绘制细胞类型的umap图（SVG格式）
svg(file = paste0(output, "/ann_umap.svg"), width = 4, height = 3)
DimPlot(object=scedata,group.by = "exp",reduction='umap',pt.size=0.8,label=TRUE,label.size = 5,repel = TRUE,cols=col)+
  theme_dr(xlength = 0.15, 
           ylength = 0.15,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_blank())
#scale_color_npg()
dev.off()

pdf(paste0(output, "/ann-diff-umap.pdf"),width=10,height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  scale_color_npg()
dev.off()

svg(paste0(output, "/ann-diff-umap.svg"),width=10,height=5)
DimPlot(scedata,reduction = "umap",split.by = "treatment",pt.size=,label=FALSE,label.size = 5,repel = TRUE,cols=col)+
  scale_color_npg()
dev.off()


# 保存更新后的数据对象
saveRDS(scedata,"celltype(Notch2).rds")




####### 计算细胞比例 ###########
col <- c("#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8', "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

output <- paste(outdir,'celltype(Notch2)', sep='/')
dir.create(output)

file_path <- file.path(outdir, "celltype(Notch2).rds")
scedata <- readRDS(file_path)

table(scedata$seurat_clusters)


######## 计算所有样本不同细胞群的细胞数
# 计算Cbr2+和Cbr2-的细胞数
cell_counts <- as.data.frame(table(scedata@meta.data$exp))
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


##################### 计算 Uchl1+ / Uchl1- 在不同 treatment 中的细胞数与比例 #####################
# 从meta.data中提取信息
meta <- scedata@meta.data

# 统计 Counts
cell_counts_exp_treat <- meta %>%
  group_by(treatment, exp) %>%
  summarise(Counts = n(), .groups = "drop")

# 计算比例
cell_counts_exp_treat <- cell_counts_exp_treat %>%
  group_by(treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

# 确保 treatment 顺序一致
cell_counts_exp_treat$treatment <- factor(cell_counts_exp_treat$treatment,levels = c("0d", "3d", "7d", "14d"))

################## 绘制 Counts 堆叠图 ##################
p1 <- ggplot(cell_counts_exp_treat, aes(x = treatment, y = Counts, fill = exp)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = '', y = 'Counts') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))

# 保存 Counts 图
ggsave(paste0(output, "/genecount_treatment.pdf"), plot = p1, width = 7, height = 6, dpi = 800)
ggsave(paste0(output, "/genecount_treatment.svg"), plot = p1, width = 7, height = 6, dpi = 800)


################## 绘制 Ratio 堆叠图 ##################
p2 <- ggplot(cell_counts_exp_treat, aes(x = treatment, y = Ratio, fill = exp)) + 
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = '', y = 'Ratio') +
  scale_fill_manual(values = col) +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5),
        axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 20)) +
  geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 7)


# 保存 Ratio 图
ggsave(paste0(output, "/geneRatio_treatment.pdf"), plot = p2, width = 7, height = 6, dpi = 800)
ggsave(paste0(output, "/geneRatio_treatment.svg"), plot = p2, width = 7, height = 6, dpi = 800)




















