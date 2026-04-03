




##可视化
rm(list=ls())


############## GSVA分析 ########
library(ggplot2)
library(dplyr)
library(msigdbr)
library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)

#BiocManager::install("GSVA")

setwd("/fastdata/outdir-1/R/gsva分析/20251222-WF-LC-转移癌-cancer cell/")
outdir <- "/fastdata/outdir-1/R/gsva分析/20251222-WF-LC-转移癌-cancer cell/"


## 输出目录
output <- file.path(outdir, "GSVA")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

## 读取对象
file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)
summary(scRNAsub$Cancer_3class)
view(scRNAsub@meta.data)

#################### 1. 准备表达矩阵与meta ####################
# expr：基因 x 细胞
expr <- as.data.frame(scRNAsub@assays$RNA@data)
expr <- expr[rowSums(expr)>0,]  #选取非零基因

# meta：需要的分组信息，这里假定 meta.data 中已有 celltype / stage / treatment
meta <- scRNAsub@meta.data[, c("Metastatic_3class", "Digestive_tract_3class","Ovarian_3class","Cancer_3class","treatment")]


#  这里的分组会决定接下来的分析哦！
#Idents(scRNAsub) <- scRNAsub$celltype
#expr <- AverageExpression(scRNAsub, assays = "RNA", slot = "data")[[1]]
#expr <- expr[rowSums(expr)>0,]  #选取非零基因
#expr <- as.matrix(expr)
#head(expr)


#################### 2. 提取 KEGG 基因集并做 GSVA ####################
# 先拿到 C2 的所有基因集，再筛选 KEGG 相关的子集合
m_all_c2 <- msigdbr(
  species    = "Homo sapiens",
  collection = "C2"               # 代替原来的 category = "C2"
)

# 保留 KEGG 相关的两个子集合：CP:KEGG_LEGACY 和 CP:KEGG_MEDICUS
m_df <- dplyr::filter(
  m_all_c2,
  gs_subcollection %in% c("CP:KEGG_LEGACY")
)

msigdbr_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)

expr_mat <- as.matrix(expr)
# 注意：如为 UMI-data，严格来说 kcdf="Poisson" 更常见，这里参考示例保持 Gaussian
# 1) 构建参数对象
gsva_par <- gsvaParam(
  exprData = expr_mat,
  geneSets = msigdbr_list,
  kcdf     = "Gaussian"
  # 或 "Poisson" 视你的数据类型而定
  # 其他参数如 minSize / maxSize / tau 等也可以在这里设置
)

# 2) 并行参数（Linux/服务器推荐 MulticoreParam）
bp <- MulticoreParam(workers = 20)

# 3) 运行 GSVA
kegg <- gsva(gsva_par, verbose = TRUE)

rownames(kegg) <- gsub("^KEGG_", "", rownames(kegg))

write.csv(
  kegg,
  file = file.path(output, "GSVA_scores.csv"),
  quote = FALSE
)


pdf(file.path(output, "GSVA_all_cells_heatmap.pdf"), width = 15, height = 12)

pheatmap(kegg,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = meta,
         fontsize_row = 5)

dev.off()








#################### 5. 按不同分组求平均 GSVA、绘制热图（celltype / stage / treatment） ####################

# 封装一个通用函数：按某个分组求平均 GSVA 并画热图
make_group_gsva_heatmap <- function(score_mat,  # kegg, pathway x cells
                                    group_vec,  # 一个分组向量（长度 = 细胞数）
                                    top_n = 30,
                                    filename = "group_gsva.pdf",
                                    main = NULL) {
  # 确保顺序一致
  group_vec <- as.factor(group_vec)
  
  # 按分组求每个通路在该组的平均 score
  data <- score_mat
  # 结果矩阵：通路 x 分组
  group_levels <- levels(group_vec)
  data1 <- NULL
  
  for (lvl in group_levels) {
    idx <- which(group_vec == lvl)
    dat <- rowMeans(data[, idx, drop = FALSE])
    data1 <- cbind(data1, dat)
  }
  
  colnames(data1) <- group_levels
  rownames(data1) <- rownames(score_mat)
  
  # 标准化并截断
  result <- t(scale(t(data1)))
  result[result > 2]  <- 2
  result[result < -2] <- -2
  
  # 选取前 top_n 个通路（这里可以根据方差排序以更有代表性）
  if (nrow(result) > top_n) {
    vars <- apply(result, 1, var)
    sel  <- order(vars, decreasing = TRUE)[1:top_n]
    result_plot <- result[sel, ]
  } else {
    result_plot <- result
  }
  
  p <- pheatmap(result_plot,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                color = colorRampPalette(c('#3399CC', "white","#FF3366"))(100),
                fontsize = 10,
                main = main)
  
  pdf(filename, width = 10, height =6)
  print(p)
  dev.off()
}

## 5.1 按 celltype 分组
meta_ct <- meta[order(meta$celltype), , drop = FALSE]
kegg_ct <- kegg[, rownames(meta_ct)]
group_ct <- meta_ct$celltype

make_group_gsva_heatmap(
  score_mat = kegg_ct,
  group_vec = group_ct,
  top_n = 50,
  filename = file.path(output, "GSVA_celltype.pdf"),
  main = "GSVA by celltype"
)



## 5.2 按 stage 分组
meta_treatment <- meta[order(meta$treatment), , drop = FALSE]
kegg_treatment <- kegg[, rownames(meta_treatment)]
group_treatment <- meta_treatment$treatment

make_group_gsva_heatmap(
  score_mat = kegg_treatment,
  group_vec = group_treatment,
  top_n = 50,
  filename = file.path(output, "GSVA_treatment.pdf"),
  main = ""
)










#################### 4. 将 GSVA 结果加入 Seurat 对象并画 FeaturePlot ####################
es <- data.frame(t(kegg), stringsAsFactors = FALSE)
scRNAsub <- AddMetaData(scRNAsub, es)

# 挑选几个感兴趣的 KEGG 通路名称（请按自己数据替换）
# 可以用: rownames(kegg) 查看所有通路名称
feat1 <- "KEGG_APOPTOSIS"
feat2 <- "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY"
feat3 <- "KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY"
feat4 <- "KEGG_AXON_GUIDANCE"

p1 <- FeaturePlot(scRNAsub, features = feat1, reduction = "umap",cols = c("#0066CC", "#CCFFFF","#FF0066"))
p2 <- FeaturePlot(scRNAsub, features = feat2, reduction = "umap",cols = c("#0066CC", "#CCFFFF","#FF0066"))
p3 <- FeaturePlot(scRNAsub, features = feat3, reduction = "umap",cols = c("#0066CC", "#CCFFFF","#FF0066"))
p4 <- FeaturePlot(scRNAsub, features = feat4, reduction = "umap",cols = c("#0066CC", "#CCFFFF","#FF0066"))

plotc <- (p1 | p2) / (p3 | p4)

ggsave(
  filename = file.path(output, "GSVA_featureplot.pdf"),
  plot = plotc,
  width = 10,
  height = 8
)




#################### 3. 按 treatment + celltype 组合分组求平均 GSVA 得分 ####################

# 修正：应该是 treatment 和 celltype 的组合
meta$combined_group <- paste0(meta$treatment, "_", meta$celltype)
cluster_vec <- meta$combined_group
names(cluster_vec) <- rownames(meta)

# 确保列名与 meta 行名对应
kegg11 <- kegg[, names(cluster_vec), drop = FALSE]

# 计算每个组合组的平均得分
kegg_combined <- sapply(
  X = unique(cluster_vec),
  FUN = function(cl){
    cells <- names(cluster_vec)[cluster_vec == cl]
    rowMeans(kegg11[, cells, drop = FALSE])
  }
)

# 转成矩阵，行为通路，列为组合组
kegg_combined <- as.matrix(kegg_combined)

#################### 4. 重新排列列顺序，确保每个treatment内的celltype顺序一致 ####################

# 提取所有唯一的celltype并排序，确保一致的顺序
all_celltypes <- sort(unique(meta$celltype))

# 提取所有唯一的treatment并排序
all_treatments <- sort(unique(meta$treatment))

# 创建新的列顺序：先按treatment排序，然后在每个treatment内按一致的celltype顺序排序
new_col_order <- c()
for(treat in all_treatments) {
  for(celltype in all_celltypes) {
    group_name <- paste0(treat, "_", celltype)
    if(group_name %in% colnames(kegg_combined)) {
      new_col_order <- c(new_col_order, group_name)
    }
  }
}

# 重新排列矩阵的列
kegg_combined <- kegg_combined[, new_col_order, drop = FALSE]

#################### 5. 选出前 50 个最有变异的通路 ####################
pathway_var <- apply(kegg_combined, 1, var, na.rm = TRUE)
pathway_var <- sort(pathway_var, decreasing = TRUE)

top50_pathways <- names(pathway_var)[1:6]
kegg_top50 <- kegg_combined[top50_pathways, , drop = FALSE]

#################### 6. 准备热图注释信息 ####################
# 提取分组信息用于列注释
col_annot <- data.frame(
  treatment = factor(gsub("_.*", "", colnames(kegg_top50))),
  celltype = factor(gsub(".*_", "", colnames(kegg_top50)))
)
rownames(col_annot) <- colnames(kegg_top50)

# 设置颜色方案
library(RColorBrewer)

# 为treatment生成颜色
treatment_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(col_annot$treatment)))
names(treatment_colors) <- unique(col_annot$treatment)

# 为celltype生成颜色
celltype_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(col_annot$celltype)))
names(celltype_colors) <- sort(unique(col_annot$celltype))

# 创建注释颜色列表
annotation_colors <- list(
  treatment = treatment_colors,
  celltype = celltype_colors
)

#################### 7. 画热图（展示前 50 通路） ####################
library(pheatmap)

# 对数据进行标准化（按行标准化）
kegg_top50_scaled <- t(scale(t(kegg_top50)))

# 确定分组分隔位置 - 在treatment变化的地方添加分隔线
treatment_labels <- col_annot$treatment
gaps_positions <- c()
for(i in 1:(length(treatment_labels)-1)) {
  if(treatment_labels[i] != treatment_labels[i+1]) {
    gaps_positions <- c(gaps_positions, i)
  }
}

# 确保gaps_positions不超过矩阵列数
if(length(gaps_positions) > 0) {
  gaps_positions <- gaps_positions[gaps_positions < ncol(kegg_top50_scaled)]
}

# 强烈建议显式命名空间，避免被其他包的同名函数/翻译机制影响
library(pheatmap)

# 计算画布尺寸（单位：inch）
pdf_w <- max(8, ncol(kegg_top50) * 0.5 + 4)
pdf_h <- max(6, nrow(kegg_top50) * 0.15 + 2)

pdf(file.path(output, "GSVA_KEGG_top50.pdf"), width = pdf_w, height = pdf_h)

p <- pheatmap::pheatmap(
  kegg_top50_scaled,
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  scale          = "none",
  annotation_col = col_annot,
  annotation_colors = annotation_colors,
  show_rownames  = TRUE,
  show_colnames  = FALSE,
  fontsize_row   = 8,
  fontsize_col   = 10,
  color          = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main           = "GSVA KEGG top50 (sorted by treatment then celltype)",
  border_color   = NA,
  cellwidth      = ifelse(ncol(kegg_top50) > 20, 15, 20),
  cellheight     = ifelse(nrow(kegg_top50) > 40, 8, 10),
  gaps_col       = if(length(gaps_positions) > 0) gaps_positions else NULL
)

dev.off()



# 可选：保存分组信息到文件
write.csv(col_annot, file.path(output, "heatmap_group_annotations.csv"), row.names = TRUE)
write.csv(kegg_top50, file.path(output, "top50_pathways_scores.csv"), row.names = TRUE)





#################### 3. 按 treatment + celltype 组合分组求平均 GSVA 得分 ####################
# 创建组合分组
meta$combined_group <- paste0(meta$treatment, "_", meta$celltype)
cluster_vec <- meta$combined_group
names(cluster_vec) <- rownames(meta)

# 确保列名与 meta 行名对应
kegg11 <- kegg[, names(cluster_vec), drop = FALSE]

# 计算每个组合组的平均得分
kegg_combined <- sapply(
  X = unique(cluster_vec),
  FUN = function(cl){
    cells <- names(cluster_vec)[cluster_vec == cl]
    rowMeans(kegg11[, cells, drop = FALSE])
  }
)

# 转成矩阵，行为通路，列为组合组
kegg_combined <- as.matrix(kegg_combined)

#################### 4. 重新排列列顺序，按celltype分组，每个celltype内按treatment排序 ####################
# 提取所有唯一的celltype和treatment
all_celltypes <- sort(unique(meta$celltype))
all_treatments <- sort(unique(meta$treatment))

# 创建新的列顺序：先按celltype排序，然后在每个celltype内按treatment排序
new_col_order <- c()
for(celltype in all_celltypes) {
  for(treat in all_treatments) {
    group_name <- paste0(treat, "_", celltype)
    if(group_name %in% colnames(kegg_combined)) {
      new_col_order <- c(new_col_order, group_name)
    }
  }
}

# 重新排列矩阵的列
kegg_combined <- kegg_combined[, new_col_order, drop = FALSE]

#################### 5. 选出前 50 个最有变异的通路 ####################
pathway_var <- apply(kegg_combined, 1, var, na.rm = TRUE)
pathway_var <- sort(pathway_var, decreasing = TRUE)

top50_pathways <- names(pathway_var)[1:6]
kegg_top50 <- kegg_combined[top50_pathways, , drop = FALSE]



# 定义感兴趣的通路列表，并删除"KEGG_"前缀
kegg_interest_ids <- c(
  "JAK_STAT_SIGNALING_PATHWAY",
  "PI3K_AKT_SIGNALING_PATHWAY", "ABC_TRANSPORTERS", "GLUTATHIONE_METABOLISM", 
  "METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450", "DRUG_METABOLISM_OTHER_ENZYMES", 
  "PATHWAYS_IN_CANCER", "MAPK_SIGNALING_PATHWAY", 
  "APOPTOSIS", "MTOR_SIGNALING_PATHWAY", "FOCAL_ADHESION", 
  "TGF_BETA_SIGNALING_PATHWAY", "HEDGEHOG_SIGNALING_PATHWAY", "WNT_SIGNALING_PATHWAY", 
  "CELL_ADHESION_MOLECULES_CAMS", "GLYCOLYSIS_GLUCONEOGENESIS","TIGHT_JUNCTION","OXIDATIVE_PHOSPHORYLATION"
  
)

# 找出实际存在于矩阵中的通路
existing_pathways <- kegg_interest_ids[kegg_interest_ids %in% rownames(kegg_combined)]

# 输出未找到的通路
kegg_interest <- kegg_combined[existing_pathways, , drop = FALSE]






#################### 6. 准备热图注释信息 ####################
# 提取分组信息用于列注释
col_annot <- data.frame(
  treatment = factor(gsub("_.*", "", colnames(kegg_top50)), levels = all_treatments),
  celltype = factor(gsub(".*_", "", colnames(kegg_top50)), levels = all_celltypes)
)
rownames(col_annot) <- colnames(kegg_top50)

# 设置颜色方案
library(RColorBrewer)

# 为treatment生成颜色 - 使用Set2调色板
treatment_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(all_treatments))
names(treatment_colors) <- all_treatments

# 为celltype生成颜色 - 使用Set3调色板（提供更多颜色）
celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_celltypes))
names(celltype_colors) <- all_celltypes

# 创建注释颜色列表
annotation_colors <- list(
  treatment = treatment_colors,
  celltype = celltype_colors
)

#################### 7. 计算每个celltype内treatment间的差异 ####################
# 可选：计算每个celltype内部不同treatment间的差异矩阵
# 这有助于识别哪些通路在特定celltype中对treatment变化最敏感

celltype_treatment_differences <- list()

for(celltype in all_celltypes) {
  # 提取当前celltype的所有列
  celltype_cols <- grep(paste0("_", celltype, "$"), colnames(kegg_interest), value = TRUE)
  
  if(length(celltype_cols) > 1) {
    # 计算当前celltype内部treatment间的平均差异
    celltype_data <- kegg_top50[, celltype_cols, drop = FALSE]
    
    # 如果有多个treatment，计算两两之间的平均绝对差异
    if(ncol(celltype_data) > 1) {
      # 计算所有通路在所有treatment组合中的平均绝对差异
      diff_matrix <- matrix(0, nrow = nrow(celltype_data), ncol = ncol(celltype_data),
                            dimnames = list(rownames(celltype_data), colnames(celltype_data)))
      
      for(i in 1:ncol(celltype_data)) {
        for(j in 1:ncol(celltype_data)) {
          if(i != j) {
            diff_matrix[, i] <- diff_matrix[, i] + abs(celltype_data[, i] - celltype_data[, j])
          }
        }
        diff_matrix[, i] <- diff_matrix[, i] / (ncol(celltype_data) - 1)
      }
      
      # 计算每个通路的平均差异
      pathway_avg_diff <- rowMeans(diff_matrix)
      pathway_avg_diff <- sort(pathway_avg_diff, decreasing = TRUE)
      
      celltype_treatment_differences[[celltype]] <- list(
        data = celltype_data,
        diff_matrix = diff_matrix,
        pathway_avg_diff = pathway_avg_diff
      )
    }
  }
}


#################### 7. 画热图（展示前 50 通路） ####################
library(ComplexHeatmap)
library(circlize)

# 对数据进行标准化（按行标准化）
kegg_top50_scaled <- t(scale(t(kegg_top50)))

# 确定分组分隔位置 - 在celltype变化的地方添加分隔线
celltype_labels <- col_annot$celltype
gaps_positions <- c()
for(i in 1:(length(celltype_labels)-1)) {
  if(celltype_labels[i] != celltype_labels[i+1]) {
    gaps_positions <- c(gaps_positions, i)
  }
}

# 创建列注释
col_ha <- HeatmapAnnotation(
  df = col_annot,
  col = annotation_colors,
  annotation_name_side = "left",
  gap = unit(2, "mm")
)

# 创建行注释（可选，显示通路分类）
# row_ha <- rowAnnotation(
#   pathway = rownames(kegg_top50_scaled),
#   show_annotation_name = FALSE
# )

# 绘制热图
ht <- Heatmap(
  kegg_top50_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  
  # 行设置
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 12),
  row_names_max_width = unit(10, "cm"),
  row_title = "KEGG Pathways",
  
  # 列设置
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = col_annot$celltype,
  column_gap = unit(3, "mm"),
  column_title = "Cell Types",
  column_title_side = "top",
  
  # 注释
  top_annotation = col_ha,
  
  # 热图参数
  border = TRUE,
  border_gp = gpar(col = "gray", lty = 1),
  heatmap_legend_param = list(
    title = "Z-score",
    title_position = "leftcenter-rot"
  ),
  
  # 行聚类参数
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  # 显示参数
  show_heatmap_legend = TRUE,
  use_raster = ifelse(ncol(kegg_top50_scaled) * nrow(kegg_top50_scaled) > 20000, TRUE, FALSE)
)

# 绘制并保存
pdf(file.path(output, "GSVA_KEGG_top50_by_celltype_treatment.pdf"), 
    width = max(10, ncol(kegg_top50) * 0.4 + 4),
    height = max(8, nrow(kegg_top50) * 0.15 + 3))

draw(ht, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)

dev.off()


#################### 7. 画热图（展示前 50 通路）- 优化版本 ####################
library(ComplexHeatmap)
library(circlize)

# 对数据进行标准化（按行标准化）
kegg_top50_scaled <- t(scale(t(kegg_top50)))

# 创建颜色映射函数
col_fun <- colorRamp2(
  seq(min(kegg_top50_scaled, na.rm = TRUE), 
      max(kegg_top50_scaled, na.rm = TRUE), 
      length = 100),
  colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)

# 创建列注释 - 优化图例参数
col_ha <- HeatmapAnnotation(
  df = col_annot,
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_legend_param = list(
    # 设置图例为两行显示
    nrow = 2,
    ncol = ceiling(length(annotation_colors) / 2),
    direction = "horizontal",
    
    # 增大图例标题大小
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    # 增大图例标签大小
    labels_gp = gpar(fontsize = 12),
    
    # 增大图例图形大小
    legend_height = unit(1.2, "cm"),
    legend_width = unit(1.2, "cm"),
    
    # 设置图例之间的间距
    gap = unit(0.3, "cm"),
    
    title_position = "topcenter"
  ),
  gap = unit(2, "mm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12)
)

# 绘制热图
ht <- Heatmap(
  kegg_top50_scaled,
  name = "Z-score",
  col = col_fun,
  
  # 行设置
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_names_max_width = unit(14, "cm"),
  row_title = "KEGG Pathways",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # 列设置
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = col_annot$celltype,
  column_gap = unit(3, "mm"),
  column_title = "Cell Types",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # 注释
  top_annotation = col_ha,
  
  # 热图参数
  border = TRUE,
  border_gp = gpar(col = "gray", lty = 1),
  
  # 热图图例设置 - 优化版本
  heatmap_legend_param = list(
    title = "GSVA Z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(10, "cm"),
    legend_height = unit(1.5, "cm"),
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    at = round(seq(min(kegg_top50_scaled), max(kegg_top50_scaled), length.out = 5), 1),
    grid_border = "black",
    color_bar = "continuous",
    labels_rot = 0,
    border = "black",
    padding = unit(c(8, 8, 8, 8), "mm")
  ),
  
  # 行聚类参数
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  # 显示参数
  show_heatmap_legend = TRUE,
  use_raster = ifelse(ncol(kegg_top50_scaled) * nrow(kegg_top50_scaled) > 20000, TRUE, FALSE),
  
  # 整体尺寸设置 - 修改了height参数，增大了热图高度
  width = unit(max(15, ncol(kegg_top50_scaled) * 0.4 + 4), "cm"),
  height = unit(max(4, nrow(kegg_top50_scaled) * 0.1 + 4), "cm")  # 从0.15增加到0.25，基础值从3增加到4
)

# 绘制并保存 - 修改了PDF高度设置
pdf(file.path(output, "GSVA_KEGG_top50_by_celltype_treatment.pdf"), 
    width = max(10, ncol(kegg_top50) * 0.3 + 6),
    height = max(6, nrow(kegg_top50) * 0.1 + 5))  # 从0.12增加到0.2，基础值从4增加到5

# 绘制热图，图例在底部，分两行显示
draw(ht, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legends = TRUE,
     ht_gap = unit(2, "cm"),
     # 增大底部边距以容纳两行图例
     padding = unit(c(20, 20, 60, 20), "mm")  # 上、右、下、左边距
)

dev.off()




#################### 7. 只选择 Treg / Tex / CTL 三种 celltype ####################


# 需要保留的 celltype
selected_celltypes <- c("PNEC")

# 基于已有的 col_annot 来筛选列（样本）
keep_cols <- col_annot$celltype %in% selected_celltypes

# 筛选后的 KEGG GSVA 矩阵
kegg_interest_sel <- kegg_top50[, keep_cols, drop = FALSE]

# 筛选后的列注释
col_annot_sel <- col_annot[keep_cols, , drop = FALSE]
col_annot_sel$celltype <- droplevels(col_annot_sel$celltype)
col_annot_sel$treatment <- droplevels(col_annot_sel$treatment)

# 颜色也只保留这三种 celltype 对应的颜色
annotation_colors_sel <- annotation_colors
annotation_colors_sel$celltype <- annotation_colors_sel$celltype[names(annotation_colors_sel$celltype) %in% selected_celltypes]


#################### 8. 画 Treg / Tex / CTL 的 KEGG 热图（Top 50，优化版本） ####################
library(ComplexHeatmap)
library(circlize)

# 对数据进行标准化（按行标准化）
kegg_top50_scaled_sel <- t(scale(t(kegg_interest_sel)))

# 创建颜色映射函数
col_fun <- colorRamp2(
  seq(min(kegg_top50_scaled_sel, na.rm = TRUE), 
      max(kegg_top50_scaled_sel, na.rm = TRUE), 
      length = 100),
  colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)

# 创建列注释 - 优化图例参数
col_ha_sel <- HeatmapAnnotation(
  df = col_annot_sel,
  col = annotation_colors_sel,
  annotation_name_side = "left",
  annotation_legend_param = list(
    nrow = 2,
    direction = "horizontal",
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    legend_height = unit(1.2, "cm"),
    legend_width = unit(1.2, "cm"),
    gap = unit(0.3, "cm"),
    title_position = "topcenter"
  ),
  gap = unit(2, "mm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12)
)

# 热图对象
ht_sel <- Heatmap(
  kegg_top50_scaled_sel,
  name = "Z-score",
  col = col_fun,
  
  # 行设置
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_names_max_width = unit(14, "cm"),
  row_title = "KEGG Pathways",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # 列设置
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = col_annot_sel$celltype,   # 这里只会有 Treg / Tex / CTL
  column_gap = unit(3, "mm"),
  column_title = "Cell Types",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # 注释
  top_annotation = col_ha_sel,
  
  # 热图边框
  border = TRUE,
  border_gp = gpar(col = "gray", lty = 1),
  
  # 热图图例设置
  heatmap_legend_param = list(
    title = "GSVA Z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(10, "cm"),
    legend_height = unit(1.5, "cm"),
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    at = round(seq(min(kegg_top50_scaled_sel, na.rm = TRUE), 
                   max(kegg_top50_scaled_sel, na.rm = TRUE), 
                   length.out = 5), 1),
    grid_border = "black",
    color_bar = "continuous",
    labels_rot = 0,
    border = "black",
    padding = unit(c(8, 8, 8, 8), "mm")
  ),
  
  # 行聚类参数
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  # 显示参数
  show_heatmap_legend = TRUE,
  use_raster = ifelse(ncol(kegg_top50_scaled_sel) * nrow(kegg_top50_scaled_sel) > 20000, TRUE, FALSE),
  
  # 整体尺寸设置
  width = unit(max(6, ncol(kegg_top50_scaled_sel) * 0.4 + 4), "cm"),
  height = unit(max(6, nrow(kegg_top50_scaled_sel) * 0.3 + 4), "cm")
)

# 保存 PDF （建议单独命名，避免覆盖全部 celltype 的结果）
pdf(file.path(output, "GSVA_KEGG_top50_PNEC.pdf"), 
    width = max(6, ncol(kegg_top50_scaled_sel) * 0.3 + 6),
    height = max(5, nrow(kegg_top50_scaled_sel) * 0.2 + 5))

draw(ht_sel, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legends = TRUE,
     ht_gap = unit(2, "cm"),
     padding = unit(c(20, 20, 60, 20), "mm")  # 上、右、下、左边距
)

dev.off()









#################### 设置分组变量列表 ####################
group_vars <- c("stage", "treatment", "Type", "BRCA1/2", "P53", "PD-L1", "HER2")

# 检查这些变量是否存在于meta数据中
existing_vars <- group_vars[group_vars %in% colnames(meta)]
missing_vars <- setdiff(group_vars, existing_vars)

if (length(missing_vars) > 0) {
  warning(paste("以下变量在meta数据中不存在:", paste(missing_vars, collapse = ", ")))
}

# 只使用存在的变量
group_vars <- existing_vars

if (length(group_vars) == 0) {
  stop("没有找到任何指定的分组变量！")
}

cat("将处理以下分组变量:", paste(group_vars, collapse = ", "), "\n")

#################### 定义自定义通路列表 ####################
kegg_interest_ids <- c(
  "JAK_STAT_SIGNALING_PATHWAY",
  "PI3K_AKT_SIGNALING_PATHWAY", "ABC_TRANSPORTERS", "GLUTATHIONE_METABOLISM", 
  "METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450", "DRUG_METABOLISM_OTHER_ENZYMES", 
  "PATHWAYS_IN_CANCER", "MAPK_SIGNALING_PATHWAY", 
  "APOPTOSIS", "MTOR_SIGNALING_PATHWAY", "FOCAL_ADHESION", 
  "TGF_BETA_SIGNALING_PATHWAY", "HEDGEHOG_SIGNALING_PATHWAY", "WNT_SIGNALING_PATHWAY", 
  "CELL_ADHESION_MOLECULES_CAMS", "GLYCOLYSIS_GLUCONEOGENESIS","TIGHT_JUNCTION","OXIDATIVE_PHOSPHORYLATION",
  
  "T_CELL_RECEPTOR_SIGNALING_PATHWAY", "CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "TRYPTOPHAN_METABOLISM", "ARGININE_AND_PROLINE_METABOLISM", 
  "PURINE_METABOLISM", "PYRIMIDINE_METABOLISM",
  "TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","PPAR_SIGNALING_PATHWAY",
  "NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY", "ANTIGEN_PROCESSING_AND_PRESENTATION"
  
)

#################### 循环处理每个分组变量 ####################
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 创建结果目录
output_dir <- file.path(output, "celltype_comparisons")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 存储所有结果
all_results <- list()

for (group_var in group_vars) {
  cat("\n正在处理分组变量:", group_var, "\n")
  
  # 检查分组变量是否有足够的非NA值
  if (sum(!is.na(meta[[group_var]])) == 0) {
    warning(paste("分组变量", group_var, "所有值都是NA，跳过..."))
    next
  }
  
  # 去除分组变量为NA的样本
  valid_samples <- !is.na(meta[[group_var]])
  meta_sub <- meta[valid_samples, ]
  kegg_sub <- kegg[, rownames(meta_sub), drop = FALSE]
  
  # 确保至少有两个组别
  group_levels <- unique(meta_sub[[group_var]])
  group_levels <- group_levels[!is.na(group_levels)]
  
  if (length(group_levels) < 2) {
    warning(paste("分组变量", group_var, "只有", length(group_levels), "个水平，跳过..."))
    next
  }
  
  cat("分组水平:", paste(group_levels, collapse = ", "), "\n")
  
  #################### 1. 创建组合分组 ####################
  meta_sub$combined_group <- paste0(meta_sub$celltype, "_", meta_sub[[group_var]])
  cluster_vec <- meta_sub$combined_group
  names(cluster_vec) <- rownames(meta_sub)
  
  # 确保列名与meta行名对应
  kegg_sub_filtered <- kegg_sub[, names(cluster_vec), drop = FALSE]
  
  #################### 2. 计算每个组合组的平均得分 ####################
  combined_groups <- unique(cluster_vec)
  cat("共创建", length(combined_groups), "个组合组\n")
  
  if (length(combined_groups) < 2) {
    warning(paste("分组变量", group_var, "只产生了", length(combined_groups), "个组合组，跳过..."))
    next
  }
  
  # 计算平均得分
  kegg_combined <- sapply(
    X = combined_groups,
    FUN = function(cl){
      cells <- names(cluster_vec)[cluster_vec == cl]
      rowMeans(kegg_sub_filtered[, cells, drop = FALSE])
    }
  )
  
  kegg_combined <- as.matrix(kegg_combined)
  
  #################### 3. 重新排列列顺序 ####################
  # 按celltype分组，每个celltype内按分组变量排序
  all_celltypes <- sort(unique(meta_sub$celltype))
  all_group_levels <- sort(unique(meta_sub[[group_var]]))
  
  # 创建新的列顺序
  new_col_order <- c()
  for(celltype in all_celltypes) {
    for(group_val in all_group_levels) {
      group_name <- paste0(celltype, "_", group_val)
      if(group_name %in% colnames(kegg_combined)) {
        new_col_order <- c(new_col_order, group_name)
      }
    }
  }
  
  # 重新排列矩阵的列
  if(length(new_col_order) > 0) {
    kegg_combined <- kegg_combined[, new_col_order, drop = FALSE]
  }
  
  #################### 4. 提取感兴趣的通路 ####################
  # 找出实际存在于矩阵中的通路
  existing_pathways <- kegg_interest_ids[kegg_interest_ids %in% rownames(kegg_combined)]
  cat("找到", length(existing_pathways), "个感兴趣的通路\n")
  
  if(length(existing_pathways) < 5) {
    warning(paste("分组变量", group_var, "只有", length(existing_pathways), "个感兴趣的通路，跳过..."))
    next
  }
  
  kegg_interest <- kegg_combined[existing_pathways, , drop = FALSE]
  
  #################### 5. 准备热图注释信息 ####################
  # 提取分组信息用于列注释
  col_annot <- data.frame(
    celltype = factor(gsub("_.*", "", colnames(kegg_interest)), levels = all_celltypes),
    group_var = factor(gsub(".*_", "", colnames(kegg_interest)), levels = all_group_levels)
  )
  colnames(col_annot)[2] <- group_var  # 使用实际的分组变量名
  rownames(col_annot) <- colnames(kegg_interest)
  
  #################### 6. 设置颜色方案 ####################
  # 为celltype生成颜色
  n_celltypes <- length(all_celltypes)
  if(n_celltypes <= 8) {
    celltype_colors <- brewer.pal(max(3, n_celltypes), "Set1")[1:n_celltypes]
  } else if(n_celltypes <= 12) {
    celltype_colors <- brewer.pal(n_celltypes, "Set3")
  } else {
    celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_celltypes)
  }
  names(celltype_colors) <- all_celltypes
  
  # 为分组变量生成颜色
  n_groups <- length(all_group_levels)
  if(n_groups <= 8) {
    group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
  } else {
    group_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_groups)
  }
  names(group_colors) <- all_group_levels
  
  # 创建注释颜色列表
  annotation_colors <- list(
    celltype = celltype_colors
  )
  annotation_colors[[group_var]] <- group_colors
  
  #################### 7. 绘制热图 ####################
  # 对数据进行标准化（按行标准化）
  kegg_interest_scaled <- t(scale(t(kegg_interest)))
  
  # 设置颜色映射
  # 计算合理的颜色断点
  data_range <- range(kegg_interest_scaled, na.rm = TRUE)
  breaks <- seq(data_range[1], data_range[2], length.out = 100)
  col_fun <- colorRamp2(breaks, colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  
  # 创建列注释
  col_ha <- HeatmapAnnotation(
    df = col_annot,
    col = annotation_colors,
    annotation_name_side = "left",
    annotation_legend_param = list(
      nrow = 2,
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    ),
    gap = unit(2, "mm")
  )
  
  # 绘制热图
  ht <- Heatmap(
    kegg_interest_scaled,
    name = "Z-score",
    col = col_fun,
    
    # 行设置
    cluster_rows = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 9),
    row_names_max_width = unit(12, "cm"),
    row_title = paste("KEGG Pathways"),
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    # 列设置
    cluster_columns = FALSE,
    show_column_names = FALSE,
    column_split = col_annot$celltype,
    column_gap = unit(2, "mm"),
    column_title = "Cell Types",
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    
    # 注释
    top_annotation = col_ha,
    
    # 热图参数
    border = TRUE,
    border_gp = gpar(col = "gray", lty = 1),
    
    # 热图图例设置
    heatmap_legend_param = list(
      title = "GSVA Z-score",
      title_position = "topcenter",
      legend_direction = "horizontal",
      legend_width = unit(8, "cm"),
      title_gp = gpar(fontsize = 11),
      labels_gp = gpar(fontsize = 10),
      at = round(seq(data_range[1], data_range[2], length.out = 5), 1)
    ),
    
    # 行聚类参数
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    
    # 显示参数
    show_heatmap_legend = TRUE,
    use_raster = ifelse(ncol(kegg_interest_scaled) * nrow(kegg_interest_scaled) > 10000, TRUE, FALSE)
  )
  
  #################### 8. 保存热图 ####################
  # 计算合适的图形尺寸
  plot_width <- max(10, ncol(kegg_interest_scaled) * 0.3 + 4)
  plot_height <- max(8, nrow(kegg_interest_scaled) * 0.2 + 3)
  
  # 安全地处理文件名中的特殊字符
  # 将斜杠和其他特殊字符替换为下划线
  safe_group_var <- gsub("[\\/\\\\:*?\"<>|]", "_", group_var)
  
  # 创建安全的文件名
  pdf_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", safe_group_var, ".pdf"))
  
  # 检查路径是否有效
  cat("保存文件到:", pdf_file, "\n")
  
  pdf(pdf_file, width = plot_width, height = plot_height)
  
  draw(ht, 
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom",
       merge_legends = TRUE,
       padding = unit(c(20, 20, 40, 20), "mm"))  # 上、右、下、左边距
  
  dev.off()
  
  cat("已保存热图:", pdf_file, "\n")
  
  #################### 9. 保存数据 ####################
  # 保存结果数据
  result_data <- list(
    matrix_raw = kegg_interest,
    matrix_scaled = kegg_interest_scaled,
    annotation = col_annot,
    colors = annotation_colors
  )
  
  all_results[[group_var]] <- result_data
  
  # 保存数据到文件 - 使用安全的文件名
  data_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", safe_group_var, "_data.rds"))
  saveRDS(result_data, file = data_file)
  
  # 保存CSV格式的数据
  csv_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", safe_group_var, "_scores.csv"))
  write.csv(kegg_interest, file = csv_file)
  
  # 保存注释信息
  annot_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", safe_group_var, "_annotation.csv"))
  write.csv(col_annot, file = annot_file)
}









#################### 设置分组变量列表 ####################
group_vars <- c("stage", "treatment", "Type", "BRCA1/2", "P53", "PD-L1", "HER2")

# 检查这些变量是否存在于meta数据中
existing_vars <- group_vars[group_vars %in% colnames(meta)]
missing_vars <- setdiff(group_vars, existing_vars)

if (length(missing_vars) > 0) {
  warning(paste("以下变量在meta数据中不存在:", paste(missing_vars, collapse = ", ")))
}

# 只使用存在的变量
group_vars <- existing_vars

if (length(group_vars) == 0) {
  stop("没有找到任何指定的分组变量！")
}

cat("将处理以下分组变量:", paste(group_vars, collapse = ", "), "\n")

#################### 定义自定义通路列表 ####################
kegg_interest_ids <- c(
  
  "JAK_STAT_SIGNALING_PATHWAY",
  "PI3K_AKT_SIGNALING_PATHWAY", "ABC_TRANSPORTERS", "GLUTATHIONE_METABOLISM", 
  "METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450", "DRUG_METABOLISM_OTHER_ENZYMES", 
  "PATHWAYS_IN_CANCER", "MAPK_SIGNALING_PATHWAY", 
  "APOPTOSIS", "MTOR_SIGNALING_PATHWAY", "FOCAL_ADHESION", 
  "TGF_BETA_SIGNALING_PATHWAY", "HEDGEHOG_SIGNALING_PATHWAY", "WNT_SIGNALING_PATHWAY", 
  "CELL_ADHESION_MOLECULES_CAMS", "GLYCOLYSIS_GLUCONEOGENESIS","TIGHT_JUNCTION","OXIDATIVE_PHOSPHORYLATION",
  
  "T_CELL_RECEPTOR_SIGNALING_PATHWAY", "CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "TRYPTOPHAN_METABOLISM", "ARGININE_AND_PROLINE_METABOLISM", 
  "PURINE_METABOLISM", "PYRIMIDINE_METABOLISM",
  "TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","PPAR_SIGNALING_PATHWAY",
  "NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY", "ANTIGEN_PROCESSING_AND_PRESENTATION"
  
)

#################### 循环处理每个分组变量（只画感兴趣 celltype） ####################
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# 创建结果目录
output_dir <- file.path(output, "selected_celltype_comparisons")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 存储所有结果
all_results <- list()

# 感兴趣的 celltype（CAF 亚群）
selected_celltypes <- c("Treg", "Tex", "CTL")

for (group_var in group_vars) {
  cat("\n正在处理分组变量:", group_var, "\n")
  
  # 检查分组变量是否有足够的非NA值
  if (sum(!is.na(meta[[group_var]])) == 0) {
    warning(paste("分组变量", group_var, "所有值都是NA，跳过..."))
    next
  }
  
  # 去除分组变量为NA的样本
  valid_samples <- !is.na(meta[[group_var]])
  meta_sub <- meta[valid_samples, ]
  kegg_sub <- kegg[, rownames(meta_sub), drop = FALSE]
  
  # 确保至少有两个组别
  group_levels <- unique(meta_sub[[group_var]])
  group_levels <- group_levels[!is.na(group_levels)]
  
  if (length(group_levels) < 2) {
    warning(paste("分组变量", group_var, "只有", length(group_levels), "个水平，跳过..."))
    next
  }
  
  cat("分组水平:", paste(group_levels, collapse = ", "), "\n")
  
  #################### 1. 创建组合分组 ####################
  meta_sub$combined_group <- paste0(meta_sub$celltype, "_", meta_sub[[group_var]])
  cluster_vec <- meta_sub$combined_group
  names(cluster_vec) <- rownames(meta_sub)
  
  # 确保列名与meta行名对应
  kegg_sub_filtered <- kegg_sub[, names(cluster_vec), drop = FALSE]
  
  #################### 2. 计算每个组合组的平均得分 ####################
  combined_groups <- unique(cluster_vec)
  cat("共创建", length(combined_groups), "个组合组\n")
  
  if (length(combined_groups) < 2) {
    warning(paste("分组变量", group_var, "只产生了", length(combined_groups), "个组合组，跳过..."))
    next
  }
  
  # 计算平均得分
  kegg_combined <- sapply(
    X = combined_groups,
    FUN = function(cl){
      cells <- names(cluster_vec)[cluster_vec == cl]
      rowMeans(kegg_sub_filtered[, cells, drop = FALSE])
    }
  )
  
  kegg_combined <- as.matrix(kegg_combined)
  
  #################### 3. 重新排列列顺序 ####################
  # 按celltype分组，每个celltype内按分组变量排序
  all_celltypes <- sort(unique(meta_sub$celltype))
  all_group_levels <- sort(unique(meta_sub[[group_var]]))
  
  # 创建新的列顺序
  new_col_order <- c()
  for(celltype in all_celltypes) {
    for(group_val in all_group_levels) {
      group_name <- paste0(celltype, "_", group_val)
      if(group_name %in% colnames(kegg_combined)) {
        new_col_order <- c(new_col_order, group_name)
      }
    }
  }
  
  # 重新排列矩阵的列
  if(length(new_col_order) > 0) {
    kegg_combined <- kegg_combined[, new_col_order, drop = FALSE]
  }
  
  #################### 4. 提取感兴趣的通路 ####################
  existing_pathways <- kegg_interest_ids[kegg_interest_ids %in% rownames(kegg_combined)]
  cat("找到", length(existing_pathways), "个感兴趣的通路\n")
  
  if(length(existing_pathways) < 5) {
    warning(paste("分组变量", group_var, "只有", length(existing_pathways), "个感兴趣的通路，跳过..."))
    next
  }
  
  kegg_interest <- kegg_combined[existing_pathways, , drop = FALSE]
  
  #################### 5. 准备热图注释信息（所有 celltype） ####################
  col_annot <- data.frame(
    celltype = factor(gsub("_.*", "", colnames(kegg_interest)), levels = all_celltypes),
    group_var = factor(gsub(".*_", "", colnames(kegg_interest)), levels = all_group_levels)
  )
  colnames(col_annot)[2] <- group_var  # 第二列命名为实际分组变量名
  rownames(col_annot) <- colnames(kegg_interest)
  
  #################### 6. 只保留感兴趣的 celltype ####################
  keep_cols <- col_annot$celltype %in% selected_celltypes
  
  if (sum(keep_cols) == 0) {
    warning(paste("分组变量", group_var, "在感兴趣的 celltype 中没有样本，跳过..."))
    next
  }
  
  # 过滤后的 KEGG 矩阵和注释
  kegg_interest_sel <- kegg_interest[, keep_cols, drop = FALSE]
  col_annot_sel <- col_annot[keep_cols, , drop = FALSE]
  
  # 去掉没出现过的水平
  col_annot_sel$celltype <- droplevels(col_annot_sel$celltype)
  col_annot_sel[[group_var]] <- droplevels(col_annot_sel[[group_var]])
  
  # 更新 celltype 和 group 水平
  sel_celltypes_levels <- levels(col_annot_sel$celltype)
  sel_group_levels <- levels(col_annot_sel[[group_var]])
  
  #################### 7. 设置颜色方案（只对筛选后的 celltype） ####################
  # 为所有 celltype 先生成颜色
  n_celltypes <- length(all_celltypes)
  if(n_celltypes <= 8) {
    celltype_colors <- brewer.pal(max(3, n_celltypes), "Set1")[1:n_celltypes]
  } else if(n_celltypes <= 12) {
    celltype_colors <- brewer.pal(n_celltypes, "Set3")
  } else {
    celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_celltypes)
  }
  names(celltype_colors) <- all_celltypes
  
  # 为分组变量生成颜色
  n_groups <- length(all_group_levels)
  if(n_groups <= 8) {
    group_colors <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
  } else {
    group_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_groups)
  }
  names(group_colors) <- all_group_levels
  
  # 只保留感兴趣 celltype 和实际出现的 group 水平对应的颜色
  annotation_colors_sel <- list(
    celltype = celltype_colors[names(celltype_colors) %in% sel_celltypes_levels]
  )
  annotation_colors_sel[[group_var]] <- group_colors[names(group_colors) %in% sel_group_levels]
  
  #################### 8. 标准化数据（行标准化） ####################
  kegg_interest_scaled_sel <- t(scale(t(kegg_interest_sel)))
  
  # 颜色映射函数
  data_min <- min(kegg_interest_scaled_sel, na.rm = TRUE)
  data_max <- max(kegg_interest_scaled_sel, na.rm = TRUE)
  col_fun <- colorRamp2(
    seq(data_min, data_max, length = 100),
    colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  )
  
  #################### 9. 列注释（优化版本） ####################
  col_ha_sel <- HeatmapAnnotation(
    df = col_annot_sel,
    col = annotation_colors_sel,
    annotation_name_side = "left",
    annotation_legend_param = list(
      nrow = 2,
      direction = "horizontal",
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 12),
      legend_height = unit(1.2, "cm"),
      legend_width = unit(1.2, "cm"),
      gap = unit(0.3, "cm"),
      title_position = "topcenter"
    ),
    gap = unit(2, "mm"),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 12)
  )
  
  #################### 10. 构建热图对象（只含 FAP+ THY1+ CAFs / iCAFs / apCAFs） ####################
  ht_sel <- Heatmap(
    kegg_interest_scaled_sel,
    name = "Z-score",
    col = col_fun,
    
    # 行设置
    cluster_rows = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    row_names_max_width = unit(14, "cm"),
    row_title = "KEGG Pathways",
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    # 列设置
    cluster_columns = FALSE,
    show_column_names = FALSE,
    column_split = col_annot_sel$celltype,   # 这里只会有 CAF 亚群
    column_gap = unit(3, "mm"),
    column_title = "Cell Types",
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    # 注释
    top_annotation = col_ha_sel,
    
    # 热图边框
    border = TRUE,
    border_gp = gpar(col = "gray", lty = 1),
    
    # 热图图例设置（参考你的优化版本）
    heatmap_legend_param = list(
      title = "GSVA Z-score",
      title_position = "topcenter",
      legend_direction = "horizontal",
      legend_width = unit(10, "cm"),
      legend_height = unit(1.5, "cm"),
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 12),
      at = round(seq(data_min, data_max, length.out = 5), 1),
      grid_border = "black",
      color_bar = "continuous",
      labels_rot = 0,
      border = "black",
      padding = unit(c(8, 8, 8, 8), "mm")
    ),
    
    # 行聚类参数
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    
    # 显示参数
    show_heatmap_legend = TRUE,
    use_raster = ifelse(ncol(kegg_interest_scaled_sel) * nrow(kegg_interest_scaled_sel) > 20000, TRUE, FALSE),
    
    # 整体尺寸设置
    width = unit(max(6, ncol(kegg_interest_scaled_sel) * 0.4 + 4), "cm"),
    height = unit(max(6, nrow(kegg_interest_scaled_sel) * 0.3 + 4), "cm")
  )
  
  #################### 11. 保存 PDF/数据（文件名包含分组变量 + selectedCelltypes） ####################
  # 安全文件名
  safe_group_var <- gsub("[\\/\\\\:*?\"<>|]", "_", group_var)
  tag <- "selectedCelltypes"
  
  pdf_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", tag, "_", safe_group_var, ".pdf"))
  cat("保存文件到:", pdf_file, "\n")
  
  pdf(pdf_file,
      width  = max(6, ncol(kegg_interest_scaled_sel) * 0.3 + 6),
      height = max(5, nrow(kegg_interest_scaled_sel) * 0.2 + 5))
  
  draw(
    ht_sel,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends = TRUE,
    ht_gap = unit(2, "cm"),
    padding = unit(c(20, 20, 60, 20), "mm")  # 上、右、下、左
  )
  
  dev.off()
  
  
  #################### 12. 保存数据对象 ####################
  result_data <- list(
    matrix_raw = kegg_interest_sel,
    matrix_scaled = kegg_interest_scaled_sel,
    annotation = col_annot_sel,
    colors = annotation_colors_sel
  )
  
  all_results[[group_var]] <- result_data
  
  # RDS
  data_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", tag, "_", safe_group_var, "_data.rds"))
  saveRDS(result_data, file = data_file)
  
  # CSV（原始 GSVA 分数）
  csv_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", tag, "_", safe_group_var, "_scores.csv"))
  write.csv(kegg_interest_sel, file = csv_file)
  
  # 注释信息
  annot_file <- file.path(output_dir, paste0("GSVA_KEGG_celltype_", tag, "_", safe_group_var, "_annotation.csv"))
  write.csv(col_annot_sel, file = annot_file)
}








#################### 3. 按多维度组合分组求平均 GSVA 得分 ####################

# 创建多维度组合分组
meta$combined_group <- paste0(
  meta$celltype, "_",      # 将celltype移到第一位
  meta$stage, "_",
  meta$treatment, "_",
  meta$Type, "_",
  meta$`BRCA1/2`, "_",  
  meta$P53, "_",
  meta$`PD-L1`, "_",    
  meta$HER2
)

# 简化组合组名，便于展示
cluster_vec <- meta$combined_group
names(cluster_vec) <- rownames(meta)

# 确保列名与 meta 行名对应
kegg11 <- kegg[, names(cluster_vec), drop = FALSE]

# 计算每个组合组的平均得分
kegg_combined <- sapply(
  X = unique(cluster_vec),
  FUN = function(cl){
    cells <- names(cluster_vec)[cluster_vec == cl]
    rowMeans(kegg11[, cells, drop = FALSE])
  }
)

# 转成矩阵，行为通路，列为组合组
kegg_combined <- as.matrix(kegg_combined)

#################### 4. 重新排列列顺序，按照指定顺序排序 ####################

# 解析组合组名，提取各个维度信息
parsed_groups <- strsplit(colnames(kegg_combined), "_")

# 创建数据框存储解析后的信息
group_info <- data.frame(
  matrix(unlist(parsed_groups), ncol = 8, byrow = TRUE),
  stringsAsFactors = FALSE
)
colnames(group_info) <- c("celltype", "stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")
group_info$full_name <- colnames(kegg_combined)

# 按照指定的顺序排序：celltype > stage > treatment > Type > BRCA1_2 > P53 > PD_L1 > HER2
group_info <- group_info[order(
  group_info$celltype,
  group_info$stage,
  group_info$treatment,
  group_info$Type,
  group_info$BRCA1_2,
  group_info$P53,
  group_info$PD_L1,
  group_info$HER2
), ]

# 重新排列矩阵的列
kegg_combined <- kegg_combined[, group_info$full_name, drop = FALSE]

# 定义感兴趣的通路列表（保持原样）
kegg_interest_ids <- c(
  "JAK_STAT_SIGNALING_PATHWAY",
  "PI3K_AKT_SIGNALING_PATHWAY", "ABC_TRANSPORTERS", "GLUTATHIONE_METABOLISM", 
  "METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450", "DRUG_METABOLISM_OTHER_ENZYMES", 
  "PATHWAYS_IN_CANCER", "MAPK_SIGNALING_PATHWAY", 
  "APOPTOSIS", "MTOR_SIGNALING_PATHWAY", "FOCAL_ADHESION", 
  "TGF_BETA_SIGNALING_PATHWAY", "HEDGEHOG_SIGNALING_PATHWAY", "WNT_SIGNALING_PATHWAY", 
  "CELL_ADHESION_MOLECULES_CAMS", "GLYCOLYSIS_GLUCONEOGENESIS","TIGHT_JUNCTION","OXIDATIVE_PHOSPHORYLATION",
  
  "T_CELL_RECEPTOR_SIGNALING_PATHWAY", "CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "TRYPTOPHAN_METABOLISM", "ARGININE_AND_PROLINE_METABOLISM", 
  "PURINE_METABOLISM", "PYRIMIDINE_METABOLISM",
  "TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","PPAR_SIGNALING_PATHWAY",
  "NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY", "ANTIGEN_PROCESSING_AND_PRESENTATION"
  
)

# 找出实际存在于矩阵中的通路
existing_pathways <- kegg_interest_ids[kegg_interest_ids %in% rownames(kegg_combined)]

# 输出未找到的通路
missing_pathways <- setdiff(kegg_interest_ids, rownames(kegg_combined))
if (length(missing_pathways) > 0) {
  cat("以下通路在矩阵中未找到:\n")
  print(missing_pathways)
}

# 使用实际存在的通路创建子集矩阵
if (length(existing_pathways) > 0) {
  kegg_interest <- kegg_combined[existing_pathways, , drop = FALSE]
  cat("成功选择了", length(existing_pathways), "个感兴趣的通路\n")
} else {
  stop("没有找到任何感兴趣的通路，请检查通路名称")
}


#################### 5. 准备热图注释信息 ####################
# 使用解析后的group_info数据框
col_annot <- group_info[, c("celltype", "stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")]
rownames(col_annot) <- colnames(kegg_interest)

# 设置颜色方案
library(RColorBrewer)

# 为每个分组变量生成颜色
generate_colors <- function(categories, palette_name = "Set3") {
  n_categories <- length(categories)
  if (n_categories <= 12) {
    colors <- brewer.pal(max(3, n_categories), palette_name)[1:n_categories]
  } else {
    colors <- colorRampPalette(brewer.pal(12, palette_name))(n_categories)
  }
  names(colors) <- categories
  return(colors)
}

# 优化celltype颜色生成：使用更丰富的颜色方案
generate_celltype_colors <- function(celltypes) {
  n_celltypes <- length(celltypes)
  sorted_celltypes <- sort(celltypes)  # 排序以确保颜色一致性
  
  # 预定义多个颜色调色板，用于组合
  color_palettes <- list(
    c("#DC050C",  "#1965B0","#FB8072", "#7BAFDE", "#882E72",
      "#B17BA6", "#FF7F00", "#FDB462", "#E7298A","#9999FF","#00CC66",
      "#A4CDE1",'#FF9999',"#66CCCC",'#4F6272',"#FF3366","#CC0066","#CC99CC","#FFCCCC","#CCFFCC","#FFFFCC",'#E5D2DD','#58A4C3',
      '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
      "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
      "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
      "#6699CC","#9999FF","#CCCCFF","#FF6699","#6699CC","#FFFFCC"),
    c("#FF6B6B",  "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD",  "#F7DC6F"),
    c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6", "#1ABC9C", "#D35400", "#34495E"),
    c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F"),
    c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00"),
    c("#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4","#4ECDC4"),
    c("#FED9A6", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC")
  )
  
  if (n_celltypes <= 8) {
    # 如果细胞类型较少，使用第一种调色板
    colors <- color_palettes[[1]][1:n_celltypes]
  } else if (n_celltypes <= 16) {
    # 中等数量的细胞类型，组合两个调色板
    colors <- c(color_palettes[[1]], color_palettes[[2]])[1:n_celltypes]
  } else if (n_celltypes <= 24) {
    # 较多细胞类型，组合三个调色板
    colors <- c(color_palettes[[1]], color_palettes[[2]], color_palettes[[3]])[1:n_celltypes]
  } else if (n_celltypes <= 32) {
    # 很多细胞类型，组合四个调色板
    colors <- c(color_palettes[[1]], color_palettes[[2]], color_palettes[[3]], color_palettes[[4]])[1:n_celltypes]
  } else if (n_celltypes <= 40) {
    # 非常多细胞类型，组合五个调色板
    colors <- c(color_palettes[[1]], color_palettes[[2]], color_palettes[[3]], 
                color_palettes[[4]], color_palettes[[5]])[1:n_celltypes]
  } else {
    # 极端情况：使用彩虹色
    colors <- rainbow(n_celltypes)
  }
  
  names(colors) <- sorted_celltypes
  return(colors)
}

# 为其他分组变量生成颜色（使用原函数）
generate_other_colors <- function(categories, palette_name = "Set3") {
  n_categories <- length(categories)
  if (n_categories <= 12) {
    colors <- brewer.pal(max(3, n_categories), palette_name)[1:n_categories]
  } else {
    colors <- colorRampPalette(brewer.pal(12, palette_name))(n_categories)
  }
  names(colors) <- categories
  return(colors)
}

# 获取所有唯一的细胞类型
unique_celltypes <- sort(unique(col_annot$celltype))
n_celltypes <- length(unique_celltypes)

# 为每个分组变量生成颜色
annotation_colors <- list(
  celltype = generate_celltype_colors(unique_celltypes),  # 使用优化的celltype颜色生成函数
  stage = generate_other_colors(unique(col_annot$stage), "Set3"),
  treatment = generate_other_colors(unique(col_annot$treatment), "Set1"),
  Type = generate_other_colors(unique(col_annot$Type), "Pastel1"),
  BRCA1_2 = generate_other_colors(unique(col_annot$BRCA1_2), "Pastel2"),
  P53 = generate_other_colors(unique(col_annot$P53), "Accent"),
  PD_L1 = generate_other_colors(unique(col_annot$PD_L1), "Dark2"),
  HER2 = generate_other_colors(unique(col_annot$HER2), "Paired")
)


#################### 6. 画热图（展示感兴趣的通路） ####################
library(pheatmap)

# 对数据进行标准化（按行标准化）
kegg_interest_scaled <- t(scale(t(kegg_interest)))

# 确定celltype分组分隔位置
celltype_labels <- col_annot$celltype
gaps_positions <- c()
for(i in 1:(length(celltype_labels)-1)) {
  if(celltype_labels[i] != celltype_labels[i+1]) {
    gaps_positions <- c(gaps_positions, i)
  }
}

# 确保gaps_positions不超过矩阵列数
if(length(gaps_positions) > 0) {
  gaps_positions <- gaps_positions[gaps_positions < ncol(kegg_interest_scaled)]
}

# 计算合适的图形尺寸
n_pathways <- nrow(kegg_interest)
n_samples <- ncol(kegg_interest)

# 绘制主热图
p <- pheatmap(
  kegg_interest_scaled,
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,  # 不聚类列，保持排序顺序
  scale          = "none",
  annotation_col = col_annot,
  annotation_colors = annotation_colors,
  show_rownames  = TRUE,
  show_colnames  = FALSE,
  fontsize_row   = ifelse(n_pathways > 30, 8, 10),
  fontsize_col   = 12,
  color          = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main           = paste0("GSVA KEGG Selected Pathways (", n_pathways, " pathways)\nSorted by: celltype > stage > treatment > Type > BRCA1_2 > P53 > PD_L1 > HER2"),
  border_color   = NA,
  cellwidth      = ifelse(n_samples > 30, 8, ifelse(n_samples > 20, 12, 15)),
  cellheight     = ifelse(n_pathways > 40, 8, 10),
  gaps_col       = if(length(gaps_positions) > 0) gaps_positions else NULL,
  filename       = file.path(output, "GSVA_KEGG_multidimensional_analysis.pdf"),
  # 整体尺寸设置 - 修改了height参数，增大了热图高度
  width = unit(max(15, ncol(kegg_top50_scaled) * 0.4 + 4), "cm"),
  height = unit(max(6, nrow(kegg_top50_scaled) * 0.3 + 4), "cm")  # 从0.15增加到0.25，基础值从3增加到4
)




#################### 7. 使用ComplexHeatmap绘制主热图（图例在底部） ####################

# 确保列注释的顺序与排序顺序一致
col_annot_ordered <- col_annot[, c("celltype", "stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")]

# 加载ComplexHeatmap包
if (!require(ComplexHeatmap)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
  library(ComplexHeatmap)
}

# 加载circlize包用于颜色映射
if (!require(circlize)) {
  install.packages("circlize")
  library(circlize)
}

# 创建颜色映射函数
col_fun <- colorRamp2(
  seq(min(kegg_interest_scaled, na.rm = TRUE), 
      max(kegg_interest_scaled, na.rm = TRUE), 
      length = 100),
  colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)

# 创建列注释 - 修改图例参数
ha_column <- HeatmapAnnotation(
  df = col_annot,
  col = annotation_colors,
  annotation_name_side = "left",
  annotation_legend_param = list(
    # 设置图例为两行显示
    nrow = 3,
    ncol = ceiling(length(annotation_colors) / 3),  # 计算每列显示的数量
    direction = "horizontal",
    
    # 增大图例标题大小
    title_gp = gpar(fontsize = 16, fontface = "bold"),  # 从10增大到14
    
    # 增大图例标签大小
    labels_gp = gpar(fontsize = 14),  # 从9增大到12
    
    # 增大图例图形大小
    legend_height = unit(1.5, "cm"),  # 增大图例高度
    legend_width = unit(1.5, "cm"),   # 增大图例宽度
    
    # 设置图例之间的间距
    gap = unit(0.5, "cm"),
    
    title_position = "topcenter",
    
    # 为每个图例单独设置布局
    celltype = list(
      title = "Cell Type",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    stage = list(
      title = "Stage",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    treatment = list(
      title = "Treatment",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    Type = list(
      title = "Type",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    BRCA1_2 = list(
      title = "BRCA1/2",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    P53 = list(
      title = "P53",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    PD_L1 = list(
      title = "PD-L1",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    ),
    HER2 = list(
      title = "HER2",
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    )
  ),
  gap = unit(1, "mm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16, fontface = "bold")  # 增大注释名称大小
)

# 确定分组分隔位置（按celltype变化）
celltype_labels <- col_annot$celltype
column_gaps <- c()
for(i in 1:(length(celltype_labels)-1)) {
  if(celltype_labels[i] != celltype_labels[i+1]) {
    column_gaps <- c(column_gaps, i)
  }
}

# 创建热图 - 修改热图图例参数
ht_main <- Heatmap(
  kegg_interest_scaled,
  name = "Z-score",
  col = col_fun,
  
  # 行设置
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = ifelse(n_pathways > 30, 16, 18), fontface = "bold"),  # 增大行名字体
  row_dend_width = unit(25, "mm"),  # 增大行聚类宽度
  row_title = "Pathways",
  row_title_gp = gpar(fontsize = 18, fontface = "bold"),  # 增大行标题
  
  # 列设置
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_title = "Samples",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),  # 增大列标题
  
  # 列分组
  column_split = if(length(column_gaps) > 0) {
    # 创建分组因子
    split_factor <- rep(1, ncol(kegg_interest_scaled))
    for(i in seq_along(column_gaps)) {
      split_factor[(column_gaps[i]+1):length(split_factor)] <- i + 1
    }
    as.factor(split_factor)
  } else {
    NULL
  },
  column_gap = if(length(column_gaps) > 0) unit(3, "mm") else unit(0, "mm"),  # 增大列间距
  
  # 注释设置
  top_annotation = ha_column,
  
  # 热图主体设置
  rect_gp = gpar(col = NA),
  border = FALSE,
  use_raster = ifelse(ncol(kegg_interest_scaled) > 50, TRUE, FALSE),
  raster_quality = 2,
  
  # 热图图例设置 - 增大尺寸
  heatmap_legend_param = list(
    title = "GSVA Z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(14, "cm"),  # 增大图例宽度
    legend_height = unit(1.8, "cm"),  # 增大图例高度
    title_gp = gpar(fontsize = 18, fontface = "bold"),  # 显著增大标题
    labels_gp = gpar(fontsize = 16),  # 显著增大标签
    at = round(seq(min(kegg_interest_scaled), max(kegg_interest_scaled), length.out = 5), 1),
    grid_border = "black",
    # 设置颜色条大小
    color_bar = "continuous",
    # 设置图例标签位置
    labels_rot = 0,
    # 设置图例边框
    border = "black",
    # 设置图例内部间距
    padding = unit(c(10, 10, 10, 10), "mm")
  ),
  
  # 整体尺寸设置 - 修改了height参数，增大了热图高度
  width = unit(max(60, ncol(kegg_top50_scaled) * 0.4 + 4), "cm"),
  height = unit(max(20, nrow(kegg_top50_scaled) * 0.3 + 4), "cm")  # 从0.15增加到0.25，基础值从3增加到4
)

# 绘制热图并保存为PDF
pdf(file.path(output, "GSVA_KEGG_multidimensional_analysis_complexheatmap.pdf"),
    width = max(14, n_samples * 0.2 + 10),  # 增大PDF宽度
    height = max(16, n_pathways * 0.12 + 12))  # 增大PDF高度

# 绘制热图，将图例分为两行
draw(ht_main,
     heatmap_legend_side = "bottom",  # 热图图例在底部
     annotation_legend_side = "bottom",  # 注释图例在底部
     merge_legends = TRUE,  # 合并图例
     # 设置图例包为两行
     ht_gap = unit(3, "cm"),  # 图例之间的间距
     # 增大底部边距以容纳两行图例
     padding = unit(c(25, 25, 80, 25), "mm"),  # 上右下左（底部从60增大到80）
     # 调整图例排列
     legend_grouping = "adjusted"
)

dev.off()





#################### 7. 分celltype绘制子热图进行比较分析 ####################

# 为了更清晰地比较每个celltype中其他维度的变化，我们为每个celltype单独绘制热图
unique_celltypes <- unique(col_annot$celltype)

for (celltype in unique_celltypes) {
  # 选择当前celltype的样本
  celltype_samples <- rownames(col_annot)[col_annot$celltype == celltype]
  
  if (length(celltype_samples) > 1) {  # 至少需要2个样本才能绘制热图
    # 提取当前celltype的数据
    celltype_data <- kegg_interest_scaled[, celltype_samples, drop = FALSE]
    celltype_annot <- col_annot[celltype_samples, , drop = FALSE]
    
    # 重新排列样本顺序，以便更好地展示其他维度的变化
    # 按treatment, stage, Type, BRCA1_2, P53, PD_L1, HER2的顺序排序
    celltype_annot$order <- with(celltype_annot, 
                                 paste(treatment, stage, Type, BRCA1_2, P53, PD_L1, HER2))
    celltype_data <- celltype_data[, order(celltype_annot$order), drop = FALSE]
    celltype_annot <- celltype_annot[order(celltype_annot$order), , drop = FALSE]
    
    # 确定分组分隔位置（按treatment变化）
    treatment_labels <- celltype_annot$treatment
    celltype_gaps <- c()
    for(i in 1:(length(treatment_labels)-1)) {
      if(treatment_labels[i] != treatment_labels[i+1]) {
        celltype_gaps <- c(celltype_gaps, i)
      }
    }
    
    # 使用ComplexHeatmap绘图
    pdf(file.path(output, paste0("GSVA_KEGG_celltype_", gsub("[^[:alnum:]]", "_", celltype), ".pdf")),
        width = max(8, length(celltype_samples) * 0.3 + 6),
        height = max(10, n_pathways * 0.25 + 2))
    
    # 创建热图
    ht <- Heatmap(
      celltype_data,
      name = "Expression",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_split = factor(celltype_annot$treatment, levels = unique(celltype_annot$treatment)),
      top_annotation = HeatmapAnnotation(
        df = celltype_annot[, c("stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")],
        col = annotation_colors,
        show_annotation_name = TRUE
      ),
      show_row_names = TRUE,
      show_column_names = FALSE,
      column_names_gp = gpar(fontsize = 8),
      col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      column_title = paste0(celltype),
      border = FALSE,
      # 整体尺寸设置 - 修改了height参数，增大了热图高度
      width = unit(max(15, ncol(kegg_top50_scaled) * 0.4 + 4), "cm"),
      height = unit(max(12, nrow(kegg_top50_scaled) * 0.3 + 4), "cm")  # 从0.15增加到0.25，基础值从3增加到4
      
    )
    
    # 使用draw函数绘制热图，并指定图例位置
    draw(ht, 
         heatmap_legend_side = "bottom",  # 热图颜色条图例放在右侧
         annotation_legend_side = "bottom",  # 注释图例放在右侧
         merge_legend = TRUE  # 合并图例
    )
    dev.off()
  }
}



#################### 7. 分celltype绘制子热图进行比较分析 ####################

# 为了更清晰地比较每个celltype中其他维度的变化，我们为每个celltype单独绘制热图
unique_celltypes <- unique(col_annot$celltype)

for (celltype in unique_celltypes) {
  # 选择当前celltype的样本
  celltype_samples <- rownames(col_annot)[col_annot$celltype == celltype]
  
  if (length(celltype_samples) > 1) {  # 至少需要2个样本才能绘制热图
    # 提取当前celltype的数据
    celltype_data <- kegg_interest_scaled[, celltype_samples, drop = FALSE]
    celltype_annot <- col_annot[celltype_samples, , drop = FALSE]
    
    # 重新排列样本顺序，以便更好地展示其他维度的变化
    # 按treatment, stage, Type, BRCA1_2, P53, PD_L1, HER2的顺序排序
    celltype_annot$order <- with(celltype_annot, 
                                 paste(treatment, stage, Type, BRCA1_2, P53, PD_L1, HER2))
    celltype_data <- celltype_data[, order(celltype_annot$order), drop = FALSE]
    celltype_annot <- celltype_annot[order(celltype_annot$order), , drop = FALSE]
    
    # 确定分组分隔位置（按treatment变化）
    treatment_labels <- celltype_annot$treatment
    celltype_gaps <- c()
    for(i in 1:(length(treatment_labels)-1)) {
      if(treatment_labels[i] != treatment_labels[i+1]) {
        celltype_gaps <- c(celltype_gaps, i)
      }
    }
    
    # 计算动态尺寸以确保完整显示
    # 1. 计算行名所需的最大宽度（基于字符数和字体大小）
    max_rowname_length <- max(nchar(rownames(celltype_data)))
    
    
    # 使用ComplexHeatmap绘图，确保行名完整显示
    pdf(file.path(output, paste0("GSVA_KEGG_celltype_", gsub("[^[:alnum:]]", "_", celltype), ".pdf")),
        width = max(8, length(celltype_samples) * 0.4 + 6),
        height = max(10, n_pathways * 0.1 + 2))
    
    # 创建热图，确保行名完整显示
    ht <- Heatmap(
      celltype_data,
      name = "Expression",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_title = paste0(celltype),
      column_split = factor(celltype_annot$treatment, levels = unique(celltype_annot$treatment)),
      top_annotation = HeatmapAnnotation(
        df = celltype_annot[, c("stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")],
        col = annotation_colors,
        show_annotation_name = TRUE,
        annotation_name_gp = gpar(fontsize = 12)
      ),
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 10),  # 设置行名字体大小
      #row_names_max_width = unit(row_name_width, "cm"),  # 确保行名区域足够宽
      show_column_names = FALSE,
      col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      border = TRUE,
      border_gp = gpar(col = "gray50", lwd = 0.5),
      heatmap_legend_param = list(
        title = "GSVA Score",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        legend_width = unit(3, "cm"),
        legend_height = unit(2, "cm"),
        direction = "horizontal"
      ),
      # 整体尺寸设置 - 修改了height参数，增大了热图高度
      width = unit(max(15, ncol(kegg_top50_scaled) * 0.4 + 4), "cm"),
      height = unit(max(12, nrow(kegg_top50_scaled) * 0.3 + 4), "cm")  # 从0.15增加到0.25，基础值从3增加到4
      
    )
    
    # 使用draw函数绘制热图，确保所有元素完整显示
    draw(ht, 
         heatmap_legend_side = "bottom",
         annotation_legend_side = "bottom",
         merge_legend = TRUE,
         padding = unit(c(2, 2, 2, 2), "cm")  # 增加边距确保完整显示
    )
    dev.off()
    
  }
}








################## 挑选细胞 ##########################
#################### 5. 准备热图注释信息 ####################
# 使用解析后的group_info数据框
col_annot <- group_info[, c("celltype", "stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")]
rownames(col_annot) <- colnames(kegg_interest)


#################### 5. 准备热图注释信息 ####################
# 使用解析后的group_info数据框
col_annot <- group_info[, c("celltype", "stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")]
rownames(col_annot) <- colnames(kegg_interest)

## ★ 只保留 Treg / Tex / CTL 三种 celltype -----------------
selected_celltypes <- c("Treg", "Tex", "CTL")

keep_cols <- col_annot$celltype %in% selected_celltypes

# 子集化矩阵和注释
kegg_interest_sel <- kegg_interest[, keep_cols, drop = FALSE]
col_annot_sel    <- droplevels(col_annot[keep_cols, , drop = FALSE])

# 如果你希望在每个 celltype 内继续按 stage/treatment/... 排一下列顺序，可以加：
col_order <- order(
  col_annot_sel$celltype,
  col_annot_sel$stage,
  col_annot_sel$treatment,
  col_annot_sel$Type,
  col_annot_sel$BRCA1_2,
  col_annot_sel$P53,
  col_annot_sel$PD_L1,
  col_annot_sel$HER2
)
kegg_interest_sel <- kegg_interest_sel[, col_order, drop = FALSE]
col_annot_sel    <- col_annot_sel[col_order, , drop = FALSE]


# 获取所有唯一的细胞类型
unique_celltypes <- sort(unique(col_annot_sel$celltype))
n_celltypes <- length(unique_celltypes)

# 为每个分组变量生成颜色
annotation_colors_sel <- list(
  celltype  = generate_celltype_colors(unique_celltypes),  # 使用优化的celltype颜色生成函数
  stage     = generate_other_colors(unique(col_annot_sel$stage), "Set3"),
  treatment = generate_other_colors(unique(col_annot_sel$treatment), "Set1"),
  Type      = generate_other_colors(unique(col_annot_sel$Type), "Pastel1"),
  BRCA1_2   = generate_other_colors(unique(col_annot_sel$BRCA1_2), "Pastel2"),
  P53       = generate_other_colors(unique(col_annot_sel$P53), "Accent"),
  PD_L1     = generate_other_colors(unique(col_annot_sel$PD_L1), "Dark2"),
  HER2      = generate_other_colors(unique(col_annot_sel$HER2), "Paired")
)


#################### 6. 画热图（只展示 Treg / Tex / CTL 中的感兴趣通路） ####################
library(pheatmap)

# 对数据进行标准化（按行标准化）
kegg_interest_scaled_sel <- t(scale(t(kegg_interest_sel)))

# 计算维度
n_pathways_sel <- nrow(kegg_interest_sel)
n_samples_sel  <- ncol(kegg_interest_sel)

# 确定 celltype 分组分隔位置
celltype_labels_sel <- col_annot_sel$celltype
gaps_positions_sel <- c()
for(i in 1:(length(celltype_labels_sel) - 1)) {
  if(celltype_labels_sel[i] != celltype_labels_sel[i + 1]) {
    gaps_positions_sel <- c(gaps_positions_sel, i)
  }
}
if(length(gaps_positions_sel) > 0) {
  gaps_positions_sel <- gaps_positions_sel[gaps_positions_sel < ncol(kegg_interest_scaled_sel)]
}

# 绘制主热图（只包含 Treg / Tex / CTL）
p <- pheatmap(
  kegg_interest_scaled_sel,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,  # 不聚类列，保持排序顺序
  scale             = "none",
  annotation_col    = col_annot_sel,
  annotation_colors = annotation_colors_sel,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = ifelse(n_pathways_sel > 30, 8, 10),
  fontsize_col      = 12,
  color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main              = paste0(
    "GSVA KEGG Selected Pathways in Treg / Tex / CTL (",
    n_pathways_sel, " pathways)\n",
    "Sorted by: celltype > stage > treatment > Type > BRCA1_2 > P53 > PD_L1 > HER2"
  ),
  border_color      = NA,
  cellwidth         = ifelse(n_samples_sel > 30, 8, ifelse(n_samples_sel > 20, 12, 15)),
  cellheight        = ifelse(n_pathways_sel > 40, 8, 10),
  gaps_col          = if(length(gaps_positions_sel) > 0) gaps_positions_sel else NULL,
  filename          = file.path(output, "GSVA_KEGG_multidimensional_Treg_Tex_CTL_pheatmap.pdf"),
  # 整体尺寸设置 - 修改了height参数，增大了热图高度
  width = unit(max(15, ncol(kegg_top50_scaled) * 0.4 + 4), "cm"),
  height = unit(max(12, nrow(kegg_top50_scaled) * 0.3 + 4), "cm")  # 从0.15增加到0.25，基础值从3增加到4
  
)



#################### 7. 使用ComplexHeatmap绘制 Treg / Tex / CTL 主热图 ####################

# 对筛选后的数据进行标准化
kegg_interest_scaled_sel <- t(scale(t(kegg_interest_sel)))

# 颜色映射函数
col_fun_sel <- colorRamp2(
  seq(min(kegg_interest_scaled_sel, na.rm = TRUE), 
      max(kegg_interest_scaled_sel, na.rm = TRUE), 
      length = 100),
  colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)

# 列注释（仍按原顺序）
col_annot_ordered_sel <- col_annot_sel[, c("celltype", "stage", "treatment", "Type", "BRCA1_2", "P53", "PD_L1", "HER2")]

# 列注释对象
ha_column_sel <- HeatmapAnnotation(
  df  = col_annot_ordered_sel,
  col = annotation_colors_sel,
  annotation_name_side = "left",
  annotation_legend_param = list(
    nrow = 3,
    direction = "horizontal",
    title_gp  = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    legend_height = unit(1.5, "cm"),
    legend_width  = unit(1.5, "cm"),
    gap = unit(0.5, "cm"),
    title_position = "topcenter"
  ),
  gap = unit(1, "mm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16, fontface = "bold")
)

# 按 celltype 计算列分隔位置
celltype_labels_sel <- col_annot_sel$celltype
column_gaps_sel <- c()
for(i in 1:(length(celltype_labels_sel) - 1)) {
  if(celltype_labels_sel[i] != celltype_labels_sel[i + 1]) {
    column_gaps_sel <- c(column_gaps_sel, i)
  }
}

n_pathways_sel <- nrow(kegg_interest_scaled_sel)
n_samples_sel  <- ncol(kegg_interest_scaled_sel)

# 热图对象
ht_main_sel <- Heatmap(
  kegg_interest_scaled_sel,
  name = "Z-score",
  col  = col_fun_sel,
  
  # 行设置
  cluster_rows   = TRUE,
  show_row_names = TRUE,
  row_names_gp   = gpar(fontsize = ifelse(n_pathways_sel > 30, 16, 18), fontface = "bold"),
  row_dend_width = unit(25, "mm"),
  row_title      = "Pathways",
  row_title_gp   = gpar(fontsize = 18, fontface = "bold"),
  
  # 列设置
  cluster_columns   = FALSE,
  show_column_names = FALSE,
  column_title      = "Samples",
  column_title_gp   = gpar(fontsize = 14, fontface = "bold"),
  
  # 列分组（按 celltype 分块）
  column_split = if(length(column_gaps_sel) > 0) {
    split_factor <- rep(1, ncol(kegg_interest_scaled_sel))
    for(i in seq_along(column_gaps_sel)) {
      split_factor[(column_gaps_sel[i] + 1):length(split_factor)] <- i + 1
    }
    as.factor(split_factor)
  } else {
    NULL
  },
  column_gap = if(length(column_gaps_sel) > 0) unit(3, "mm") else unit(0, "mm"),
  
  # 注释
  top_annotation = ha_column_sel,
  
  rect_gp      = gpar(col = NA),
  border       = FALSE,
  use_raster   = ifelse(ncol(kegg_interest_scaled_sel) > 50, TRUE, FALSE),
  raster_quality = 2,
  
  heatmap_legend_param = list(
    title           = "GSVA Z-score",
    title_position  = "topcenter",
    legend_direction = "horizontal",
    legend_width    = unit(14, "cm"),
    legend_height   = unit(1.8, "cm"),
    title_gp        = gpar(fontsize = 18, fontface = "bold"),
    labels_gp       = gpar(fontsize = 16),
    at              = round(seq(min(kegg_interest_scaled_sel), 
                                max(kegg_interest_scaled_sel), 
                                length.out = 5), 1),
    grid_border     = "black",
    color_bar       = "continuous",
    labels_rot      = 0,
    border          = "black",
    padding         = unit(c(10, 10, 10, 10), "mm")
  ),
  
  width  = unit(max(28, n_samples_sel * 0.5), "cm"),
  height = unit(max(25, n_pathways_sel * 0.35), "cm")
)

# 保存 PDF
pdf(file.path(output, "GSVA_KEGG_multidimensional_Treg_Tex_CTL_complexheatmap.pdf"),
    width  = max(28, n_samples_sel * 0.2 + 10),
    height = max(16, n_pathways_sel * 0.12 + 12))

draw(ht_main_sel,
     heatmap_legend_side      = "bottom",
     annotation_legend_side   = "bottom",
     merge_legends            = TRUE,
     ht_gap                   = unit(3, "cm"),
     padding                  = unit(c(25, 25, 80, 25), "mm")
)

dev.off()










