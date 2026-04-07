############################# Integration ####################################
# Load Seurat objects for 8C and 4C
seurat_8C <- readRDS("/out(WF)/celltype(8C).rds")
seurat_4C <- readRDS("/out/sce_female(e4CL-primed).rds")
seurat_E0 <- readRDS("/out/sce_female(E0-E6).rds")
seurat_E3 <- readRDS("/out/sce_female(E3-E7).rds")

# Modify orig.ident and treatment fields in seurat_8C
seurat_8C$orig.ident <- "8CLC"
seurat_8C$treatment <- "8CLC"

# Extract cells with celltype == "8CLC"
seurat_8C_8CLC <- subset(seurat_8C, subset = celltype == "8CLC")
table(seurat_8C_8CLC@meta.data$celltype)

# Add a new column "source"
seurat_8C_8CLC$source <- "This study"
seurat_4C$source <- "Mazid et al."
seurat_E0$source <- "Yan et al."
seurat_E3$source <- "Petropoulos et al."

# Create object list
seurat.list <- list(seurat_8C_8CLC, seurat_4C, seurat_E3)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)

# Perform data integration
combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:20)

# Normalization and dimensional reduction
DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(combined_seurat, npcs = 30, verbose = FALSE)
combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:20)

# Add celltype column based on treatment
combined_seurat@meta.data$celltype <- combined_seurat@meta.data$treatment
combined_seurat@meta.data$seurat_clusters <- combined_seurat@meta.data$treatment

# Remove columns containing NA values
combined_seurat@meta.data <- combined_seurat@meta.data[, colSums(is.na(combined_seurat@meta.data)) == 0]

View(combined_seurat@meta.data)

# Save integrated Seurat object
saveRDS(combined_seurat, file = "/out/combined_seurat_8c.rds")

####################### Seurat Analysis #####################

# Set working directory and output path
setwd("/out/")
outdir <- "/out/"

MergeOUT <- paste(outdir, 'Merge', sep='/')
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, 'Multiple_', sep='/')

# Load integrated object
file_path <- file.path(outdir, "combined_seurat_8c.rds")
ScRNA <- readRDS(file_path)

# Reorder celltype
celltype_order <- c("E3","E4","E5","E6","E7",
                    "direct e4CL-D3","direct e4CL-D5","direct e4CL-D7",
                    "4CL-D1","4CL-D2","4CL-D3","4CL-D5","4CL-D8","4CL-D12",
                    "e4CL-D1","e4CL-D2","e4CL-D3","e4CL-D5",
                    "8CLC")

ScRNA$celltype <- factor(ScRNA$celltype, levels = celltype_order)

# Reorder source
source_order <- c("This study", "Mazid et al.", "Petropoulos et al.")
ScRNA$source <- factor(ScRNA$source, levels = source_order)

View(ScRNA@meta.data)

# Violin plot for QC metrics
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 12, height = 6)
VlnPlot(ScRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rps"), ncol = 4, group.by = "source", pt.size = 0)
dev.off()

# Split UMAP by source
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 20, height = 6)
DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, split.by = "source", group.by = "celltype")
dev.off()

# UMAP plots
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 8, height = 6)
DimPlot(ScRNA, reduction = "umap", group.by = "celltype")
dev.off()

# Highlight cells by source
p1 <- DimPlot(ScRNA, reduction = "umap", group.by = "celltype", cells.highlight = list("This study" = WhichCells(ScRNA, expression = source == "This study")), cols = "lightgrey")
p2 <- DimPlot(ScRNA, reduction = "umap", group.by = "celltype", cells.highlight = list("Mazid et al." = WhichCells(ScRNA, expression = source == "Mazid et al.")), cols = "lightgrey")
p3 <- DimPlot(ScRNA, reduction = "umap", group.by = "celltype", cells.highlight = list("Petropoulos et al." = WhichCells(ScRNA, expression = source == "Petropoulos et al.")), cols = "lightgrey")

combined_plot <- p1 + p2 + p3

pdf(paste0(OUTPUT, "cluster_umap_by_source.pdf"), width = 12, height = 4.5)
print(combined_plot)
dev.off()

########### Manual cell annotation ########

# Load data
file_path <- file.path(outdir, "ScRNA（分群后）.rds")
scedata <- readRDS(file_path)

# Extract TPRX1 expression
TPRX1_expression <- scedata[["RNA"]]@data["TPRX1", ]

# Define celltype based on expression
scedata$celltype <- ifelse(TPRX1_expression > 0, "8CLC", "non-8CLC")

# Plot UMAP
pdf(paste(outdir, "celltype/ann_umap_8c.pdf", sep='/'), width = 6, height = 5)
DimPlot(object=scedata, group.by = "celltype", reduction='umap')
dev.off()

# Calculate expression ratio
expressed_cells <- sum(scedata$celltype == "8CLC")
total_cells <- nrow(scedata@meta.data)
expression_ratio <- expressed_cells / total_cells * 100

# Plot FeaturePlot
pdf(paste0(outdir, "/celltype/TPRX1_FeaturePlot_umap.pdf"), width = 4, height = 4)
FeaturePlot(scedata, features = "TPRX1", reduction = "umap")
dev.off()

# Marker gene analysis
cellmarker <- c("TPRX1","SOX2","NANOG","KLF4")
cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

plot <- DotPlot(scedata, features = unique(cellmarker), group.by = "celltype")

ggsave(filename = paste(outdir, "celltype/marker_DotPlot.pdf", sep='/'), plot = plot, width = 10, height = 8)

# Violin plot
library(reshape2)
vln.df <- as.data.frame(scedata[["RNA"]]@data[cellmarker,])
vln.df$gene <- rownames(vln.df)
vln.df <- melt(vln.df, id = "gene")
colnames(vln.df)[c(2,3)] <- c("CB", "exp")

anno <- scedata@meta.data[, c("CB", "celltype")]
vln.df <- inner_join(vln.df, anno, by = "CB")

plot <- ggplot(vln.df, aes(exp, celltype)) + geom_violin()

ggsave(filename = paste(outdir, "celltype/marker_ViolinPlot.pdf", sep='/'), plot = plot, width = 12, height = 6)


