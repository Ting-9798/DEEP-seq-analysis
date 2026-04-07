# Clear environment
rm(list = ls())

# Load required libraries
library(ggplot2)
library(Matrix)
library(dplyr)

# Set working directory
setwd("/data/")
folders <- c("A", "B","C")
sample_names <- c("A", "B","C")  # Corresponding sample names

outdir <- "/out/"

# Define thresholds for each sample (thresholds you provided)
thresholds <- c("A" = 5000, "B" = 1000, "C" = 10000)

# Specify sample highlight colors (consistent with original)
sample_colors <- c("A" = "#882E72", "B" = "#D0F199", "C" = "#31CDEE")
grey_color <- "grey80"

# Store data for all samples
barcode_rank_all <- data.frame()

for (i in seq_along(folders)) {
  folder <- folders[i]
  sample_name <- sample_names[i]
  
  # Build file paths
  barcode_file <- paste0(folder, "barcodes.tsv")
  feature_file <- paste0(folder, "features.tsv")
  matrix_file <- paste0(folder, "matrix.mtx")
  
  # Read data
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)
  features <- read.delim(feature_file, header = FALSE, stringsAsFactors = FALSE)
  matrix <- readMM(matrix_file)
  
  # Add row and column names
  rownames(matrix) <- features$V1  # Gene names
  colnames(matrix) <- barcodes$V1  # Cell barcodes
  
  # Calculate UMI counts
  umi_counts <- colSums(matrix)
  
  # Generate barcode rank data
  barcode_rank_df <- data.frame(
    barcode = names(umi_counts),
    umi_count = umi_counts,
    sample = sample_name  # Add sample name
  )
  barcode_rank_df <- barcode_rank_df[order(barcode_rank_df$umi_count, decreasing = TRUE), ]
  barcode_rank_df$rank <- seq(1, nrow(barcode_rank_df))
  
  # Merge into total data
  barcode_rank_all <- rbind(barcode_rank_all, barcode_rank_df)
  
  # Save barcode_rank_df for each sample as CSV format
  barcode_rank_sample_csv <- paste0(outdir, "barcode_rank_", sample_name, ".csv")
  write.csv(barcode_rank_df, file = barcode_rank_sample_csv, row.names = FALSE)
}

# Save barcode_rank_all as CSV format
barcode_rank_all_csv <- paste0(outdir, "barcode_rank_all.csv")
write.csv(barcode_rank_all, file = barcode_rank_all_csv, row.names = FALSE)

# Add threshold column for each sample (for filtering convenience)
barcode_rank_all <- barcode_rank_all %>%
  mutate(threshold = thresholds[sample],
         highlight = umi_count > threshold)  # highlight = TRUE indicates highlighted portion


# Base plot: first draw grey curves for all samples (as background)
p <- ggplot() +
  geom_line(data = barcode_rank_all,
            aes(x = rank, y = umi_count, group = sample),
            color = grey_color, size = 1.0, alpha = 0.9)

# Then overlay the highlighted segments for each sample (umi_count > threshold)
# Filter by sample and plot using corresponding sample colors
for (s in unique(barcode_rank_all$sample)) {
  sub <- barcode_rank_all %>% filter(sample == s & highlight == TRUE)
  # Skip if no points meet the highlight condition for this sample
  if (nrow(sub) > 0) {
    p <- p + geom_line(data = sub,
                       aes(x = rank, y = umi_count, color = sample),
                       size = 1.5)
  }
}

# Get max rank value (for placing labels in the upper right corner)
max_rank <- max(barcode_rank_all$rank)

# Finally set colors, coordinate transformation, and theme (add three sample threshold lines)
barcode_rank_plot <- p +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title = "Barcode Rank",
       x = "",
       y = "UMI Counts",
       color = "Sample") +
  
  # Three threshold lines (10X_v3 = 5×10^3, InDrops-2 = 10^3, WF = 10^4)
  geom_hline(yintercept = 5e3,  linetype = "dashed", color = sample_colors["A"],     size = 0.8) +
  geom_hline(yintercept = 1e3,  linetype = "dashed", color = sample_colors["B"],  size = 0.8) +
  geom_hline(yintercept = 1e4,  linetype = "dashed", color = sample_colors["C"],         size = 0.8) +
  
  # ★★★ Place labels in the upper right corner of the threshold lines ★★★
  annotate("text", x = max_rank * 0.3, y = 5e3, 
           label = "5e3", 
           color = "black", hjust = 0, vjust = -0.5, size = 5) +
  
  annotate("text", x = max_rank * 0.3, y = 1e3, 
           label = "1e3", 
           color = "black", hjust = 0, vjust = -0.5, size = 5) +
  
  annotate("text", x = max_rank * 0.3, y = 1e4, 
           label = "1e4", 
           color = "black", hjust = 0, vjust = -0.5, size = 5) +
  
  scale_color_manual(values = sample_colors) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", size = 1.2),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save plots
barcode_rank_plot_svg <- paste0(outdir, "barcode_rank_plot_combined_highlighted.svg")
barcode_rank_plot_pdf <- paste0(outdir, "barcode_rank_plot_combined_highlighted.pdf")
ggsave(filename = barcode_rank_plot_pdf, plot = barcode_rank_plot, width = 7, height = 6)



