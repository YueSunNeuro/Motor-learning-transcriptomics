################################################################################
# Title: Differential Expression Visualization
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
#   This script visualizes differential expression results by:
#     - Lollipop plots of up- and down-regulated DEGs by cell type
#     - Heatmaps showing DEG intersections across groups
#     - Volcano plots of DEGs by cluster and comparison
#     - Heatmaps of selected genes and average expression per group
#
#   Note: The example shown uses projection neuron cell types from the motor cortex
#         (Fos-TRAP). The same workflow applies to all cell types and to striatum
#         (Arc-TRAP) samples.
################################################################################

# ----------------------- 1. Load libraries -----------------------------------
library(ggridges)
library(ggplot2)
library(dplyr)
library(cowplot)
library(pheatmap)
library(viridis)
library(Seurat)

# ----------------------- 2. Lollipop plots for DEGs --------------------------
# Read DEG summary table
data <- read.csv("Fos.N.DEG.by.groups.adjP.csv", sep =",", check.names = FALSE)

# Ensure CellType order matches CSV
# Early vs Ctrl, upregulated
data$CellType <- factor(data$CellType, levels = data$CellType)

# Create the lollipop plot with labels
p1<- ggplot(data, aes(x = CellType, y = E_vs_C_up)) +
  geom_point(size = 10, color = "#FF9F9F") +  # Points
  geom_segment(aes(xend = CellType, yend = 0), color = "#FF9F9F", size = 1.5) +
  geom_text(aes(label = E_vs_C_up), vjust = 0.5, color = "white") +  # Labels
  labs(title = "Early vs Ctrl", x = "Cell Type", y = "Number of DEGs") +  # Labels
  ylim(0, 200) +
  coord_flip() +
  #scale_y_reverse()+
  theme_minimal()  # Theme

# Late vs Ctrl, upregulated
p2<- ggplot(data, aes(x = CellType, y = L_vs_C_up)) +
  geom_point(size = 10, color = "#FF002D") +  # Points
  geom_segment(aes(xend = CellType, yend = 0), color = "#FF002D", size = 1.5) +
  geom_text(aes(label = L_vs_C_up), vjust = 0.5, color = "white") +  # Labels
  labs(title = "Late vs Ctrl", x = "Cell Type", y = "Number of DEGs") +  # Labels
  ylim(0, 200) +
  coord_flip() +
  #scale_y_reverse()+
  theme_minimal()  # Theme

# Early vs Ctrl, downregulated
p3<- ggplot(data, aes(x = CellType, y = E_vs_C_down)) +
  geom_segment(aes(xend = CellType, yend = 0), color = "#6DD5FF", size = 1.5) +
  geom_point(size = 10, color = "#6DD5FF") +  # Points
  geom_text(aes(label = E_vs_C_down), vjust = 0.5, color = "white") +  # Labels
  labs(title = "Early vs Ctrl", x = "Cell Type", y = "Number of DEGs") +  # Labels
  scale_y_reverse(limits = c(200, 0)) +  # Reverse the y-axis and set limits
  coord_flip() +
  theme_minimal()  # Theme

# Late vs Ctrl, downregulated
p4<- ggplot(data, aes(x = CellType, y = L_vs_C_down)) +
  geom_segment(aes(xend = CellType, yend = 0), color = "#4273B9", size = 1.5) +
  geom_point(size = 10, color = "#4273B9") +  # Points
  geom_text(aes(label = L_vs_C_down), vjust = 0.5, color = "white") +  # Labels
  labs(title = "Late vs Ctrl", x = "Cell Type", y = "Number of DEGs") +  # Labels
  scale_y_reverse(limits = c(200, 0))+  
  coord_flip() +
  theme_minimal()  # Theme

plot_grid(p1, p2, p3, p4, labels = c('A', 'B','C', 'D'), ncol = 2)

# Save figure
ggsave("Fig.4.1.Fos.N.DEG.by.groups.adjP.Lollipop.pdf", width=8, height=10)

# ----------------------- 3. Heatmaps of intersecting DEGs --------------------
## Upregulated DEGs
# Load data from the CSV file
data <- read.csv("All_PN_gene.list.by.groups.up.Padj.csv", header = TRUE, 
                  stringsAsFactors = FALSE)

# Convert columns to lists of genes
genes_list <- lapply(data, as.character)
genes_list <- lapply(genes_list, function(x) x[x != ""])  # Remove empty elements
names(genes_list) <- colnames(data)  # Assign column names to the list

# Compute the pairwise intersections
intersection_matrix <- sapply(genes_list, function(x) {
  sapply(genes_list, function(y) length(intersect(x, y)))
})

# Mask the lower triangle
heatmap_data <- as.matrix(intersection_matrix)
heatmap_data[lower.tri(heatmap_data, diag = T)] <- NA  

# Plot heatmap
p1<-pheatmap(
  heatmap_data,
  display_numbers = TRUE,         # Show the actual numbers on the heatmap
  number_format = "%.0f",         # Format numbers as integers
  colorRampPalette(c("white","purple4"))(60),
  breaks = seq(0, 60, length.out = 60),  # Fixed color scale
  main = "Numbers of Intersecting DEGs",
  # fontsize = 10,
  na_col = "lightgrey",           # Set NA cells to white
  cluster_rows = FALSE,           # Disable row clustering
  cluster_cols = FALSE            # Disable column clustering
)
ggsave("fig3.Heatmap_Intersections.All_PN.up.Padj2.pdf", plot = p1, 
        width = 6.5, height = 6)

## Downregulated DEGs
data <- read.csv("All_PN_gene.list.by.groups.down.Padj.csv", header = TRUE, 
                  stringsAsFactors = FALSE)

# Convert columns to lists of genes
genes_list <- lapply(data, as.character)
genes_list <- lapply(genes_list, function(x) x[x != ""])  # Remove empty elements
names(genes_list) <- colnames(data)  # Assign column names to the list

# Compute the pairwise intersections
intersection_matrix <- sapply(genes_list, function(x) {
  sapply(genes_list, function(y) length(intersect(x, y)))
})

# Mask the lower triangle
heatmap_data <- as.matrix(intersection_matrix)
heatmap_data[lower.tri(heatmap_data, diag = F)] <- NA  # Mask lower triangle and diagonal

# Plot the heatmap
p2<-pheatmap(
  heatmap_data,
  display_numbers = TRUE,         # Show the actual numbers on the heatmap
  number_format = "%.0f",         # Format numbers as integers
  colorRampPalette(c("white","purple4"))(170),  # Color gradient color = magma(170), 
  breaks = seq(0, 170, length.out = 170),  # Fixed color scale
  main = "Numbers of Intersecting DEGs",
  # fontsize = 10,
  na_col = "lightgrey",           # Set NA cells to white
  cluster_rows = FALSE,           # Disable row clustering
  cluster_cols = FALSE            # Disable column clustering
)
ggsave("fig3.Heatmap_Intersections.All_PN.down.Padj.pdf", plot = p2, width = 6.5, height = 6)

# --------------------------------- 3. Volcano Plots ------------------------------
# Define the genes to exclude from plotting
excluded_genes <- c("CT010467.1", "Hbb-bs")

combined_data <- combined_data %>%
  filter(!(gene %in% excluded_genes))

# Add a combined significance score (|log2FC| × -log10 adjusted p-value)
combined_data <- combined_data %>%
  mutate(
    significance_score = abs(avg_log2FC) * -log10(p_val_adj)  # Combine fold change and p-value
  )

# Select top genes per cluster and comparison based on significance score
top_genes <- combined_data %>%
  group_by(cluster, Comparison) %>%
  arrange(desc(significance_score)) %>%
  slice_head(n = 20)  # Select top 10 genes per group

# Create the volcano plots
p <- ggplot(combined_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Regulation)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_grid(cluster~ Comparison, scales = "free_y") +  # separate panels per cluster & comparison
  theme_miminal () +
  labs(
    title = "DEGs in Motor cortex",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  # Significance and fold-change cutoffs
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Horizontal line at p = 0.05
  geom_vline(xintercept = c(-log2(1.2),log2(1.2)), linetype = "dashed", color = "black") +  # Vertical line at fold change = 0
  # Aesthetics 
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.ticks.x = element_line(size = 0.5)
  ) +
  xlim(-1.2, 1.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + # start y-axis at 0

  # Define custom colors for each regulation category
  scale_color_manual(values = c(
    "Early Up" = "#FF9F9F", 
    "Early Down" = "#6DD5FF", 
    "Late Up" = "#FF002D", 
    "Late Down" = "#4273B9"
  )) +
  
  # Label only the top 10 high and low expressed genes
  geom_text_repel(data = top_genes, aes(label = gene), color = "black", segment.size = 0.1,
                  size = 2, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.3)
ggsave("SmallVolcanoPlots.by.groups.padj'.pdf", plot = p, width=15, height=10)

# ----------- 4. Heatmaps of DEGs (average expression per group) ------------------------
Fos.Neuron <- readRDS("Fos.Neuron.TRAP.labeled.rds")
Fos.Neuron
# 49962 features across 20344 samples within 3 assays 
# Active assay: RNA (23590 features)
# 2 other assays present: SCT, Neuron
# 3 dimensional reductions calculated: pca, umap, tsne

# Normalize data back to log-normalization (required for correct fold-change)
Fos.Neuron <- NormalizeData(Fos.Neuron,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)

# Subset TRAP-labeled neurons and projection neurons
Idents(Fos.Neuron) <- "TRAP"
TRAP <- subset(Fos.Neuron, idents = "P")  # TRAP+ cells only
All_PN <- subset(TRAP,
                 idents = c("L23_IT-12","L23_IT-3","L45_IT","L56_IT-1","L56_IT-2","L6_CT"))
All_PN
# 49962 features across 2955 samples within 3 assays 
# Active assay: RNA (23590 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

# Create combined cell-group column
All_PN$Cell.groups <- paste(Idents(All_PN), All_PN$groups, sep = "_")

# Average expression heatmap of selected genes
Gene <- c("Cacna1a","Cacna2d1","Camk2a","Cnksr2","Eif4g1","Gria2","Nptn","Nrcam","Sv2a",
          "Napa","Atp6v1g2","Slc17a7","Slc24a2","Syp","Stxbp1","Syt1","Syn1","Cplx2",
          "Camk2b","Nrxn2","Grin1","Shank1","Syt7","Vamp2","Cplx1","Cplx2")

data <- subset(All_PN, features = Gene)
Idents(data) <- "CellType"

# Calculate average expression in each group
avg_expr <- AverageExpression(data, add.ident ="groups", assay = "RNA", slot = "data")
avg_expr_matrix <- avg_expr$RNA

write.csv(avg_expr_matrix, "Figure4c_source_data.csv", row.names = T)

p <- pheatmap(avg_expr_matrix, scale = "row", cluster_cols = F, cluster_rows = T,
             fontsize_row = 5, border_color = "grey", gaps_col= c(3,6,9,12,15),
             color = viridis(50)) 
ggsave("figx.0.Fos.PN.Neuron.Ave.expr.Heatmap.pdf", plot = p, width=6, height=10)
