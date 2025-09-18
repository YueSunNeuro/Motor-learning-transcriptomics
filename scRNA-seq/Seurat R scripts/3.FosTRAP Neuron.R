################################################################################
# Title: Identify Neuronal Clusters & TRAP Cells – Fos-TRAP Motor Cortex
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
# This script labels tdTomato-expressing (TRAP) neurons and annotates neuronal
# clusters from Fos-TRAP single-cell RNA-seq data.
# Steps:
# 1. Load neuron-only Seurat object
# 2. Annotate neuronal clusters by cell type
# 3. Visualize marker genes across cell types
# 4. Normalize counts to CPM
# 5. Label tdTomato-positive (TRAP+) neurons
# 6. Visualize tdTomato expression by group and cell type
# 7. Export source data for figures
#
# Input:
# - Fos.Neuron.rds
#
# Output:
# - Fos.Neuron.TRAP.labeled.rds
# - Violin plots of tdTomato and marker genes
# - Source data CSV files for figures
#
# Note:
# tdTomato⁺ cells are identified by reads aligning to “wpre,” a unique 3′ sequence
# appended to the mm10 reference from the Ai9 (Rosa-CAG-LSL-tdTomato-WPRE) allele
# when building the custom Cell Ranger reference.
################################################################################

# Load required packages -------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(reticulate)
library(ggpubr)
library(chameleon) 

################################################################################
# 1. Load neuron-only Seurat object and inspect clusters
################################################################################
Fos.Neuron <- readRDS("Fos.Neuron.rds")
head(Fos.Neuron@meta.data)

Idents(Fos.Neuron) <- "seurat_clusters"
table(Idents(Fos.Neuron))

# UMAP plots
pdf("fig2.1.1.Fos.Neuron.UMAP.pdf", width=12, height=10)
DimPlot(Fos.Neuron, reduction="umap", pt.size=0.1, label=TRUE)
dev.off()

pdf("fig2.1.2.Fos.Neuron.UMAP.NOLABEL.pdf", width=12, height=10)
DimPlot(Fos.Neuron, reduction="umap", pt.size=0.1, label=FALSE)
dev.off()

pdf("fig2.1.3.Fos.Neuron.UMAP.split.pdf", width=15, height=5)
DimPlot(Fos.Neuron, reduction="umap", label=TRUE, split.by="groups")
dev.off()

################################################################################
# 2. Annotate neuronal clusters by cell type
################################################################################
Idents(Fos.Neuron) <- "seurat_clusters"
Fos.Neuron <- RenameIdents(Fos.Neuron,
                           `6`="L23_IT", `20`="L23_IT", `9`="L45_IT", `14`="L45_IT",
                           `23`="L56_IT-1", `32`="L56_IT-2", `10`="L56_IT-3", `52`="L6_IT",
                           `18`="L5_NP", `50`="L5_PT", `2`="L6_CT", `4`="L6_CT", `24`="L6b",
                           `40`="PV", `35`="Sst", `31`="Lamp5", `26`="Lamp5", `42`="Vip_Enpp2",
                           `38`="Vip_Htr3a", `43`="Deptor", `27`="Pthlh", `25`="Crh", `33`="Crh",
                           `47`="Crhbp")

Fos.Neuron$CellType <- Idents(Fos.Neuron)
table(Idents(Fos.Neuron))

# t-SNE plots of annotated cell types
pdf("fig2.1.3.Fos.Neuron.celltype.tsne.small.pdf", width = 6.5, height = 5)
DimPlot(object = Fos.Neuron, reduction = "tsne", pt.size= 0.5, label = T) + 
scale_color_chameleon(minimal_saturation = 10, minimal_lightness = 40, maximal_lightness = 80) 
dev.off()

pdf("fig2.4.4.Fos.Neuron.celltype.tsne.small.split.pdf", width = 15, height = 5)
DimPlot(object = Fos.Neuron, reduction = "tsne", pt.size= 0.5, label = TRUE, split.by = "groups")  + 
  scale_color_chameleon(minimal_saturation = 10, minimal_lightness = 40, maximal_lightness = 80)
dev.off()

################################################################################
# 3. Visualize marker genes by cell type
################################################################################
marker_genes <- c("Slc17a7","Slc30a3","Calb1","Rorb","Tnnc1","Deptor","Tesc","Bmp3","Oprk1",
                  "Rprm","Tshz2","Etv1","Lypd1","Nts","Foxp2","Ctgf",
                  "Gad1","Pvalb","Sst","Lamp5","Vip","Ndnf","Enpp2","Htr3a","Pthlh","Crh","Crhbp")

DefaultAssay(Fos.Neuron) <- "RNA"
Idents(Fos.Neuron) <- "CellType"

pdf("fig3.0.Fos.Neuron.all.VlnPlot-3.pdf", width=10, height=50)
VlnPlot(Fos.Neuron, features = marker_genes, pt.size = 0, ncol = 1)
dev.off()

################################################################################
# 4. Normalize RNA counts to CPM (RC method, counts per million)
################################################################################
DefaultAssay(Fos.Neuron) <- "RNA"
Fos.Neuron <- NormalizeData(Fos.Neuron, normalization.method = "RC", scale.factor = 1e6)

################################################################################
# 5. Label tdTomato-positive (TRAP+) neurons
################################################################################
# Visualize tdTom expression
pdf("fig2.6.1.Fos.Neuron.sub.tdTom.VlnPlot.RNA_CPM.pdf", width=10, height=4)
VlnPlot(Fos.Neuron, features = "wpre", pt.size = 0.1)
dev.off()

pdf("fig2.6.1.Fos.Neuron.sub.tdTom.VlnPlot.RNA_logCPM.pdf", width=10, height=4)
VlnPlot(Fos.Neuron, features = "wpre", pt.size = 0.1, log = TRUE)
dev.off()

# Subset TRAP+ neurons: CPM > 10
Fos.tdTom <- subset(Fos.Neuron, subset = wpre > 10)
Fos.tdTom$TRAP <- "P"  # tdTom positive
Fos.tdTom$rowname <- rownames(Fos.tdTom@meta.data)

# Subset TRAP- neurons
Fos.nontdTom <- subset(Fos.Neuron, subset = wpre < 10)
Fos.nontdTom$TRAP <- "N"  # tdTom negative
Fos.nontdTom$rowname <- rownames(Fos.nontdTom@meta.data)

# Merge TRAP labels back into main neuron metadata
Fos.tdTom.list <- merge(Fos.tdTom@meta.data, Fos.nontdTom@meta.data, all = TRUE)
rownames(Fos.tdTom.list) <- Fos.tdTom.list$rowname
Fos.tdTom.list <- subset(Fos.tdTom.list, select = "TRAP")

Fos.Neuron@meta.data <- merge(Fos.Neuron@meta.data, Fos.tdTom.list, by = "row.names", all = TRUE, sort = FALSE)
rownames(Fos.Neuron@meta.data) <- Fos.Neuron@meta.data$Row.names
Fos.Neuron@meta.data <- subset(Fos.Neuron@meta.data, select = -Row.names)

# Save updated Seurat object
saveRDS(Fos.Neuron, "Fos.Neuron.TRAP.labeled.rds")

################################################################################
# 6. Visualize tdTomato expression by group and cell type
################################################################################
Idents(Fos.Neuron) <- "groups"
pdf("fig3.0.Fos.Neuron.VlnPlot.tdTom.by.groups.CPM.pdf", width=6, height=6)
VlnPlot(Fos.Neuron, features = "wpre", pt.size = 0.1,
        cols = c("grey", "red"), group.by = "groups", split.by = "TRAP", log = TRUE)
dev.off()

Idents(Fos.Neuron)<- "CellType"
pdf("fig3.0.Fos.Neuron.VlnPlot.tdTom.by.CellTypes.CPM.pdf", width=10, height=4)
VlnPlot(Fos.Neuron, features = "wpre", pt.size = 0.1, split.by = "TRAP", log = T)
dev.off()

# TRAP Dimplot
Idents(Fos.Neuron)<- "TRAP"
pdf("fig2.5.3.Fos.Neuron.tdTom.FeaturePlot.split.small.pdf", width = 10.5, height = 4)
DimPlot(object = Fos.Neuron, reduction = "tsne", pt.size= 0.5, label = TRUE, cols=c("grey","red"),split.by = "groups", shuffle = T)
dev.off()
################################################################################
# 7. Export source data for figures
################################################################################
# t-SNE coordinates + expression for all marker genes
df <- as.data.frame(Embeddings(Fos.Neuron, "tsne"))
expr <- FetchData(Fos.Neuron, vars = marker_genes)
df <- cbind(df, expr)
df$CellType <- Idents(Fos.Neuron)
df$cluster <- Fos.Neuron$seurat_clusters
df$groups <- Fos.Neuron$groups
df$TRAP <- Fos.Neuron$TRAP
df$Cell_ID <- rownames(df)

# Add tdTom expression
expr_wpre <- FetchData(Fos.Neuron, vars = "wpre")
df <- cbind(df, expr_wpre)

write.csv(df, "Figure-tsne_source_data.csv", row.names = FALSE)

# VIP interneuron-specific source data
vip_genes <- c("Vip","Htr3a","Cck","Calb2","Npy","Chat")
df_vip <- as.data.frame(Embeddings(Fos.Neuron, "tsne"))
expr_vip <- FetchData(Fos.Neuron, vars = vip_genes)
df_vip <- cbind(df_vip, expr_vip)
df_vip$CellType <- Idents(Fos.Neuron)
df_vip$cluster <- Fos.Neuron$seurat_clusters
df_vip$groups <- Fos.Neuron$groups
df_vip$Cell_ID <- rownames(df_vip)

write.csv(df_vip, "Figure_VIP-IN_source_data.csv", row.names = FALSE)

################################################################################
# 8. Export number of TRAP vs non-TRAP cells across groups by cell type
################################################################################
# Build summary table
trap_summary <- Fos.Neuron@meta.data %>%
  dplyr::group_by(CellType, groups, TRAP) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  tidyr::unite(col = "Group_TRAP", groups, TRAP, sep = "-") %>%  # e.g. Ctrl-P, Ctrl-N
  tidyr::pivot_wider(names_from = "Group_TRAP", values_from = "n", values_fill = 0)

# Save as CSV
write.csv(trap_summary, "TRAP_nonTRAP_counts_byCellType.csv", row.names = FALSE)
