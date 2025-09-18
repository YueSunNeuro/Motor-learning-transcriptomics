################################################################################
# Title: Identify Neuronal Clusters & TRAP Cells – Arc-TRAP Striatum
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
# This script labels tdTomato-expressing (TRAP) neurons and annotates neuronal
# clusters from Arc-TRAP single-cell RNA-seq data.
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
# - Arc.Neuron.rds
#
# Output:
# - Arc.Neuron.TRAP.labeled.rds
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
Arc.Neuron <- readRDS("Arc.Neuron.rds")
head(Arc.Neuron@meta.data)

Idents(Arc.Neuron) <- "seurat_clusters"
table(Idents(Arc.Neuron))

# UMAP plots
pdf("fig2.1.1.Arc.Neuron.UMAP.pdf", width=12, height=10)
DimPlot(Arc.Neuron, reduction="umap", pt.size=0.1, label=TRUE)
dev.off()

pdf("fig2.1.2.Arc.Neuron.UMAP.NOLABEL.pdf", width=12, height=10)
DimPlot(Arc.Neuron, reduction="umap", pt.size=0.1, label=FALSE)
dev.off()

pdf("fig2.1.3.Arc.Neuron.UMAP.split.pdf", width=15, height=5)
DimPlot(Arc.Neuron, reduction="umap", label=TRUE, split.by="groups")
dev.off()

################################################################################
# 2. Annotate neuronal clusters by cell type
################################################################################
Idents(Arc.Neuron) <- "seurat_clusters"
Arc.Neuron <-RenameIdents(Arc.Neuron, `6`="dSPN",`17`="dSPN",`5`="dSPN",`8`="dSPN",`19`="dSPN",
                          `27`="eSPN",`20`="eSPN",`32`="eSPN",
                          `3`="iSPN",`12`="iSPN",`14`="iSPN",`7`="iSPN",`15`="iSPN",`21`="iSPN",
                          `23`="Sst",`33`="PV",`36`="Th",`42`="Chat")

Arc.Neuron$CellType <- Idents(Arc.Neuron)
table(Idents(Arc.Neuron))

# t-SNE plots of annotated cell types
pdf("fig2.1.3.Arc.Neuron.celltype.tsne.small.pdf", width = 6.5, height = 5)
DimPlot(object = Arc.Neuron, reduction = "tsne", pt.size= 0.5, label = T) + 
scale_color_chameleon(minimal_saturation = 10, minimal_lightness = 20, maximal_lightness = 80) 
dev.off()

pdf("fig2.4.4.Arc.Neuron.celltype.tsne.small.split.pdf", width = 15, height = 5)
DimPlot(object = Arc.Neuron, reduction = "tsne", pt.size= 0.5, label = TRUE, split.by = "groups")  + 
scale_color_chameleon(minimal_saturation = 10, minimal_lightness = 20, maximal_lightness = 80) 
dev.off()

################################################################################
# 3. Visualize marker genes by cell type
################################################################################
marker_genes <- c("Tac1","Drd1","Pdyn","Id4","Rn7sk","Wfs1","Lypd1","Tshz1","Otof",
                  "Penk","Drd2","Adora2a","Arc","Me2","Camk2n2","Sst","Npy","Pvalb",
                  "Pthlh","Kit","Th","Htr3a","Chat","Slc17a8")

DefaultAssay(Arc.Neuron) <- "RNA"
Idents(Arc.Neuron) <- "CellType"

pdf("fig3.0.Arc.Neuron.all.VlnPlot-3.pdf", width=10, height=50)
VlnPlot(Arc.Neuron, features = marker_genes, pt.size = 0, ncol = 1)
dev.off()

################################################################################
# 4. Normalize RNA counts to CPM (RC method, counts per million)
################################################################################
DefaultAssay(Arc.Neuron) <- "RNA"
Arc.Neuron <- NormalizeData(Arc.Neuron, normalization.method = "RC", scale.factor = 1e6)

################################################################################
# 5. Label tdTomato-positive (TRAP+) neurons
################################################################################
# Visualize tdTom expression
pdf("fig2.6.1.Arc.Neuron.sub.tdTom.VlnPlot.RNA_CPM.pdf", width=10, height=4)
VlnPlot(Arc.Neuron, features = "wpre", pt.size = 0.1)
dev.off()

pdf("fig2.6.1.Arc.Neuron.sub.tdTom.VlnPlot.RNA_logCPM.pdf", width=10, height=4)
VlnPlot(Arc.Neuron, features = "wpre", pt.size = 0.1, log = TRUE)
dev.off()

# Subset TRAP+ neurons: CPM > 10
Arc.tdTom <- subset(Arc.Neuron, subset = wpre > 10)
Arc.tdTom$TRAP <- "P"  # tdTom positive
Arc.tdTom$rowname <- rownames(Arc.tdTom@meta.data)

# Subset TRAP- neurons
Arc.nontdTom <- subset(Arc.Neuron, subset = wpre < 10)
Arc.nontdTom$TRAP <- "N"  # tdTom negative
Arc.nontdTom$rowname <- rownames(Arc.nontdTom@meta.data)

# Merge TRAP labels back into main neuron metadata
Arc.tdTom.list <- merge(Arc.tdTom@meta.data, Arc.nontdTom@meta.data, all = TRUE)
rownames(Arc.tdTom.list) <- Arc.tdTom.list$rowname
Arc.tdTom.list <- subset(Arc.tdTom.list, select = "TRAP")

Arc.Neuron@meta.data <- merge(Arc.Neuron@meta.data, Arc.tdTom.list, by = "row.names", all = TRUE, sort = FALSE)
rownames(Arc.Neuron@meta.data) <- Arc.Neuron@meta.data$Row.names
Arc.Neuron@meta.data <- subset(Arc.Neuron@meta.data, select = -Row.names)

# Save updated Seurat object
saveRDS(Arc.Neuron, "Arc.Neuron.TRAP.labeled.rds")

################################################################################
# 6. Visualize tdTomato expression by group and cell type
################################################################################
Idents(Arc.Neuron) <- "groups"
pdf("fig3.0.Arc.Neuron.VlnPlot.tdTom.by.groups.CPM.pdf", width=6, height=6)
VlnPlot(Arc.Neuron, features = "wpre", pt.size = 0.1,
        cols = c("grey", "red"), group.by = "groups", split.by = "TRAP", log = TRUE)
dev.off()

Idents(Arc.Neuron)<- "CellType"
pdf("fig3.0.Arc.Neuron.VlnPlot.tdTom.by.CellTypes.CPM.pdf", width=10, height=4)
VlnPlot(Arc.Neuron, features = "wpre", pt.size = 0.1, split.by = "TRAP", log = T)
dev.off()

# TRAP Dimplot
Idents(Arc.Neuron)<- "TRAP"
pdf("fig2.5.3.Arc.Neuron.tdTom.FeaturePlot.split.small.pdf", width = 10.5, height = 4)
DimPlot(object = Arc.Neuron, reduction = "tsne", pt.size= 0.5, label = TRUE, cols=c("grey","red"),split.by = "groups", shuffle = T)
dev.off()
################################################################################
# 7. Export source data for figures
################################################################################
# t-SNE coordinates + expression for all marker genes
df <- as.data.frame(Embeddings(Arc.Neuron, "tsne"))
expr <- FetchData(Arc.Neuron, vars = marker_genes)
df <- cbind(df, expr)
df$CellType <- Idents(Arc.Neuron)
df$cluster <- Arc.Neuron$seurat_clusters
df$groups <- Arc.Neuron$groups
df$TRAP <- Arc.Neuron$TRAP
df$Cell_ID <- rownames(df)

# Add tdTom expression
expr_wpre <- FetchData(Arc.Neuron, vars = "wpre")
df <- cbind(df, expr_wpre)

write.csv(df, "Figure-tsne_source_data.csv", row.names = FALSE)

################################################################################
# 8. Export number of TRAP vs non-TRAP cells across groups by cell type
################################################################################
# Build summary table
trap_summary <- Arc.Neuron@meta.data %>%
  dplyr::group_by(CellType, groups, TRAP) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  tidyr::unite(col = "Group_TRAP", groups, TRAP, sep = "-") %>%  # e.g. Ctrl-P, Ctrl-N
  tidyr::pivot_wider(names_from = "Group_TRAP", values_from = "n", values_fill = 0)

# Save as CSV
write.csv(trap_summary, "TRAP_nonTRAP_counts_byCellType.csv", row.names = FALSE)
