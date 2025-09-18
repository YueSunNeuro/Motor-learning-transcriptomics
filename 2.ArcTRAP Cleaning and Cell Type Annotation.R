################################################################################
# Title: Cluster Cleaning and Cell Type Annotation – Striatum
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
# This script performs post-integration processing of the Arc-TRAP motor cortex
# single-cell RNA-seq dataset. It includes:
# - Normalization and scaling of RNA counts
# - Cluster cleaning and renaming
# - Visualization (UMAP, t-SNE, Heatmaps, DotPlots, FeaturePlots)
# - Extraction of neurons and glia subsets
# - Generating source data for figures
#
# Inputs:
# - Arc.integrated.2.0.rds (integrated Seurat object from previous workflow)
#
# Outputs:
# - Arc.integrated.normalized.rds
# - Arc.clean.rds
# - UMAP/t-SNE plots
# - Heatmaps and DotPlots
# - Subset Seurat objects for neurons and glia
# - CSV files for source data and cluster-averaged expression
################################################################################

# Load required packages -------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(viridis)      # For heatmap color scaling
library(chameleon)    # Alternative color scales
library(RColorBrewer) # Color palettes for plots

################################################################################
# 1. Load integrated Seurat object and inspect
################################################################################
Arc.integrated <- readRDS("Arc.integrated.2.0.rds")
Arc.integrated
# 51149 features across 58378 samples within 3 assays 
# Active assay: RNA (24311 features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

head(Arc.integrated@meta.data)

# Inspect cell counts per cluster
table(Idents(Arc.integrated))

################################################################################
# 2. Normalize and scale RNA assay
################################################################################
DefaultAssay(Arc.integrated) <- "RNA"

# Normalize counts
Arc.integrated <- NormalizeData(Arc.integrated)

# Scale data and regress out mitochondrial content
Arc.integrated <- ScaleData(Arc.integrated, vars.to.regress = "percent.mt", verbose = TRUE)

# Save normalized object
saveRDS(Arc.integrated, "Arc.integrated.normalized.rds")

################################################################################
# 3. Build cluster tree and clean clusters
################################################################################
Arc.integrated <- BuildClusterTree(Arc.integrated, dims = 1:50)

# Visualize cluster hierarchy
pdf("fig.4.Arc.integrated.ClusterTree.pdf", width = 20, height = 10)
PlotClusterTree(Arc.integrated)
dev.off()

# Subset to retain only selected clusters
selected_clusters <- c("6","17","5","8","19","27","20","32","3","12","14",
                       "7","15","21","23","33","36","42","0","49","30","39",
                       "2","4","13","25","26","29","31","9","22","37","41",
                       "11","16","18","38","1","10","34","43","28")

Arc.clean <- subset(Arc.integrated, idents = selected_clusters)
Arc.clean
# An object of class Seurat 
# 51149 features across 55517 samples within 3 assays 
# Active assay: RNA (24311 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

table(Idents(Arc.clean))

# Save cleaned object
saveRDS(Arc.clean, "Arc.clean.rds")

################################################################################
# 4. Visualize cleaned clusters
################################################################################
# Rebuild cluster tree for cleaned object
Arc.clean <- BuildClusterTree(Arc.clean, dims = 1:50)
pdf("fig.4.Arc.clean.ClusterTree.pdf", width = 20, height = 10)
PlotClusterTree(Arc.clean)
dev.off()

# UMAP and t-SNE plots
pdf("fig2.1.Arc.clean.UMAP.NOLABEL.for.figure.pdf", width = 12, height = 10)
DimPlot(Arc.clean, reduction = "umap", pt.size = 0.2, label = FALSE)
dev.off()

pdf("fig2.1.Arc.clean.tsne.NOLABEL.for.figure.pdf", width = 12, height = 10)
DimPlot(Arc.clean, reduction = "tsne", pt.size = 0.2, label = FALSE)
dev.off()

################################################################################
# 5. Save gene list
################################################################################
all_genes <- rownames(Arc.clean)
write.table(all_genes, "Table0.Arc.clean.allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

Idents(Arc.clean) <- "groups"
ctrl.set <-subset(Arc.clean, idents = "ctrl")
Idents(ctrl.set) <- "batch"
ctrl.markers <- FindAllMarkers(ctrl.set, min.pct = 0.25, logfc.threshold = 0.4, only.pos = TRUE)
write.table(ctrl.markers, file="Table.Arc.clean.ctrl.batch.sensitive.genes.txt", sep="\t", quote=F, row.names=T, col.names=T)

################################################################################
# 6. Rename clusters as hierarchical sub-clusters for figure
################################################################################
Arc.clean<- RenameIdents(Arc.clean, `6` = "1.1",`17` = "1.2",`5` = "1.3",`8` = "1.4",
                         `19` = "1.5",`27` = "2.1",`20` = "2.2",`32` = "2.3",`3` = "3.1",
                         `12` = "3.2",`14` = "3.3",`7` = "3.4",`15` = "3.5",`21` = "3.6",
                         `23` = "4.1",`33` = "4.2",`36` = "4.3",`42` = "4.4",`2` = "5.1",
                         `26` = "5.2",`4` = "5.3",`25` = "5.4",`13` = "5.5",`29` = "5.6",
                         `31` = "5.7",`0` = "6.1",`49` = "6.2",`30` = "6.3",`39` = "6.4",
                         `34` = "7.1",`43` = "7.2",`1` = "7.3",`10` = "7.4",`41` = "8.1",
                         `37` = "8.2",`22` = "8.3",`9` = "8.4",`28` = "9",`18` = "10.1",
                         `16` = "10.2",`11` = "10.3",`38` = "11")

# Store as sub_clusters metadata
Arc.clean$clusters <- Idents(Arc.clean)
head(Arc.clean@meta.data)

# find markers for every cluster compared to all remaining cells, report only the positive ones
Arc.markers <- FindAllMarkers(Arc.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(Arc.markers, "Arc.markers.rds")

top5 <- Arc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

# Heatmap
pdf("fig3.0.Arc.clean.Heatmap.pdf", width=10, height=10)
DoHeatmap(Arc.clean, features = top5$gene) + scale_fill_viridis(limits=c(-1.5, 2.5))
dev.off()

# DotPlot
pdf("fig3.0.Arc.clean.Dotplot.pdf", width=10, height=12)
DotPlot(Arc.clean, features = c("Ccdc153", "Hexb","C1qb", "Dcn", "Vtn", "Rgs5", "Ly6c1", "Aldoc","Gja1","Top2a","Dcx","Mog","Olig1","Pdgfra",
                                "Chat","Th","Pthlh","Sst","Camk2n2","Drd2","Drd1"),  cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
dev.off()

################################################################################
# 7. Cell type annotation
################################################################################
Idents(Arc.clean) <- "seurat_clusters"
Arc.clean <- RenameIdents(Arc.clean, `6` = "dSPN",`17` = "dSPN",`5` = "dSPN",`8` = "dSPN",
                         `19` = "dSPN",`27` = "eSPN",`20` = "eSPN",`32` = "eSPN",`3` = "iSPN",
                         `12` = "iSPN",`14` = "iSPN",`7` = "iSPN",`15` = "iSPN",`21` = "iSPN",
                         `23` = "Int",`33` = "Int",`36` = "Int",`42` = "Int", `0` = "OPC", 
                         `49` = "OPC", `30` = "OL", `39` = "OL",`2` = "NG", `4` = "NG", 
                         `13` = "NG", `25` = "NG", `26` = "NG", `29` = "NG", `31` = "NG", 
                         `9` = "AS", `22` = "AS", `37` = "AS", `41` = "AS", `11` = "BV", 
                         `16` = "BV", `18` = "BV",`38` = "Fibroblast_like", `1` = "MG", 
                         `10` = "MG", `34` = "MG", `43` = "MG", `28`="Ependyma")

# Confirm cell type counts
table(Idents(Arc.clean))
# dSPN            eSPN            iSPN             Int             OPC              OL              NG              AS              BV Fibroblast_like              MG        Ependyma 
# 9272            2735           10479            2516            3825            1130           10129            3901            4415             387            5837             891 

Arc.clean$CellType <- Idents(Arc.clean)

Idents(Arc.clean) <- "seurat_clusters"
Arc.clean <- RenameIdents(Arc.clean, `6` = "dSPN",`17` = "dSPN",`5` = "dSPN",`8` = "dSPN",
                         `19` = "dSPN",`27` = "eSPN",`20` = "eSPN",`32` = "eSPN",`3` = "iSPN",
                         `12` = "iSPN",`14` = "iSPN",`7` = "iSPN",`15` = "iSPN",`21` = "iSPN",
                         `23`="Sst",`33`="PV",`36`="Th",`42`="Chat",`42` = "Int", `0` = "OPC", 
                         `49` = "OPC", `30` = "OL", `39` = "OL",`2` = "NG", `4` = "NG", 
                         `13` = "NG", `25` = "NG", `26` = "NG", `29` = "NG", `31` = "NG", 
                         `9` = "AS", `22` = "AS", `37` = "AS", `41` = "AS", `11` = "BV", 
                         `16` = "BV", `18` = "BV",`38` = "Fibroblast_like", `1` = "MG", 
                         `10` = "MG", `34` = "MG", `43` = "MG", `28`="Ependyma")
Arc.clean$SubType <- Idents(Arc.clean)

# Save labeled object
saveRDS(Arc.clean, "Arc.clean.celltype.rds")

################################################################################
# 8. Visualization: UMAP, t-SNE, DotPlots, FeaturePlots
################################################################################
# UMAP colored by cell type
pdf("fig2.1.Arc.clean.UMAP.celltype.pdf", width = 10, height = 8)
DimPlot(Arc.clean, label = TRUE, pt.size = 0.1) + scale_color_brewer(palette = "Paired")
dev.off()

# UMAP split by groups
pdf("fig2.4.Arc.clean.UMAP.celltype.split.pdf", width = 15, height = 5)
DimPlot(Arc.clean, label = TRUE, split.by = "groups") + scale_color_brewer(palette = "Paired")
dev.off()

# t-SNE for nFeature_RNA
pdf("fig2.Arc.clean.tsne.nFeature.FeaturePlot.pdf", width = 10, height = 8)
FeaturePlot(Arc.clean, features = "nFeature_RNA", reduction = "tsne",
            label = TRUE, min.cutoff = "q10", order = TRUE, cols = c("ivory2", "blue"), pt.size = 0.5)
dev.off()

# Marker gene feature plots
pdf("fig3.0.Arc.clean.Cell.FeaturePlot.small.pdf", width = 16, height = 14)
FeaturePlot(Arc.clean, features = c("Slc17a7","Gad2","Olig1","Aldh1l1","Ly6c1","Rgs5","Dcn","Hexb"),
            min.cutoff = "q10", order = TRUE, pt.size = 8, raster = TRUE, raster.dpi = c(3000, 3000))
dev.off()

################################################################################
# 9. Subset neurons and glia
################################################################################
# Neurons: PN + Int
Idents(Arc.clean) <- "CellType"
Arc.Neuron <- subset(Arc.clean, idents = c("PN","Int"))
saveRDS(Arc.Neuron, "Arc.Neuron.rds")

# Glia: OL, OPC, AS, MG
Arc.Glia <- subset(Arc.clean, idents = c("OL","OPC","AS","MG"))
saveRDS(Arc.Glia, "Arc.Glia.rds")

################################################################################
# 10. Generate source data for figures
################################################################################
# Extract UMAP coordinates and selected gene expression
genes <- c("Drd1","Drd2","Nkx2-1","Olig1","Aldh1l1","Hexb","Ly6c1","Rgs5","Dcn")
df <- as.data.frame(Embeddings(Arc.clean, "umap"))
expr <- FetchData(Arc.clean, vars = genes)
df <- cbind(df, expr)

df$CellType <- Idents(Arc.clean)
df$cluster <- Arc.clean$seurat_clusters
df$groups <- Arc.clean$groups
df$Cell_ID <- rownames(df)
write.csv(df, "Figure1_source_data.csv", row.names = FALSE)

# Cluster marker genes
Arc.markers <- FindAllMarkers(Arc.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(Arc.markers, "Arc.markers.by.CellType.rds")
write.csv(Arc.markers, "Arc.markers.by.CellType.csv")

# Top genes per cluster
top5 <- Arc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top20 <- Arc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Scaled expression for top genes
top_genes <- unique(top5$gene)
scaled_data <- Arc.clean[["RNA"]]@scale.data[top_genes, ]
scaled_df <- as.data.frame(t(scaled_data))
scaled_df$CellType <- Idents(Arc.clean)
scaled_df$Cell_ID <- rownames(scaled_df)
scaled_df <- scaled_df %>% select(Cell_ID, CellType, everything()) %>% arrange(CellType)
write.csv(scaled_df, "Ext.Fig.2_Heatmap_source_data.csv", row.names = FALSE)

# Cluster-averaged scaled expression
avg_expr <- AverageExpression(Arc.clean, features = unique(top20$gene),
                              group.by = "ident", slot = "scale.data")
avg_mat <- avg_expr$RNA
write.csv(avg_mat, "Ext.Fig.2_markergene_avg_expr_data.csv")
