################################################################################
# Title: Cluster Cleaning and Cell Type Annotation – Motor Cortex
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
# This script performs post-integration processing of the Fos-TRAP motor cortex
# single-cell RNA-seq dataset. It includes:
# - Normalization and scaling of RNA counts
# - Cluster cleaning and renaming
# - Visualization (UMAP, t-SNE, Heatmaps, DotPlots, FeaturePlots)
# - Extraction of neurons and glia subsets
# - Generating source data for figures
#
# Inputs:
# - Fos.integrated.2.0.rds (integrated Seurat object from previous workflow)
#
# Outputs:
# - Fos.integrated.normalized.rds
# - Fos.clean.rds
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
Fos.integrated <- readRDS("Fos.integrated.2.0.rds")
Fos.integrated
# 49962 features across 49779 samples within 3 assays 
# Active assay: RNA (23590 features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

head(Fos.integrated@meta.data)

# Inspect cell counts per cluster
table(Idents(Fos.integrated))

################################################################################
# 2. Normalize and scale RNA assay
################################################################################
DefaultAssay(Fos.integrated) <- "RNA"

# Normalize counts
Fos.integrated <- NormalizeData(Fos.integrated)

# Scale data and regress out mitochondrial content
Fos.integrated <- ScaleData(Fos.integrated, vars.to.regress = "percent.mt", verbose = TRUE)

# Save normalized object
saveRDS(Fos.integrated, "Fos.integrated.normalized.rds")

################################################################################
# 3. Build cluster tree and clean clusters
################################################################################
Fos.integrated <- BuildClusterTree(Fos.integrated, dims = 1:50)

# Visualize cluster hierarchy
pdf("fig.4.Fos.integrated.ClusterTree.pdf", width = 20, height = 10)
PlotClusterTree(Fos.integrated)
dev.off()

# Subset to retain only selected clusters
selected_clusters <- c("2","4","6","9","10","14","18","20","23","24","32","50","52","25","26","27",
                       "31","33","35","38","40","42","43","47","3","5","7","12","16","19","13",
                       "15","0","11","17","34","22","44","1","8","21","28","36","39")
Fos.clean <- subset(Fos.integrated, idents = selected_clusters)
Fos.clean
# An object of class Seurat 
# 49962 features across 47519 samples within 3 assays 
# Active assay: RNA (23590 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

table(Idents(Fos.clean))

# Save cleaned object
saveRDS(Fos.clean, "Fos.clean.rds")

################################################################################
# 4. Visualize cleaned clusters
################################################################################
# Rebuild cluster tree for cleaned object
Fos.clean <- BuildClusterTree(Fos.clean, dims = 1:50)
pdf("fig.4.Fos.clean.ClusterTree.pdf", width = 20, height = 10)
PlotClusterTree(Fos.clean)
dev.off()

# UMAP and t-SNE plots
pdf("fig2.1.Fos.clean.UMAP.NOLABEL.for.figure.pdf", width = 12, height = 10)
DimPlot(Fos.clean, reduction = "umap", pt.size = 0.2, label = FALSE)
dev.off()

pdf("fig2.1.Fos.clean.tsne.NOLABEL.for.figure.pdf", width = 12, height = 10)
DimPlot(Fos.clean, reduction = "tsne", pt.size = 0.2, label = FALSE)
dev.off()

################################################################################
# 5. Save gene list
################################################################################
all_genes <- rownames(Fos.clean)
write.table(all_genes, "Table0.Fos.clean.allgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

Idents(Fos.clean) <- "groups"
ctrl.set <-subset(Fos.clean, idents = "ctrl")
Idents(ctrl.set) <- "batch"
ctrl.markers <- FindAllMarkers(ctrl.set, min.pct = 0.25, logfc.threshold = 0.4, only.pos = TRUE)
write.table(ctrl.markers, file="Table.Fos.clean.ctrl.batch.sensitive.genes.txt", sep="\t", quote=F, row.names=T, col.names=T)

################################################################################
# 6. Rename clusters as hierarchical sub-clusters for figure
################################################################################
Fos.clean <- RenameIdents(Fos.clean,
                          `2` = "1.1.1",`4` = "1.1.2",`50` = "1.1.3",`24` = "1.2.1",
                          `18` = "1.3.1",`9` = "1.4.1",`14` = "1.4.2",`6` = "1.5.1",`20` = "1.5.2",
                          `10` = "1.6.1",`52` = "1.6.2",`23` = "1.7.1",`32` = "1.7.2",
                          `40` = "2.1.1",`35` = "2.1.2",`31` = "2.2.1",`26` = "2.2.2",
                          `42` = "2.3.1",`38` = "2.3.2",`43` = "2.3.3",`27` = "2.4.1",
                          `25` = "2.4.2",`33` = "2.4.3",`47` = "2.5.1",`22` = "3.1.1",
                          `44` = "3.1.2",`0` = "3.2.1",`17` = "3.2.2",`11` = "3.2.3",
                          `34` = "3.2.4",`39` = "4.1.1",`13` = "5.1.1",`15` = "5.1.2",
                          `8` = "6.1.1",`28` = "6.1.2",`21` = "6.2.1",`1` = "6.3.1",
                          `36` = "7.1.1",`3` = "8.1.1",`5` = "8.1.2",`7` = "8.1.3",
                          `19` = "8.1.4",`12` = "8.2.1",`16` = "8.2.2")

# Store as sub_clusters metadata
Fos.clean$sub_clusters <- Idents(Fos.clean)
head(Fos.clean@meta.data)

# find markers for every cluster compared to all remaining cells, report only the positive ones
Fos.markers <- FindAllMarkers(Fos.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(Fos.markers, "Fos.markers.rds")

top5 <- Fos.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

# Heatmap
pdf("fig3.0.Fos.clean.Heatmap.top5.NOLINE.pdf", width=10, height=10)
DoHeatmap(Fos.clean, features = top5$gene, group.by = "sub_clusters") + scale_fill_viridis(limits=c(-2.5, 2.5))
dev.off()

# DotPlot
pdf("fig3.0.Fos.clean.Dotplot.pdf", width=6.5, height=5.7)
DotPlot(Fos.clean, features = c("C1qb","Hexb","Dcn","Ly6c1","Rgs5","Vtn","Gja1","Aldoc","Top2a","Pdgfra","Olig1",
                                "Gad2","Gad1", "Slc17a8","Slc17a6", "Slc17a7"),  cols = c("blue", "red"), dot.scale = 6) + RotatedAxis()
dev.off()

################################################################################
# 7. Cell type annotation
################################################################################
Idents(Fos.clean) <- "seurat_clusters"
Fos.clean <- RenameIdents(Fos.clean,
                          `2` = "PN", `4` = "PN",`6` = "PN",`9` = "PN",`10` = "PN",`14` = "PN",
                          `18` = "PN",`20` = "PN",`23` = "PN",`24` = "PN",`32` = "PN",`50` = "PN",
                          `52` = "PN", `25` = "Int",`26` = "Int",`27` = "Int",`31` = "Int",
                          `33` = "Int",`35` = "Int",`38` = "Int",`40` = "Int",`42` = "Int",
                          `43` = "Int",`47` = "Int", `22` = "OL",`44` = "OL",`0` = "OPC",
                          `11` = "OPC",`17` = "OPC",`34` = "OPC",`39` = "NG",`13` = "AS",
                          `15` = "AS",`28` = "BV",`21` = "BV",`8` = "BV",`1` = "BV",
                          `36` = "Fibroblast_like",`3` = "MG",`5` = "MG",`7` = "MG",
                          `12` = "MG",`16` = "MG",`19` = "MG")

# Confirm cell type counts
table(Idents(Fos.clean))
#             PN             Int              OL             OPC              NG              AS              BV Fibroblast_like              MG 
#          16022            4322             964            6161             271            3311            5679             339           10450 

Fos.clean$CellType <- Idents(Fos.clean)
Idents(Fos.clean) <- "seurat_clusters"
Fos.clean <- RenameIdents(Fos.clean,
                          `6`="L23_IT", `20`="L23_IT", `9`="L45_IT", `14`="L45_IT",
                          `23`="L56_IT-1", `32`="L56_IT-2", `10`="L56_IT-3", `52`="L6_IT",
                          `18`="L5_NP", `50`="L5_PT", `2`="L6_CT", `4`="L6_CT", `24`="L6b",
                          `40`="PV", `35`="Sst", `31`="Lamp5", `26`="Lamp5", `42`="Vip_Enpp2",
                          `38`="Vip_Htr3a", `43`="Deptor", `27`="Pthlh", `25`="Crh", `33`="Crh",
                          `47`="Crhbp", `22` = "OL",`44` = "OL",`0` = "OPC",
                          `11` = "OPC",`17` = "OPC",`34` = "OPC",`39` = "NG",`13` = "AS",
                          `15` = "AS",`28` = "BV",`21` = "BV",`8` = "BV",`1` = "BV",
                          `36` = "Fibroblast_like",`3` = "MG",`5` = "MG",`7` = "MG",
                          `12` = "MG",`16` = "MG",`19` = "MG")

Fos.clean$SubType <- Idents(Fos.clean)

# Save labeled object
saveRDS(Fos.clean, "Fos.clean.celltype.rds")

################################################################################
# 8. Visualization: UMAP, t-SNE, DotPlots, FeaturePlots
################################################################################
# UMAP colored by cell type
pdf("fig2.1.Fos.clean.UMAP.celltype.pdf", width = 10, height = 8)
DimPlot(Fos.clean, label = TRUE, pt.size = 0.1) + scale_color_brewer(palette = "Paired")
dev.off()

# UMAP split by groups
pdf("fig2.4.Fos.clean.UMAP.celltype.split.pdf", width = 15, height = 5)
DimPlot(Fos.clean, label = TRUE, split.by = "groups") + scale_color_brewer(palette = "Paired")
dev.off()

# t-SNE for nFeature_RNA
pdf("fig2.Fos.clean.tsne.nFeature.FeaturePlot.pdf", width = 10, height = 8)
FeaturePlot(Fos.clean, features = "nFeature_RNA", reduction = "tsne",
            label = TRUE, min.cutoff = "q10", order = TRUE, cols = c("ivory2", "blue"), pt.size = 0.5)
dev.off()

# Marker gene feature plots
pdf("fig3.0.Fos.clean.Cell.FeaturePlot.small.pdf", width = 16, height = 14)
FeaturePlot(Fos.clean, features = c("Slc17a7","Gad2","Olig1","Aldh1l1","Ly6c1","Rgs5","Dcn","Hexb"),
            min.cutoff = "q10", order = TRUE, pt.size = 8, raster = TRUE, raster.dpi = c(3000, 3000))
dev.off()

################################################################################
# 9. Subset neurons and glia
################################################################################
# Neurons: PN + Int
Fos.Neuron <- subset(Fos.clean, idents = c("PN","Int"))
saveRDS(Fos.Neuron, "Fos.Neuron.rds")

# Glia: OL, OPC, AS, MG
Fos.Glia <- subset(Fos.clean, idents = c("OL","OPC","AS","MG"))
saveRDS(Fos.Glia, "Fos.Glia.rds")

################################################################################
# 10. Generate source data for figures
################################################################################
# Extract UMAP coordinates and selected gene expression
genes <- c("Rbfox3","Slc17a7","Gad2","Aldh1l1","Olig1","Hexb","Ly6c1","Rgs5","Dcn")
df <- as.data.frame(Embeddings(Fos.clean, "umap"))
expr <- FetchData(Fos.clean, vars = genes)
df <- cbind(df, expr)

df$CellType <- Idents(Fos.clean)
df$cluster <- Fos.clean$seurat_clusters
df$groups <- Fos.clean$groups
df$Cell_ID <- rownames(df)
write.csv(df, "Figure1_source_data.csv", row.names = FALSE)

# Cluster marker genes
Fos.markers <- FindAllMarkers(Fos.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(Fos.markers, "Fos.markers.by.CellType.rds")
write.csv(Fos.markers, "Fos.markers.by.CellType.csv")

# Top genes per cluster
top5 <- Fos.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top20 <- Fos.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Scaled expression for top genes
top_genes <- unique(top5$gene)
scaled_data <- Fos.clean[["RNA"]]@scale.data[top_genes, ]
scaled_df <- as.data.frame(t(scaled_data))
scaled_df$CellType <- Idents(Fos.clean)
scaled_df$Cell_ID <- rownames(scaled_df)
scaled_df <- scaled_df %>% select(Cell_ID, CellType, everything()) %>% arrange(CellType)
write.csv(scaled_df, "Ext.Fig.2_Heatmap_source_data.csv", row.names = FALSE)

# Cluster-averaged scaled expression
avg_expr <- AverageExpression(Fos.clean, features = unique(top20$gene),
                              group.by = "ident", slot = "scale.data")
avg_mat <- avg_expr$RNA
write.csv(avg_mat, "Ext.Fig.2_markergene_avg_expr_data.csv")
