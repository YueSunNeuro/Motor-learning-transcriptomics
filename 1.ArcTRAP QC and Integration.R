################################################################################
# Title: Single-cell RNA-seq Processing and Integration – Striatum (Arc-TRAP)
# Author: Yue Sun
# Date: 2025-09-12
#
# Description:
# This script processes 10x Genomics single-cell RNA-seq data from multiple
# mouse striatum samples (control, early learning, late learning).
# It creates Seurat objects, performs quality control, merges and normalizes
# the data, and integrates the datasets using SCTransform.
#
# Inputs:
# - Raw 10x Genomics matrix files (available from GEO accession GSM9059721)
# - Striatum data are processed and integrated here.
# - The same workflow can be applied to motor cortex (Arc-TRAP) data.
#
# Outputs:
# - Individual and merged Seurat objects saved as .rds files
# - QC violin plots saved as .pdf files
# - Integrated Seurat object with PCA/UMAP clustering
################################################################################

# Load required packages -------------------------------------------------------
library(Seurat)       
library(ggplot2)      
library(cowplot)      # multi-panel plots
library(grid)         # low-level grid graphics

################################################################################
# 1. Define input directories (10x Genomics filtered_feature_bc_matrix folders)
################################################################################

# Control samples (update with your actual data paths)
DIR_Arc_Ctrl_1 = "path/to/your/data/yue_YS05-1/outs/filtered_feature_bc_matrix"
DIR_Arc_Ctrl_2 = "path/to/your/data/yue_YS05-2/outs/filtered_feature_bc_matrix"
DIR_Arc_Ctrl_3 = "path/to/your/data/yue_YS08-1/outs/filtered_feature_bc_matrix"
DIR_Arc_Ctrl_4 = "path/to/your/data/yue_YS08-2/outs/filtered_feature_bc_matrix"

# Early learning samples
DIR_Arc_Early_C1 = "path/to/your/data/yue_YS09-1/outs/filtered_feature_bc_matrix"
DIR_Arc_Early_I1 = "path/to/your/data/yue_YS09-2/outs/filtered_feature_bc_matrix"
DIR_Arc_Early_C2 = "path/to/your/data/yue_YS11-1/outs/filtered_feature_bc_matrix"
DIR_Arc_Early_I2 = "path/to/your/data/yue_YS11-2/outs/filtered_feature_bc_matrix"

# Late learning samples
DIR_Arc_trained_C1 = "path/to/your/data/yue_YS04-1/outs/filtered_feature_bc_matrix"
DIR_Arc_trained_I1 = "path/to/your/data/yue_YS04-2/outs/filtered_feature_bc_matrix"
DIR_Arc_trained_C2 = "path/to/your/data/yue_YS13-1/outs/filtered_feature_bc_matrix"
DIR_Arc_trained_I2 = "path/to/your/data/yue_YS13-2/outs/filtered_feature_bc_matrix"

################################################################################
# 2. Read 10X data -------------------------------------------------------------
# Each call to Read10X returns a sparse matrix of counts.
################################################################################
ctrl_1.data <- Read10X(data.dir = DIR_Arc_Ctrl_1)
ctrl_2.data <- Read10X(data.dir = DIR_Arc_Ctrl_2)
ctrl_3.data <- Read10X(data.dir = DIR_Arc_Ctrl_3)
ctrl_4.data <- Read10X(data.dir = DIR_Arc_Ctrl_4)

Early_C1.data <- Read10X(data.dir = DIR_Arc_Early_C1)
Early_I1.data <- Read10X(data.dir = DIR_Arc_Early_I1)
Early_C2.data <- Read10X(data.dir = DIR_Arc_Early_C2)
Early_I2.data <- Read10X(data.dir = DIR_Arc_Early_I2)

late_I1.data <- Read10X(data.dir = DIR_Arc_trained_I1)
late_C1.data <- Read10X(data.dir = DIR_Arc_trained_C1)
late_I2.data <- Read10X(data.dir = DIR_Arc_trained_I2)
late_C2.data <- Read10X(data.dir = DIR_Arc_trained_C2)

################################################################################
# 3. Create Seurat objects and annotate metadata --------------------------------
# Each object represents one sample. Add group, lateral side, batch information.
################################################################################

## Control samples
ctrl_1 <- CreateSeuratObject(counts = ctrl_1.data, project = "C1", min.cells = 5)
ctrl_2 <- CreateSeuratObject(counts = ctrl_2.data, project = "C2", min.cells = 5)
ctrl_3 <- CreateSeuratObject(counts = ctrl_3.data, project = "C3", min.cells = 5)
ctrl_4 <- CreateSeuratObject(counts = ctrl_4.data, project = "C4", min.cells = 5)

# Annotate group (condition), lateral side, and batch
ctrl_1$groups <-"ctrl"; ctrl_1$lateral <-"ctrl_L"; ctrl_1$batch <-"ctrl_B1"
ctrl_2$groups <-"ctrl"; ctrl_2$lateral <-"ctrl_R"; ctrl_2$batch <-"ctrl_B1"
ctrl_3$groups <-"ctrl"; ctrl_3$lateral <-"ctrl_R"; ctrl_3$batch <-"ctrl_B2"
ctrl_4$groups <-"ctrl"; ctrl_4$lateral <-"ctrl_L"; ctrl_4$batch <-"ctrl_B2"

# Calculate mitochondrial (%) and ribosomal (%) content per cell
ctrl_1[["percent.mt"]] <- PercentageFeatureSet(object = ctrl_1, pattern = "^mt-")
ctrl_2[["percent.mt"]] <- PercentageFeatureSet(object = ctrl_2, pattern = "^mt-")
ctrl_3[["percent.mt"]] <- PercentageFeatureSet(object = ctrl_3, pattern = "^mt-")
ctrl_4[["percent.mt"]] <- PercentageFeatureSet(object = ctrl_4, pattern = "^mt-")

ctrl_1[["percent.ribo"]] <- PercentageFeatureSet(object = ctrl_1, pattern = "^Rna")
ctrl_2[["percent.ribo"]] <- PercentageFeatureSet(object = ctrl_2, pattern = "^Rna")
ctrl_3[["percent.ribo"]] <- PercentageFeatureSet(object = ctrl_3, pattern = "^Rna")
ctrl_4[["percent.ribo"]] <- PercentageFeatureSet(object = ctrl_4, pattern = "^Rna")

# Repeat the same block for early_I1, early_C1, early_I2, early_C2 and
# for late_I1, late_C1, late_I2, late_C2 to ensure consistent metadata.

## Early samples
early_I1 <- CreateSeuratObject(counts = Early_I1.data, project = "EI1", min.cells = 5)
early_C1 <- CreateSeuratObject(counts = Early_C1.data, project = "EC1", min.cells = 5)
early_I2 <- CreateSeuratObject(counts = Early_I2.data, project = "EI2", min.cells = 5)
early_C2 <- CreateSeuratObject(counts = Early_C2.data, project = "EC2", min.cells = 5)

# Annotate group (condition), lateral side, and batch
early_I1$groups <-"early"; early_I1$lateral <-"early_I"; early_I1$batch <-"early_B1"
early_C1$groups <-"early"; early_C1$lateral <-"early_C"; early_C1$batch <-"early_B1"
early_I2$groups <-"early"; early_I2$lateral <-"early_I"; early_I2$batch <-"early_B2"
early_C2$groups <-"early"; early_C2$lateral <-"early_C"; early_C2$batch <-"early_B2"

# Calculate mitochondrial (%) and ribosomal (%) content per cell
early_I1[["percent.mt"]] <- PercentageFeatureSet(object = early_I1, pattern = "^mt-")
early_C1[["percent.mt"]] <- PercentageFeatureSet(object = early_C1, pattern = "^mt-")
early_I1[["percent.ribo"]] <- PercentageFeatureSet(object = early_I1, pattern = "^Rna")
early_C1[["percent.ribo"]] <- PercentageFeatureSet(object = early_C1, pattern = "^Rna")

early_I2[["percent.mt"]] <- PercentageFeatureSet(object = early_I2, pattern = "^mt-")
early_C2[["percent.mt"]] <- PercentageFeatureSet(object = early_C2, pattern = "^mt-")
early_I2[["percent.ribo"]] <- PercentageFeatureSet(object = early_I2, pattern = "^Rna")
early_C2[["percent.ribo"]] <- PercentageFeatureSet(object = early_C2, pattern = "^Rna")

## Late samples
late_I1 <- CreateSeuratObject(counts = late_I1.data, project = "LI1", min.cells = 5)
late_C1 <- CreateSeuratObject(counts = late_C1.data, project = "LC1", min.cells = 5)
late_I2 <- CreateSeuratObject(counts = late_I2.data, project = "LI2", min.cells = 5)
late_C2 <- CreateSeuratObject(counts = late_C2.data, project = "LC2", min.cells = 5)

# Annotate group (condition), lateral side, and batch
late_I1$groups <-"late"; late_I1$lateral <-"late_I"; late_I1$batch <-"late_B1"
late_C1$groups <-"late"; late_C1$lateral <-"late_C"; late_C1$batch <-"late_B1"
late_I2$groups <-"late"; late_I2$lateral <-"late_I"; late_I2$batch <-"late_B2"
late_C2$groups <-"late"; late_C2$lateral <-"late_C"; late_C2$batch <-"late_B2"

# Calculate mitochondrial (%) and ribosomal (%) content per cell
late_I1[["percent.mt"]] <- PercentageFeatureSet(object = late_I1, pattern = "^mt-")
late_C1[["percent.mt"]] <- PercentageFeatureSet(object = late_C1, pattern = "^mt-")
late_I1[["percent.ribo"]] <- PercentageFeatureSet(object = late_I1, pattern = "^Rna")
late_C1[["percent.ribo"]] <- PercentageFeatureSet(object = late_C1, pattern = "^Rna")

late_I2[["percent.mt"]] <- PercentageFeatureSet(object = late_I2, pattern = "^mt-")
late_C2[["percent.mt"]] <- PercentageFeatureSet(object = late_C2, pattern = "^mt-")
late_I2[["percent.ribo"]] <- PercentageFeatureSet(object = late_I2, pattern = "^Rna")
late_C2[["percent.ribo"]] <- PercentageFeatureSet(object = late_C2, pattern = "^Rna")

################################################################################
# 4. Merge Seurat objects per condition ----------------------------------------
# Combine replicates for each condition separately.
################################################################################
ctrl.merge <- merge(ctrl_1, y = c(ctrl_2, ctrl_3),
                    add.cell.ids = c("CTRL_1","CTRL_2","CTRL_3"),
                    project = "ctrl.merge")
saveRDS(ctrl.merge, "ctrl_merge.rds")

early.merge <- merge(early_I1, y = c(early_C1, early_I2, early_C2),
                     add.cell.ids = c("early_I1","early_C1","early_I2","early_C2"),
                     project = "early.merge")
saveRDS(early.merge, "early_merge.rds")

late.merge <- merge(late_I1, y = c(late_C1, late_I2, late_C2),
                    add.cell.ids = c("late_I1","late_C1","late_I2","late_C2"),
                    project = "late.merge")
saveRDS(late.merge, "late_merge.rds")

################################################################################
# 5. Merge all groups together -------------------------------------------------
# Combine control, early, and late samples into one Seurat object.
################################################################################
Arc.merge <- merge(ctrl_1, y = c(ctrl_2,ctrl_3,ctrl_4,
                   early_I1, early_C1,early_I2, early_C2,
                   late_I1, late_C1,late_I2, late_C2), 
                   add.cell.ids = c("CTRL_1","CTRL_2", "CTRL_3","CTRL_4", 
                    "early_I1", "early_C1", "early_I2", "early_C2", 
                    "late_I1", "late_C1", "late_I2", "late_C2"), 
                   project = "Arc_merge")
saveRDS(Arc.merge, "Arc_merge.rds")

################################################################################
# 6. Quality Control Plots -----------------------------------------------------
# Visualize nCount_RNA, nFeature_RNA, percent.mt, percent.ribo across samples.
################################################################################
Idents(Arc.merge) <- "orig.ident"
VlnPlot(Arc.merge, features = "nCount_RNA", pt.size = 0)
ggsave("fig0.1.Arc.nCount.vlnplot.pdf", width=12, height=5)
# Repeat with Idents set to "groups","lateral","batch" and for
# nFeature_RNA, percent.mt, percent.ribo.

################################################################################
# 7. Filter cells --------------------------------------------------------------
# Example thresholds: >500 detected genes and <5% mitochondrial reads.
################################################################################
Arc.merge <- subset(Arc.merge, subset = nFeature_RNA > 500 & percent.mt < 5)
Arc.merge
# An object of class Seurat 
# 24311 features across 58378 samples within 1 assay 
# Active assay: RNA (24311 features)
################################################################################
# 8. SCTransform normalization per batch ---------------------------------------
# Normalize each batch separately before integration.
################################################################################
Arc.list <- SplitObject(Arc.merge, split.by = "batch")

for (i in 1:length(Arc.list)) {
  Arc.list[[i]] <- SCTransform(Arc.list[[i]],
                               vars.to.regress = "percent.mt",
                               return.only.var.genes = FALSE,
                               verbose = TRUE)
}
saveRDS(Arc.list, "Arc.list.rds")

################################################################################
# 9. Integration using SCTransform ---------------------------------------------
# Identify features for integration, prepare data, find anchors, then integrate.
################################################################################
Arc.features <- SelectIntegrationFeatures(Arc.list, nfeatures = 3000)
Arc.list <- PrepSCTIntegration(Arc.list, anchor.features = Arc.features)
Arc.anchors <- FindIntegrationAnchors(Arc.list,
                                      normalization.method = "SCT",
                                      anchor.features = Arc.features)
saveRDS(Arc.anchors, "Arc.anchors.rds")

Arc.integrated <- IntegrateData(anchorset = Arc.anchors, normalization.method = "SCT")

################################################################################
# 10. Dimensionality Reduction & Clustering ------------------------------------
# Run PCA, choose dimensions, cluster cells, and run UMAP/t-SNE.
################################################################################
Arc.integrated <- RunPCA(Arc.integrated, npcs = 50)

pdf("fig1.1.Arc.integrated.pca.pdf",width=5,height=6)
VizDimLoadings(Arc.integrated, dims=1:2)
dev.off()

# JackStraw and Elbow plots help choose number of PCs
Arc.integrated<- JackStraw(Arc.integrated, num.replicate = 50)
Arc.integrated <- ScoreJackStraw(Arc.integrated, dims = 1:50)

pdf("fig1.3.Arc.integrated.JackStraw.pca.pdf",width=6, height=4)
JackStrawPlot(Arc.integrated, dims = 1:50)
dev.off()

pdf("fig1.4.Arc.integrated.ElbowPlot.pca.pdf",width=6, height=4)
ElbowPlot(Arc.integrated, ndims = 100)
dev.off()

Arc.integrated <- FindNeighbors(Arc.integrated, dims = 1:50)
Arc.integrated <- FindClusters(Arc.integrated, resolution = 2)
Arc.integrated <- RunUMAP(Arc.integrated, dims = 1:50)
Arc.integrated <- RunTSNE(Arc.integrated, dims = 1:50)

# Save UMAP plots
pdf("fig2.0.Arc.UMAP.pdf", width = 12, height = 10)
DimPlot(object = Arc.integrated, reduction = "umap", pt.size= 0.2, label = TRUE)
dev.off()

pdf("fig2.1.Arc.UMAP.split.pdf", width = 12, height = 5)
DimPlot(object = Arc.integrated, label = TRUE,split.by = "groups")
dev.off()

pdf("fig2.2.Arc.UMAP.subgroups.pdf", width = 8, height = 6)
DimPlot(object = Arc.integrated, label = FALSE,group.by = "orig.ident")
dev.off()

################################################################################
# 11. Save gene list and integrated object -------------------------------------
################################################################################
write.table(rownames(Arc.integrated), "Table0.Arc.integrated.allgenes.txt",
            sep="\t", quote=FALSE, row.names=TRUE, col.names=FALSE)

saveRDS(Arc.integrated, "Arc.integrated.2.0.rds")

################################################################################
# 12. Identify cell markers per cluster ----------------------------------------
# Loop through all clusters and find conserved markers across groups.
################################################################################
dim_array_clusters =  dim(table(Idents(Arc.integrated))) - 1 
LIST.CLUSTERS.and.MARKERS <- list()

DefaultAssay(object = Arc.integrated) <- "SCT"

for (i in 0:dim_array_clusters)
{ 
  LIST.CLUSTERS.and.MARKERS[[i+1]] <-FindConservedMarkers(Arc.integrated,
                                                          ident.1 = i,
                                                          grouping.var = "groups")
  
  x <- as.data.frame(as.matrix(LIST.CLUSTERS.and.MARKERS[[i+1]]))  
  x$gene <- row.names(x)
  print(x$gene[1:20])  # Print top 20 genes to console
  
  write.table(x,
              file=paste("Table3.Arc.integrated.CONSERVED.MARKERS.cluster", i,
                         "LIST.txt", sep="."),
              sep="\t", quote=F, row.names=T, col.names=T)
  
  # DotPlot
  DotPlot(Arc.integrated, features = rev(x$gene[1:20]),
          cols = c("blue", "red"), dot.scale = 5) + RotatedAxis()
  ggsave(paste("fig4.Arc.integrated.CONSERVED.MARKERS.Cluster",i,"DotPlot.pdf",
               sep="."),width=10, height=8)
  
  # VlnPlots
  VlnPlot(Arc.integrated, features = rev(x$gene[10:1]), pt.size = 0, ncol = 2)
  ggsave(paste("fig5.Arc.integrated.CONSERVED.MARKERS.Cluster",i,"VlnPlot1.pdf",
               sep="."),width=10, height=10)
  
  VlnPlot(Arc.integrated, features = rev(x$gene[20:11]), pt.size = 0, ncol = 2)
  ggsave(paste("fig5.Arc.integrated.CONSERVED.MARKERS.Cluster",i,"VlnPlot2.pdf",
               sep="."),width=10, height=10)
  
  # FeaturePlots
  FeaturePlot(Arc.integrated, features =rev(x$gene[12:1]),
              min.cutoff = "q10", order = TRUE, ncol = 4)
  ggsave(paste("fig6.Arc.integrated.CONSERVED.MARKERS.Cluster",i,"FeaturePlot1.pdf",
               sep="."), width=15, height=9)
  
  FeaturePlot(Arc.integrated, features =rev(x$gene[20:13]),
              min.cutoff = "q10", order = TRUE, ncol = 4)
  ggsave(paste("fig6.Arc.integrated.CONSERVED.MARKERS.Cluster",i,"FeaturePlot2.pdf",
               sep="."), width=15, height=6)
}

################################################################################
# End of Script ----------------------------------------------------------------
################################################################################
