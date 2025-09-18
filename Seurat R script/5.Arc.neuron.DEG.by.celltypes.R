################################################################################
# Title: Identify DEGs in TRAP Neurons – Arc-TRAP Motor Cortex
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
#   This script performs differential expression analysis on Arc-TRAP single-cell
#   RNA-seq data (TRAP+ vs TRAP-, Early vs Ctrl, Late vs Ctrl) for each neuronal
#   cell type. Generates DEG tables and volcano plots.
################################################################################

# Load libraries ---------------------------------------------------------------
library(Seurat)
library(dplyr)
library(VennDiagram)
library(reticulate)
library(cowplot)
library(grid)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(plyr)

################################################################################
# 1. Load Seurat object 
################################################################################
Arc.Neuron <- readRDS("Arc.Neuron.TRAP.labeled.rds")
Arc.Neuron
# 51149 features across 25002 samples within 3 assays 
# Active assay: RNA (24311 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

################################################################################
# 2. Data normalization
################################################################################
Arc.Neuron <-NormalizeData(Arc.Neuron, normalization.method = "LogNormalize", scale.factor = 10000)
# FindMarkers calculates wrong fold-changes when not using Log-Normalization 
# https://github.com/satijalab/seurat/issues/3879

################################################################################
# 3. Filter genes
################################################################################
# Remove reporter and batch-sensitive genes
All.neuron.genes <- rownames(Arc.Neuron)
All.neuron.genes <- All.neuron.genes[!All.neuron.genes %in% c("wpre", "tdtomato","Neo" )]

# Remove mitochondrial genes
All.neuron.genes <- All.neuron.genes[!grepl("^mt-", All.neuron.genes)]

# Remove batch-sensitive genes identified between clean Ctrl batches
ctrl.DEG <- read.table("Table.Arc.clean.ctrl.batch.sensitive.genes.txt", header = TRUE) 
genelist <- All.neuron.genes[!All.neuron.genes %in% ctrl.DEG$gene]
write.table(genelist, file = "Table.genelist.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

################################################################################
# 4. Define TRAP groups
################################################################################
DefaultAssay(object = Arc.Neuron) <- "RNA"
Idents(Arc.Neuron) <- "CellType"
Arc.Neuron$TRAP.groups <- paste(Idents(Arc.Neuron),Arc.Neuron$groups,Arc.Neuron$TRAP, sep = "_")

Idents(Arc.Neuron) <- "TRAP.groups"

# List of neuronal cell types
CellType <- c("dSPN","iSPN","eSPN","Sst","PV","Th","Chat")

################################################################################
# 5. TRAP+ vs TRAP- in Ctrl (non-learning related activity induced DEGs)
################################################################################
ggplot_list <- list()
ctrl_P_vs_N_results <- list()

# Loop through each cell type
for (j in CellType) { 
  ctrl_P_vs_N <- FindMarkers(Arc.Neuron, 
                             ident.1 = paste(j, "_ctrl_P", sep = ""), 
                             ident.2 = paste(j, "_ctrl_N", sep = ""), 
                             features = rownames(genelist), 
                             min.pct = 0.1, 
                             logfc.threshold = 0)
  ctrl_P_vs_N <- as.data.frame(as.matrix(ctrl_P_vs_N))
  ctrl_P_vs_N$'p_val_adj' <- p.adjust(ctrl_P_vs_N$'p_val', method = "BH") # fdr
  ctrl_P_vs_N$'abs_log2FC' <- abs(ctrl_P_vs_N$avg_log2FC)
  ctrl_P_vs_N$'cluster' <- j
  ctrl_P_vs_N$'Comparison' <- "ctrl_P_vs_N"
  ctrl_P_vs_N$'gene' <- rownames(ctrl_P_vs_N)

  Gene <- subset(ctrl_P_vs_N,abs_log2FC >log2(1.2)& p_val_adj <0.05, 
                 select = c(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, Comparison))
  ctrl_P_vs_N_results[[j]] <- Gene
  write.table(Gene, file=paste("Table.9.0.Arc.Neuron.ctrl_P_vs_N",j,"DE.LIST.Padj.txt", sep="."),
              sep="\t", quote=F, row.names=T, col.names=T)
  
  # Count upregulated and downregulated genes
  upregulated_ctrl <- nrow(subset(ctrl_P_vs_N, avg_log2FC > log2(1.2) & p_val_adj < 0.05))
  downregulated_ctrl <- nrow(subset(ctrl_P_vs_N, avg_log2FC < -log2(1.2) & p_val_adj < 0.05))
  
  # Store the counts
  upregulated_counts[[paste(j, "ctrl_P_vs_N", sep = "_")]] <- upregulated_ctrl
  downregulated_counts[[paste(j, "ctrl_P_vs_N", sep = "_")]] <- downregulated_ctrl
}

ctrl_combined <- ldply(ctrl_P_vs_N_results, data.frame)
write.table(ctrl_combined, file=paste("Table.7.Arc.Neuron.ctrl_P_vs_N.DE.LIST.Padj.txt", sep="."),
            sep="\t", quote=F, row.names=T, col.names=T)

# Convert lists to data frames for easier viewing and saving
upregulated_df <- as.data.frame(do.call(rbind, upregulated_counts))
downregulated_df <- as.data.frame(do.call(rbind, downregulated_counts))

# Write the counts to files
write.table(upregulated_df, file="Table.8.Arc.Neuron.P_vs_N.Upregulated_Gene_Counts.Padj.txt", 
                            sep="\t", quote=F, col.names=NA)
write.table(downregulated_df, file="Table.8.Arc.Neuron.P_vs_N.Downregulated_Gene_Counts.Padj.txt", 
                            sep="\t", quote=F, col.names=NA) 

################################################################################
# 6. Early TRAP vs Ctrl TRAP, Late TRAP vs Ctrl TRAP
################################################################################
ggplot_list <- list()
early_vs_ctrl_results <- list()
late_vs_ctrl_results <- list()
upregulated_counts <- list()
downregulated_counts <- list()

# Loop through each cell type
for (j in CellType) { 
  ## Early vs Ctrl
  early_vs_ctrl <- FindMarkers(Arc.Neuron, 
                               ident.1 = paste(j, "_early_P", sep = ""), 
                               ident.2 = paste(j, "_ctrl_P", sep = ""), 
                               features = rownames(genelist), 
                               logfc.threshold = 0)
  early_vs_ctrl <- as.data.frame(as.matrix(early_vs_ctrl))
  early_vs_ctrl$'p_val_adj' <- p.adjust(early_vs_ctrl$'p_val', method = "BH") # fdr
  early_vs_ctrl$'abs_log2FC' <- abs(early_vs_ctrl$avg_log2FC)
  early_vs_ctrl$'cluster' <- j
  early_vs_ctrl$'Comparison' <- "early_vs_ctrl"
  early_vs_ctrl$'gene' <- rownames(early_vs_ctrl)

  Gene <- subset(early_vs_ctrl,abs_log2FC >log2(1.2)& p_val_adj <0.05, 
                 select = c(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, Comparison))
  early_vs_ctrl_results[[j]] <- Gene
  write.table(Gene, file=paste("Table.9.1.Arc.Neuron.early.vs.ctrl",j,"DE.LIST.txt", sep="."),
              sep="\t", quote=F, row.names=T, col.names=T)

  # Count upregulated and downregulated genes
  upregulated_early <- nrow(subset(early_vs_ctrl, avg_log2FC > log2(1.2) & p_val < 0.05))
  downregulated_early <- nrow(subset(early_vs_ctrl, avg_log2FC < -log2(1.2) & p_val < 0.05))
  
  # Store the counts
  upregulated_counts[[paste(j, "early_vs_ctrl", sep = "_")]] <- upregulated_early
  downregulated_counts[[paste(j, "early_vs_ctrl", sep = "_")]] <- downregulated_early

  # Find markers for late.vs.ctrl comparison
  late_vs_ctrl <- FindMarkers(Arc.Neuron, 
                              ident.1 = paste(j, "_late_P", sep = ""), 
                              ident.2 = paste(j, "_ctrl_P", sep = ""), 
                              features = rownames(genelist), 
                              logfc.threshold = 0)
  late_vs_ctrl <- as.data.frame(as.matrix(late_vs_ctrl))
  late_vs_ctrl$'p_val_adj' <- p.adjust(late_vs_ctrl$'p_val',  method = "BH") # fdr
  late_vs_ctrl$'abs_log2FC' <- abs(late_vs_ctrl$avg_log2FC)
  late_vs_ctrl$'cluster' <- j
  late_vs_ctrl$'Comparison' <- "late_vs_ctrl"
  late_vs_ctrl$'gene' <- rownames(late_vs_ctrl)

  Gene <- subset(late_vs_ctrl,abs_log2FC >log2(1.2)& p_val_adj <0.05, 
                 select = c(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj,cluster,Comparison))
  late_vs_ctrl_results[[j]] <- Gene
  write.table(Gene, file=paste("Table.9.2.Arc.Neuron.late.vs.ctrl",j,"DE.LIST.txt", sep="."),
                    sep="\t", quote=F, row.names=T, col.names=T)

  # Count upregulated and downregulated genes
  upregulated_late <- nrow(subset(late_vs_ctrl, avg_log2FC > log2(1.2) & p_val < 0.05))
  downregulated_late <- nrow(subset(late_vs_ctrl, avg_log2FC < -log2(1.2) & p_val < 0.05))
  
  # Store the counts
  upregulated_counts[[paste(j, "late_vs_ctrl", sep = "_")]] <- upregulated_late
  downregulated_counts[[paste(j, "late_vs_ctrl", sep = "_")]] <- downregulated_late

  # Volcano plots
  p1 <- EnhancedVolcano(early_vs_ctrl, lab = rownames(early_vs_ctrl), labSize = 1,
                        x = 'avg_log2FC', y = 'p_val_adj',
                        xlim = c(-1.5, 1.5), #ylim = c(0, 70),
                        title = 'TRAP Early vs ctrl',
                        subtitle = j,
                        col = c("grey30", "grey30", "grey30", "red2"),
                        pCutoff = 50e-3,
                        FCcutoff = log2(1.2),
                        pointSize = 1.0,
                        gridlines.major =F, gridlines.minor =F)
  
  p2 <- EnhancedVolcano(late_vs_ctrl, lab = rownames(late_vs_ctrl), labSize = 1,
                        x = 'avg_log2FC', y = 'p_val_adj',
                        xlim = c(-1.5, 1.5), #ylim = c(0, 15),
                        title = 'TRAP Late vs Ctrl',
                        subtitle = j,
                        col = c("grey30", "grey30", "grey30", "red2"),
                        pCutoff = 50e-3,
                        FCcutoff = log2(1.2),
                        pointSize = 1.0,
                        gridlines.major =F, gridlines.minor =F)
  ggplot_list[[j]] <- plot_grid(p1, p2, ncol = 2)
}


# Save all plots to a single PDF file
pdf("Ext-1_groups.Volcanoplots.pdf", width = 12, height = 32.5) # Adjust width and height as needed
grid.arrange(grobs = ggplot_list[1:5], ncol = 1) # Adjust ncol to control number of columns per row
dev.off()

pdf("Ext-2_groups.Volcanoplots.pdf", width = 12, height = 32.5) # Adjust width and height as needed
grid.arrange(grobs = ggplot_list[6:10], ncol = 1) # Adjust ncol to control number of columns per row
dev.off()

pdf("Int-1_groups.Volcanoplots.pdf", width = 12, height = 32.5) # Adjust width and height as needed
grid.arrange(grobs = ggplot_list[11:15], ncol = 1) # Adjust ncol to control number of columns per row
dev.off()

pdf("Int-2_groups.Volcanoplots.pdf", width = 12, height = 26) # Adjust width and height as needed
grid.arrange(grobs = ggplot_list[16:19], ncol = 1) # Adjust ncol to control number of columns per row
dev.off()

# Save all plots to a single PDF file
pdf("SPN_Volcanoplots.P.adj.ylim.pdf", width = 12, height = 19.5) # Adjust width and height as needed
grid.arrange(grobs = ggplot_list[1:3], ncol = 1) # Adjust ncol to control number of columns per row
dev.off()

pdf("Int_Volcanoplots.P.adj.ylim.pdf", width = 12, height = 25) # Adjust width and height as needed
grid.arrange(grobs = ggplot_list[4:7], ncol = 1) # Adjust ncol to control number of columns per row
dev.off()

early_combined <- ldply(early_vs_ctrl_results, data.frame)
late_combined <- ldply(late_vs_ctrl_results, data.frame)

write.table(early_combined, file=paste("Table.7.Arc.Neuron.early_vs_ctrl.DE.LIST.txt", sep="."),sep="\t", quote=F, row.names=T, col.names=T)
write.table(late_combined, file=paste("Table.7.Arc.Neuron.late_vs_ctrl.DE.LIST.txt", sep="."),sep="\t", quote=F, row.names=T, col.names=T)

# Convert lists to data frames for easier viewing and saving
upregulated_df <- as.data.frame(do.call(rbind, upregulated_counts))
downregulated_df <- as.data.frame(do.call(rbind, downregulated_counts))

# Write the counts to files
write.table(upregulated_df, file="Table.8.Arc.Neuron.groups.Upregulated_Gene_Counts.txt", sep="\t", quote=F, col.names=NA)
write.table(downregulated_df, file="Table.8.Arc.Neuron.groups.Downregulated_Gene_Counts.txt", sep="\t", quote=F, col.names=NA)


