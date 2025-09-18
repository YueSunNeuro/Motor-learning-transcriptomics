################################################################################
# Title: Identify DEGs in Glia – Motor Cortex and Striatum
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
#   This script performs differential expression analysis (Early vs Ctrl, Late vs Ctrl) 
#   for each glial cell type in the motor cortex, generating DEG tables and volcano plots. 
#   The same workflow can be applied to striatum (Arc-TRAP) data by substituting input files. 
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
Fos.glia <- readRDS("Fos.Glia.rds")
Fos.glia

################################################################################
# 2. Normalize data
################################################################################
Fos.glia <-NormalizeData(Fos.glia, normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
# FindMarkers calculates wrong fold-changes when not using Log-Normalization 
# https://github.com/satijalab/seurat/issues/3879

################################################################################
# 3. Filter genes
################################################################################
# Remove reporter and batch-sensitive genes
All.glia.genes <- rownames(Fos.glia)
All.glia.genes <- All.glia.genes[!All.glia.genes %in% c("wpre", "tdtomato","Neo" )]

# Remove mitochondrial genes
All.glia.genes <- All.glia.genes[!grepl("^mt-", All.glia.genes)]

# Remove batch-sensitive genes identified between clean Ctrl batches
ctrl.DEG <- read.table("Table.Fos.clean.ctrl.batch.sensitive.genes.txt", header = TRUE) 
genelist <- All.glia.genes[!All.glia.genes %in% ctrl.DEG$gene]

################################################################################
# 4. Define comparation groups
################################################################################
DefaultAssay(object = Fos.glia) <- "RNA"
Idents(Fos.glia) <- "CellType"

# List of glial cell types
CellType <- c("OL","OPC","AS","MG")

# Create combined group labels
Fos.glia$Cell.groups <- paste(Idents(Fos.glia),Fos.glia$groups, sep = "_")

Idents(Fos.glia) <- "Cell.groups"
table(Idents(Fos.glia))

################################################################################
# 5. Perform differential expression analysis
################################################################################
ggplot_list <- list()
early_vs_ctrl_results <- list()
late_vs_ctrl_results <- list()
upregulated_counts <- list()
downregulated_counts <- list()

# Loop through each cell type
for (j in CellType) { 
  ## Early vs Ctrl
  early_vs_ctrl <- FindMarkers(Fos.glia, 
                              ident.1 = paste(j, "early", sep = ""), 
                              ident.2 = paste(j, "ctrl", sep = ""), 
                              features = rownames(genelist), 
                              min.pct = 0.1, 
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
  write.table(Gene, file=paste("Table.9.1.Fos.glia.early.vs.ctrl",j,"DE.LIST.txt", sep="."),
              sep="\t", quote=F, row.names=T, col.names=T)

  # Count upregulated and downregulated genes
  upregulated_early <- nrow(subset(early_vs_ctrl, avg_log2FC > log2(1.2) & p_val < 0.05))
  downregulated_early <- nrow(subset(early_vs_ctrl, avg_log2FC < -log2(1.2) & p_val < 0.05))
  
  # Store the counts
  upregulated_counts[[paste(j, "early_vs_ctrl", sep = "_")]] <- upregulated_early
  downregulated_counts[[paste(j, "early_vs_ctrl", sep = "_")]] <- downregulated_early

  ## Late vs Ctrl 
  late_vs_ctrl <- FindMarkers(Fos.glia, 
                              ident.1 = paste(j, "late", sep = ""), 
                              ident.2 = paste(j, "ctrl", sep = ""), 
                              features = rownames(genelist), 
                              min.pct = 0.1, 
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
  write.table(Gene, file=paste("Table.9.2.Fos.glia.late.vs.ctrl",j,"DE.LIST.txt", sep="."),
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
                        xlim = c(-2, 2), #ylim = c(0, 70),
                        title = 'Early vs ctrl',
                        subtitle = j,
                        col = c("grey30", "grey30", "grey30", "red2"),
                        pCutoff = 50e-3,
                        FCcutoff = log2(1.2),
                        pointSize = 1.0,
                        gridlines.major =F, gridlines.minor =F)
  
  p2 <- EnhancedVolcano(late_vs_ctrl, lab = rownames(late_vs_ctrl), labSize = 1,
                        x = 'avg_log2FC', y = 'p_val_adj',
                        xlim = c(-2, 2), #ylim = c(0, 15),
                        title = 'Late vs Ctrl',
                        subtitle = j,
                        col = c("grey30", "grey30", "grey30", "red2"),
                        pCutoff = 50e-3,
                        FCcutoff = log2(1.2),
                        pointSize = 1.0,
                        gridlines.major =F, gridlines.minor =F)
  ggplot_list[[j]] <- plot_grid(p1, p2, ncol = 2)
}

# Save all plots to a single PDF file
pdf("Fos.all.glia.groups.Volcanoplots.Padj.pdf", width = 12, height = 26) 
grid.arrange(grobs = ggplot_list, ncol = 1) 
dev.off()

################################################################################
# 6. Save combined DEG tables and counts
################################################################################
early_combined <- ldply(early_vs_ctrl_results, data.frame)
late_combined <- ldply(late_vs_ctrl_results, data.frame)

write.table(early_combined, file=paste("Table.7.Fos.all.glia.early_vs_ctrl.DE.LIST.Padj.txt", sep="."),
                            sep="\t", quote=F, row.names=T, col.names=T)
write.table(late_combined, file=paste("Table.7.Fos.all.glia.late_vs_ctrl.DE.LIST.Padj.txt", sep="."),
                            sep="\t", quote=F, row.names=T, col.names=T)

# Convert lists to data frames for easier viewing and saving
upregulated_df <- as.data.frame(do.call(rbind, upregulated_counts))
downregulated_df <- as.data.frame(do.call(rbind, downregulated_counts))

# Write the counts to files
write.table(upregulated_df, file="Table.8.Fos.all.glia.groups.Upregulated_Gene_Counts.Padj.txt", 
                            sep="\t", quote=F, col.names=NA)
write.table(downregulated_df, file="Table.8.Fos.all.glia.groups.Downregulated_Gene_Counts.Padj.txt", 
                            sep="\t", quote=F, col.names=NA)

################################################################################
# 7. Repeat for striatum (replace Fos with Arc-TRAP objects)
################################################################################
