#################################################################################
# Title: Fisher’s Exact Test on TRAP vs non-TRAP Counts by CellType
# Author: Yue Sun
# Date: 2025-09-13
#
# Description:
#   This script reads in summary tables (e.g. TRAP_nonTRAP_counts_byCellType.csv)
#   and performs Fisher’s exact tests for each row (cell type or cluster) to test
#   whether the proportion of TRAP vs non-TRAP cells differs between groups.
#
# Input: TRAP_nonTRAP_counts_byCellType.csv
# Output: tab-delimited files with raw and FDR-corrected p-values
# Note: Showing motor cortex (Fos-TRAP) here. The same workflow can be applied to 
# striatum (Arc-TRAP) data. 
################################################################################

# Read CSV
data <- read.csv("TRAP_nonTRAP_counts_byCellType.csv", header = TRUE, row.names = 1, 
                  check.names = FALSE)
head(data)

##########################################
# Ctrl vs Early
##########################################
p_values <- numeric(nrow(data))  # initialize vector for p-values

for (i in 1:nrow(data)) {
  # Extract counts for Ctrl vs Early
  counts <- c(
    data[i, "Ctrl-TRAP"], data[i, "Ctrl-nonTRAP"],
    data[i, "Early-TRAP"], data[i, "Early-nonTRAP"]
  )
  # Reshape counts to 2x2 matrix
  counts_matrix <- matrix(counts, nrow = 2, byrow = TRUE,
                          dimnames = list(c("Ctrl","Early"), c("TRAP","nonTRAP")))
  
  # Fisher's exact test
  result <- fisher.test(counts_matrix)
  p_values[i] <- result$p.value
}

# FDR correction
p_values_fdr <- p.adjust(p_values, method = "fdr")

# Add to data frame
data$Ctrl_vs_Early_pvalue <- p_values
data$Ctrl_vs_Early_pvalue_fdr <- p_values_fdr

# Save results
write.table(data, file="Fisher.Ctrl_vs_Early.by.groups.result.txt",
            sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

##########################################
# Ctrl vs Late
##########################################
p_values <- numeric(nrow(data))  # reset p-values

for (i in 1:nrow(data)) {
  # Extract counts for Ctrl vs Late
  counts <- c(
    data[i, "Ctrl-TRAP"], data[i, "Ctrl-nonTRAP"],
    data[i, "Late-TRAP"], data[i, "Late-nonTRAP"]
  )
  # Reshape counts to 2x2 matrix
  counts_matrix <- matrix(counts, nrow = 2, byrow = TRUE,
                          dimnames = list(c("Ctrl","Late"), c("TRAP","nonTRAP")))
  
  # Fisher's exact test
  result <- fisher.test(counts_matrix)
  p_values[i] <- result$p.value
}

# FDR correction
p_values_fdr <- p.adjust(p_values, method = "fdr")

# Add to data frame
data$Ctrl_vs_Late_pvalue <- p_values
data$Ctrl_vs_Late_pvalue_fdr <- p_values_fdr

# Save results
write.table(data, file="Fisher.Ctrl_vs_Late.by.groups.result.txt",
            sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
