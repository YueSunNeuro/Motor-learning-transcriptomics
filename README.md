# Analysis of Neural Activity and Gene Expression During Motor Learning

This repository contains code used in the publication:  
**"Motor learning drives region-specific transcriptomic remodeling in the motor cortex and dorsal striatum"**  
By Yue Sun et al. (2025)

---

##  1. Single-Cell RNA-seq Analysis (scRNA-seq)

###  `scRNAseq/`

This directory contains code for processing and analyzing single-cell RNA sequencing data from the mouse motor cortex and dorsal striatum.

####  `1_cellranger_pipeline/`
Scripts for running [CellRanger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) to process raw sequencing data.

- `create_custom_reference.sh`: Shell script to create a custom reference genome using CellRanger `mkref`.
- `Ai9-WPRE.fa`: FASTA file containing exogenous sequence for reference building.
- `Ai9-WPRE.gtf`: GTF annotation for exogenous sequence.
- `count.sh`: Example SLURM-compatible script for running CellRanger `count`.
- `sample.csv`: Example sample sheet for CellRanger input.

####  `2_seurat_analysis/`
R scripts for downstream analysis using the [Seurat](https://satijalab.org/seurat/) package.

- `seurat_analysis.R`: Main pipeline for QC, normalization, clustering, dimensionality reduction, and differential gene expression analysis.

---

##  2. Calcium Imaging Analysis

###  `calcium_imaging/`

This directory contains MATLAB scripts for analyzing **2-photon in vivo calcium imaging data** from **HTR3a-Cre mice** performing a **head-fixed forelimb reaching task**.

---

###  File Overview

#### `Htr3a_CaImaging_analysis_part1.m`
- **Purpose:** Prepares the dataset for analysis by assembling a `Data` structure.
- **Inputs:**
  - GCaMP fluorescence traces from ImageJ/FIJI
  - DeepLabCut `.csv` files for behavioral tracking
  - `.xlsx` files from BORIS for event annotations
  - `overview.mat` file containing file paths for each session
- **Output:** `Data` structure used in subsequent steps.

#### `Htr3a_CaImaging_analysis_part2.m`
- **Purpose:** Visualizes calcium activity aligned to behavioral events.
- **Input:** `Data` structure from part 1.
- **Output:** Plots and figures of trial-aligned activity.

#### `Htr3a_Ca_analysis_f.m`
- **Purpose:** Main function to extract, align, and process calcium traces.
- **Inputs:**
  - Paths to event, calcium, and DLC files
  - Pre- and post-event time windows
  - Behavioral event type
  - Optional modifiers
- **Outputs:**
  - Raw and z-scored calcium traces aligned to events
  - Trial-averaged data
  - Aligned DeepLabCut behavioral data

#### `Htr3a_Ca_analysis_SigMod.m`
- **Purpose:** Identifies significantly modulated neurons and summarizes activity.
- **Input:** `Data` structure from part 1.
- **Output:** `results_corr` structure with:
  - Number of significantly modulated neurons
  - Average baseline and movement ΔF/F values per mouse

---

###  Dependencies

- **MATLAB** R2021b or later
- **Image Processing Toolbox** (recommended)
- Required input data exported from:
  - **ImageJ/FIJI** (ROI Manager traces)
  - **DeepLabCut** (pose tracking `.csv` files)
  - **BORIS** (behavioral event `.xlsx` annotations)

---
