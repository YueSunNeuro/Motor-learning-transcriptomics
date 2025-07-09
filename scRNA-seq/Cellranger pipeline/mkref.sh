#!/bin/bash
#SBATCH --job-name=mkfastq_sample         # SLURM job name
#SBATCH --output=mkfastq_sample.out       # Standard output log
#SBATCH --error=mkfastq_sample.err        # Standard error log
#SBATCH --time=6:00:00                 # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=16             # Number of CPU cores
#SBATCH --mem=64G                      # Total memory

# -------------------------------
# Description:
# Converts raw Illumina BCL files to FASTQ using Cell Ranger's mkfastq.
# Requires sample sheet (.csv) and run folder containing BCL files.
# -------------------------------

# Load required modules (depends on HPC environment)
module load biology
module load bcl2fastq

# Add Cell Ranger to PATH (update this path as needed)
export PATH=$PATH:/path/to/cellranger-3.0.2

# Run cellranger mkfastq
cellranger mkfastq \
  --id=sample-out \                      # Output folder name
  --run=/BH2HY7BBXY/ \                   # Path to BCL run folder (update as needed)
  --csv=/sample.csv \                    # Sample sheet file (update as needed)
  --qc \                                 # Generate quality control reports
  --localcores=16 \                      # Number of CPU cores to use
  --localmem=64 \                        # Amount of RAM to allocate (in GB)
  --ignore-dual-index                    # Skip dual index check (use if only single index is used)
