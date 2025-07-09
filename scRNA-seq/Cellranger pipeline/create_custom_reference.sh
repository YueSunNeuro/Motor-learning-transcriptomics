#!/bin/bash
#SBATCH --job-name=testgenome
#SBATCH --output=testgenome.out
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Download genome FASTA (GRCm38)
wget ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Download GTF annotation
wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gunzip Mus_musculus.GRCm38.93.gtf.gz

# Step 1: Append Ai9-WPRE.fa to genome fasta
cp Mus_musculus.GRCm38.dna.primary_assembly.fa Mus_musculus.GRCm38.dna.primary_assembly_Ai9-WPRE.fa
cat Ai9-WPRE.fa >> Mus_musculus.GRCm38.dna.primary_assembly_Ai9-WPRE.fa
tail -n 5 Mus_musculus.GRCm38.dna.primary_assembly_Ai9-WPRE.fa

# Step 2: Append tdtom.gtf to annotation GTF
cp Mus_musculus.GRCm38.93.gtf Mus_musculus.GRCm38.93.gtf
cat Ai9-WPRE.gtf >> Mus_musculus.GRCm38.93.Ai9-WPRE.gtf
tail -n 5 Mus_musculus.GRCm38.93.Ai9-WPRE.gtf

# Step 3: Run CellRanger mkref
module load biology
export PATH=$PATH:/path/to/cellranger-3.0.2
cellranger mkref --genome=cellranger_Ai9-WPRE \
                 --fasta=/home/groups/dingjun/Mus_musculus.GRCm38.dna.primary_assembly_Ai9-WPRE.fa \
                 --genes=/home/groups/dingjun/Mus_musculus.GRCm38.93.Ai9-WPRE.gtf \
                 --memgb 64