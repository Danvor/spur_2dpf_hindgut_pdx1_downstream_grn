#!/bin/bash
#
#SBATCH --job-name=TOBIAS
#SBATCH --cpus-per-task=16
#SBATCH --output=log_tobias.txt
echo "Usage: running_TOBIAS_footprinting.sh <prefix of ATACorrect and ScoreBigwig output files> <path to genome> <path to CRM bed file> <prefix of BINDetect output files"
date -R
module load TOBIAS/latest
module load bamtools/2.5.1
for fileA in nucl_free_nodup_paired*A_ATAC.bam; do
  fileB="${fileA/A_/B_}"
  srun bamtools merge -in "$fileA" -in "$fileB" -out $1_merged_ATAC.bam
done
srun TOBIAS ATACorrect --bam $1_merged_ATAC.bam --genome $2 --peaks $3 --regions-out $3 --outdir $1_corrections/ --prefix $1 --cores 16
srun TOBIAS ScoreBigwig --signal $1_corrections/$1_corrected.bw --regions $3 --output $1_footprints.bw --cores 16
srun TOBIAS BINDetect --signals $1_footprints.bw --motifs ~/Scripts/jaspar2020_vertebrate.jaspar --genome $2 --peaks $3 --outdir $4_motifs_1e-5/ --prefix $4_TF_1e-5 --cores 16 --motif-pvalue 1e-5 --naming name --cond_names $4 --skip-excel
date -R
