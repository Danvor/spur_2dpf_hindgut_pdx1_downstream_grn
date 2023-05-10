#!/bin/bash
#
#SBATCH --job-name=cellranger_count
#SBATCH --cpus-per-task=32
#SBATCH --output=log_count_cr.txt

#$1 prefix name for sample
#$2 number of cells forced based on the web summary Barcode Rank Plot
#$3 additional identifier
echo "CellRanger 3.0.2 running for S. purpuratus."

srun cellranger count --id=$1_$2_$3 --transcriptome=~/Data/Sp/SingleCell/Index/Sp_cellranger_index --fastqs=~/Data/Sp/SingleCell/FASTQs/$1 --sample=$1 --force-cells=$2 --nosecondary
cd ./$1_$2_$3/outs/filtered_feature_bc_matrix
gzip -d *.gz
mv features.tsv genes.tsv
cd ./$1_$2_$3/outs/
mkdir $1_$2
cp web_summary.html $1_$2
cp filtered_feature_bc_matrix $1_$2
