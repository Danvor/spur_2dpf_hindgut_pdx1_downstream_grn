#!/bin/bash
#
#SBATCH --job-name=salmon
#SBATCH --cpus-per-task=16
#SBATCH --output=log_salmon.txt

module load salmon/0.11.3
#$1 is the transcriptome fasta file
#$2 is the index folder name
echo "Usage: salmon index -t <path to transcriptome fasta> -i <path to index output folder> --type quasi -k 25"
srun salmon index -t $1 -i $2 --type quasi -k 25

