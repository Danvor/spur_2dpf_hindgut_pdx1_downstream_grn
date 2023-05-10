#!/bin/bash
#
#SBATCH --job-name=bowtie2_index
#SBATCH --cpus-per-task=16
#SBATCH --output=log_bowtie2.txt
module load bowtie/2.3.4.1
#$1 is genome fasta
#$2 is prefix for genome name
mkdir $2_bowtie_ind
srun bowtie2-build --threads 32 $1 $2_bowtie_ind/$2
