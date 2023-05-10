#!/bin/bash
#
#SBATCH --job-name=trimm
#SBATCH --cpus-per-task=8
#SBATCH --output=log_trimm.txt
date -R
module load trimmomatic/0.38
#Trimming RNA-seq bad quality reads
#In folder with raw RNA-seq fastq files
for r1_file in *R1_001.fastq.gz; do
  r2_file="${r1_file/_R1_/_R2_}"
  srun java -jar /opt/trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 32 -phred33 "$r1_file" "$r2_file" paired_"$r1_file" unpaired_"$r1_file" paired_"$r2_file" unpaired_"$r2_file" ILLUMINACLIP:~/Data/Sp/RNAseq/TruSeq3-PE-2-mod.fa:2:30:10:6 CROP:90 HEADCROP:10 SLIDINGWINDOW:3:25 MINLEN:25
done
srun rm unpaired_*
date -R
