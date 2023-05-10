#!/bin/bash
#
#SBATCH --job-name=salmon
#SBATCH --cpus-per-task=8
#SBATCH --output=log_salmon.txt

module load salmon/0.11.3
#In folder with 12 fastq files, 6 for 3 replicates of Sp-Pdx1 MO (R1 and R2) and 6 for 3 replicates of untreated embryos (R1 and R2)
for r1_file in paired_Sp_48hpf_*; do
  r2_file="${r1_file/_R1_/_R2_}"
  out="${r1_file/_R1_/_both_}"
srun salmon quant -i ~/Data/Sp/Sp_Transcriptome_Index_oneLox/ -l A --seqBias --gcBias -1 "$r1_file" -2 "$r2_file" -o quants_oneLox/gc_bias_corrected/quant_"$out"
done 
