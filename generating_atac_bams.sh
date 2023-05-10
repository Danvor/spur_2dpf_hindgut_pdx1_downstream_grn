#!/bin/bash
#
#SBATCH --job-name=atac_bam
#SBATCH --cpus-per-task=16
#SBATCH --output=log_atacbam.txt

module load trimmomatic/0.38
module load bowtie/2.3.4.1
module load samtools/1.7.1
module load deeptools/3.1.3
module load bamtools/2.5.1

date -R
#Trim the bad quality reads
for r1_file in *R1.fq.gz; do
  r2_file="${r1_file/_R1/_R2}"
  srun java -jar /opt/trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 16 -phred33 "$r1_file" "$r2_file" paired_"$r1_file" unpaired_"$r1_file" paired_"$r2_file" unpaired_"$r2_file" ILLUMINACLIP:~/Scripts/Adapters.fa:2:30:10:6 CROP:40 SLIDINGWINDOW:3:25 MINLEN:25
done
srun rm unpaired_*

#Map the reads and fix mates
for file1 in paired_*_R1.fq.gz; do
  file2="${file1/_R1/_R2}"
  out_file="${file1/_R1.fq.gz/.bam}"
  srun bowtie2 -t --no-unal --no-mixed -X 2000 -p 16 -x ~/Data/Sp/bowtieSpur3.1ind -1 "$file1" -2 "$file2" | samtools view -b -q 30 - | samtools sort -n - | samtools fixmate -m - - | samtools sort -@ 16 -o temp_"$out_file"
done

#Remove duplicates
for file3 in temp_*; do
  out="${file3/temp_/nodup_}"
  srun samtools markdup -r "$file3" "$out"
done

#Check the success of the markdup command and remove temporary files
for file4 in nodup_*; do
  [[ $(find "$file4" -type f -size +50M) ]] && rm temp_* || echo "Error!!!"
done

#Filter bam by size
for file5 in nodup_*; do
  srun samtools view -h "$file5" | awk -v LEN=130 '{if ($9 <= LEN && $9 >= -(LEN) && $9 != 0 || $1 ~ /^@/) print $0}' | samtools sort -O bam -o nucl_free_"$file5" -
done

#Index nucleosome free bam files
for file6 in nucl_free_*bam; do
  srun samtools index "$file6"
done

#Get tracks for nucleosome free files
#Effective genome size for for Sp3.1 is 815936258
for file7 in nucl_free_*; do
out="${file7/.bam/.bw}"
srun bamCoverage -of bigwig -p max -bs 1 --effectiveGenomeSize 815936258 -b "$file7" -o "$out"
done

date -R

#merging the bam files
for file8 in nucl_free*A_ATAC.bam; do
  file9="${file8/A_/B_}"
  out="${file8/A_/merged_}"
srun bamtools merge -in "$file8" -in "$file9" -out "$out"
done
