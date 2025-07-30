#!/bin/bash
#SBATCH --mail-type=END
#SBATCH --mail-user=c.kimballrhines001@umb.edu
#SBATCH -o bam_stats%j.out
#SBATCH -c 2
#SBATCH --mem=5G
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -t 2:00:00

# Load samtools
module load samtools/1.14

# Run stats on all alignments
for file in ./sorted.trimmed.*.bam;
do
   echo ${file}
   samtools stats -@ 2 ${file}
   sleep 1
done
