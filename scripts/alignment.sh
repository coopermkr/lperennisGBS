#!/bin/bash
#SBATCH -o align_ddrad%j.out
#SBATCH -c 5
#SBATCH --mem=8G
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -t 6:00:00

echo ${1}
ref=/project/pi_brook_moyers_umb_edu/lupine/reference/hap2/lupinehap2.fasta
trim=../0.trim

module load bwa/0.7.17
module load samtools/1.14

# Provide an index file, comment out line below once it has run once
#bwa index $ref

# Align the reads selected by the loop file
bwa mem -t 5 -P $ref $trim/${1}.1.fq.gz $trim/${1}.2.fq.gz > ${1}.sam

# Convert sam to bam
samtools view -@ 5 -F 0x4 -q 20 -Sb -o ${1}.bam ${1}.sam

# Sort bam
samtools sort -@ 5 -O bam -o sorted.${1}.bam ${1}.bam

# Calculate statistics
samtools stats -@ 5 ${1}.sam > stats.${1}.txt

echo done
