#!/bin/bash
#SBATCH -o trimmo_%j.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -t 3:00:00

module load trimmomatic/0.39

echo ${1}
dir=/project/pi_brook_moyers_umb_edu/lupine/rawData

trimmomatic PE -threads 2 -basein $dir/${1}_R1_001.fastq.gz trim.${1}.p.1.fq.gz trim.${1}.u.1.fq.gz trim.${1}.p.2.fq.gz trim.${1}.u.2.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 MINLEN:36
