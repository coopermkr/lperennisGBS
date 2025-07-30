#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=20000
#SBATCH -p cpu
#SBATCH -t 24:00:00
#SBATCH -o fastqc-reverse%j.out

module load fastqc/0.11.9

fastqc --threads 20 /project/pi_brook_moyers_umb_edu/lupine/rawData/*R2*.fastq.gz -o .
