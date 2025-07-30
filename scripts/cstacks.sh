#!/bin/bash
#SBATCH --mail-type=END
#SBATCH --mail-user=c.kimballrhines001@umb.edu
#SBATCH -o stacks_%j.out
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -t 24:00:00

module load stacks/2.65

#Build SNP catalog
#cstacks -n 6 -P ./ -M novoMap -p 32

#Match samples against catalog
#sstacks -P ./ -M novoMap -p 32

#Transpose data so it is stord by locus
tsv2bam -P ./ -M novoMap --pe-reads-dir ./ -t 32

#Build paired-end contigs from the metapopulation
gstacks -P ./ -M novoMap -t 32

#Calculate population statistics
populations -P ./ -M novoMap -t 32 -r 0.65 --vcf --genepop --structure --fstats --hwe


