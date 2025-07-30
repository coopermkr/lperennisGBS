#!/bin/bash
#SBATCH --mail-type=END
#SBATCH --mail-user=c.kimballrhines001@umb.edu
#SBATCH -o stacks_%j.out
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -t 24:00:00

module load stacks/2.65

# Run stacks specifying the paths to our ID files
#ref_map.pl -T 24 --popmap ./map -o ./albusDef/ --samples ./

# Run stacks denovo
denovo_map.pl -T 24 --popmap ./map -o ./denovo/ --samples ./ --paired
