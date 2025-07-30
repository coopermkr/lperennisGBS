#!/bin/bash
for file in ../0.trim/*.fastq.gz;
do
   echo ${file}
   name=$(cut -d '/' -f 3 <<< ${file} | cut -d '_' -f 1)
   echo ${name}
   sbatch -o align.$name.log ../scripts/alignment.sh ${name}
   sleep 3
done

