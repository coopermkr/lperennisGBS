#!/bin/bash
while read in /project/pi_brook_moyers_umb_edu/lupine/rawData/*_R1_001.fastq.gz; 
do
   echo ${file}
   name=$(cut -d '/' -f 6 <<< ${file} | cut -d '_' -f 1)
   echo ${name}
   sbatch -o trim.$name.log ../scripts/trimmo.sh ${name}
   sleep 1
done

