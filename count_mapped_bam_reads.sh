#!/bin/sh
#counts the number of mapped reads in the BAM file for each replicate 
repnames=( 10-CC_rep2 9-CC_rep1 18-H1_rep1 19-H1_rep2 20-H1_rep3 21-48hr_rep1 7-48hr_rep2 8-48hr_rep3 4-3hr_rep1 5-3hr_rep2 6-3hr_rep3 1-16hr_rep1 2-16hr_rep2 3-16hr_rep3 11-ES_rep1 12-ES_rep2 13-Hk_rep1 14-Hk_rep2 22-M5_rep1 15-M5_rep2 16-M5_rep3 17-M5_rep4 )
for rep in "${repnames[@]}"
do
echo "$rep `samtools view -F 0x40 /srv/scratch/annashch/stemcells/het/peaks/peaks/$rep/*nonchrM.bam | cut -f1 | sort | uniq | wc -l &`" 
done
