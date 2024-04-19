#!/bin/bash

# Get the coverage percentages

for SNAME in $(ls | egrep '\.bam$') 
do 
i=$SNAME 
bedtools genomecov -ibam $i > $i.coverage.txt 
done

# Manually compile 0-4X coverage over the whole genome (replace chromosomename by the name of each chromosome)

for SNAME in $(ls | egrep '\.coverage.txt$') 
do 
i=$SNAME 
grep 'chromosomename' $i | head -5 >> $i.crop.txt 
done 
