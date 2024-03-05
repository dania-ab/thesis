#!/bin/bash

SECONDS = 0

# change working directory
# cd /path/
# within your working folder, other folders: HISAT2 data quants script

# (1) Run fastqc
fastqc folder/read.fastq -o folder/

# (2) Run Hisat2

# (3) Run alignment

# (4) Run featureCounts - Quantification

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
