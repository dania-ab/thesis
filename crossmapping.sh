# https://github.com/Gabaldonlab/crossmapper
# Installation
# conda install -c gabaldonlab -c bioconda crossmapper

# Running
crossmapper RNA -t 4 -rlay PE -rlen 150 -o ~/Desktop -rc TRUE -gn CA H SC -g fasta CA.fna human.fna SC.fna  -a gff albicans.gff human.gff sc.gff
