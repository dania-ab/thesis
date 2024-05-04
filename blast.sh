#!/bin/bash

SECONDS=0
duration=$SECONDS

blastn -query sc.fasta -subject ca.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out SC_H.txt
blastn -query sc.fasta -subject human.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out SC_CA.txt
blastn -query ca.fasta -subject human.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out CA_H.txt

echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."
