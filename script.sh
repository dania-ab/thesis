#!/bin/bash

SECONDS=0

# Preparation and setup of required files 

FILES=$(pwd)
WKDIR=$(echo $FILES | sed 's:/required_files::g')

GENOME=$WKDIR/required_files/refgenomes.fasta
FEATURES_H=$WKDIR/required_files/human.gff
FEATURES_C=$WKDIR/required_files/sc.gff
FEATURES_S=$WKDIR/required_files/albicans.gff
mkdir $WKDIR/QC
PICARD=$WKDIR/required_files/picard.jar
rRNA_H=$WKDIR/required_files/rRNA_H.gff
rRNA_C=$WKDIR/required_files/rRNA_CA.gff
rRNA_S=$WKDIR/required_files/rRNA_SC.gff

# Prompts

echo 'Are the data in bam or fastq format?'
read -p 'Specify file format (bam or fastq): ' FORMAT
read -p 'Do you want to do a quality control of the raw data (yes or no): ' QCRAW
read -p 'How many threads (cores) should be used for the analysis (use 1 if you are not sure): ' THREAD

#### (1) QC of raw data ####

if [ $QCRAW == 'yes' ]
then 
mkdir $WKDIR/QC_raw
echo 'Quality control of raw data:'
if [ $FORMAT == 'bam' ]
then
for i in $WKDIR/*.bam
do
fastqc -o $WKDIR/QC_raw $i
done
else	
for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(q\.gz$)')
do
i=$WKDIR/$SNAME
fastqc -o $WKDIR/QC_raw $i
done
fi
multiqc -o $WKDIR/QC_raw $WKDIR/QC_raw
else
echo 'No QC of raw data done.'
fi

# Convert .bam to .fastq format

if [ $FORMAT == 'bam' ]
then
echo 'File format is bam.'
for i in $WKDIR/*.bam
do
bamToFastq -i $i -fq $i.fq
done
elif [ $FORMAT == 'fastq' ]
then
echo 'File format is fastq.'
else
echo 'Invalid file format! Options are "bam" or "fastq".'
exit
fi

#### (2) Mapping of all files with STAR ####

for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(L*_1\.fq\.gz$)')
do

i1=$WKDIR/$SNAME
i2=$(echo $i1 | sed 's/_1.fq.gz/_2.fq.gz/')
i=$($SNAME | sed 's/_L.*//') # test !!

star --genomeDir ~/Desktop $i1 $i2 --readFilesCommand gunzip -c --outFileNamePrefix $i --outSAMtype BAM SortedByCoordinate

#### (3) Further processing of BAM files ####

mv $i1.count.txt $WKDIR/required_files

# Removal of rRNA

intersectBed -v -abam $i.STARAligned.sortedByCoord.out.bam -b $rRNA_H > $i.rRNA_H.bam 
intersectBed -v -abam $i.rRNA_H.bam -b $rRNA_C > $i.rRNA_H_C.bam 
intersectBed -v -abam $i.rRNA_H_C.bam -b $rRNA_S > $i.rRNA.final.bam 

# Labelling of duplicated reads and removal of optical duplicates
java -jar $PICARD MarkDuplicates -REMOVE_SEQUENCING_DUPLICATES true -I $i.rRNA.final.bam -O $i.final.bam -M $WKDIR/QC/$i.markdup.metrics.txt

# Index final bam file
samtools index $i.final.bam 

# Quality control and statistics about mapped samples
samtools flagstat $i.final.bam >> $WKDIR/QC/$i.final.flagstat_analysis.txt   # flagstat analysis

fastqc $i.final.bam -o $WKDIR/QC 
done


#### (4) Count reads with HTSeq ####

mkdir $WKDIR/count
mkdir $WKDIR/diff_expr_analysis

for SNAME in $WKDIR/*.final.bam
do
i=$($SNAME | sed 's/.final.bam//')
htseq-count -f bam -r pos -s no -t gene -i ID $WKDIR/$SNAME $FEATURES_H > $i.H.count.txt
mv $i.H.count.txt $WKDIR/count
echo "Human transcripts done"
done

for SNAME in $WKDIR/*.final.bam
do
i=$($SNAME | sed 's/.final.bam//')
htseq-count -f bam -r pos -s no -t gene -i ID $WKDIR/$SNAME $FEATURES_C > $i.CA.count.txt
mv $i.CA.count.txt $WKDIR/count
echo "Candida transcripts done"
done

for SNAME in $WKDIR/*.final.bam
do
i=$($SNAME | sed 's/.final.bam//')
htseq-count -f bam -s no -t gene -i ID $WKDIR/$SNAME $FEATURES_S > $i.SC.count.txt
mv $i.SC.count.txt $WKDIR/count
echo "Saccharomyces transcripts done"
done

for i in $WKDIR/count/*.count.txt
do
sed '$d' "$i" | sed '$d' | sed '$d' | sed '$d' | sed '$d' > "$i.crop.txt"
done

duration=$SECONDS
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."
