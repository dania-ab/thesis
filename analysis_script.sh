#!/bin/bash

# Preparation and setup of required files

FILES=$(pwd)
WKDIR=$(echo $FILES | sed 's:/required_files::g')

read -p 'Do you want to retrieve genomic data from the CGD? (yes or no): ' GENEDATA

# Ask for raw data file format

echo 'Are the data in bam or fastq format?'
read -p 'Specify file format (bam or fastq): ' FORMAT

read -p 'Do you want to do a quality control of the raw data (yes or no): ' QCRAW

read -p 'Are the libraries prepared in a strand-specific way? (yes or no): ' STRANDED

read -p 'Are the data paired-end? (yes or no): ' PAIRED

read -p 'How many threads (cores) should be used for the analysis (use 1, if you are not sure): ' THREAD

if [ $GENEDATA == 'yes' ]
then
	echo 'Retrieving genomic data from CGD.'

	wget http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_chromosomes.fasta.gz  ## include WKDIR/required_files folder in wget!!!
	wget http://www.candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff
	cat C_albicans_SC5314_A22_current_features.gff | egrep -v "Ca22chr[1-7R]B" > C_albicans_SC5314_A22_current_features_haploid.gff
	gunzip C_albicans_SC5314_A22_current_chromosomes.fasta.gz
	cat C_albicans_SC5314_A22_current_chromosomes.fasta | egrep ">Ca22chr[1-7RM][A_]" | sed 's/>//g' | sed 's/(/1	/g' | sed 's/ nucleotides)//g' | sed 's/ /	/g' > chrMA.bed
	bedtools getfasta -fi C_albicans_SC5314_A22_current_chromosomes.fasta -bed chrMA.bed | fold -w 60 | sed 's/:1-[0-9]*//g' > C_albicans_SC5314_A22_current_chromosomesAM.fasta
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.gff.gz
	gunzip GCF_000182965.3_ASM18296v3_genomic.gff.gz ### contains Entrez ID to CGDID mappings
else
	echo 'No genomic data are retrieved.'
fi


GENOME=$WKDIR/required_files/C_albicans_SC5314_A22_current_chromosomesAM.fasta
FEATURES=$WKDIR/required_files/C_albicans_SC5314_A22_current_features_haploid.gff
ADAPT1=$(cat $WKDIR/required_files/config_file.txt | grep Read1: | cut -d ":" -f 2)
#ADAPT2=$(cat $WKDIR/required_files/config_file.txt | grep Read2: | cut -d ":" -f 2) # uncomment this for paired-end data
rRNA=$WKDIR/required_files/Ca_A22chrAM_rRNAloci.bed
mkdir $WKDIR/QC
PICARD=$WKDIR/required_files/picard.jar


# QC of raw data

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


# Adapter removal with cutadapt and mapping of all files with NGM

# for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(L2_1\.fq\.gz$)')
for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(R*_1\.fq\.gz$)')
do
	i1=$WKDIR/$SNAME
#	i2=$(echo $i1| sed 's/L2_1.fq.gz/L2_2.fq.gz/')
	i2=$(echo $i1| sed 's/_1.fq.gz/_2.fq.gz/')

	cutadapt --interleaved -j $THREAD -q 30 -O 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $i1.trimmed.fq.gz -p $i2.trimmed.fq.gz $i1 $i2  #2>$WKDIR/QC/Cutadapt_$SNAME.txt   # removes Illumina TrueSeq adapters from reads (change -a for different adapters), -j specifies number of cores to use, remove if not sure
	echo "it okey"
	#rm $i

	ngm -r $GENOME -1 $i1.trimmed.fq.gz -2 $i2.trimmed.fq.gz -o $i1.trimmed.fq.bam -b -p -Q 30 -t $THREAD # add -p for paired-end data, -t 6 is optional - means 6 threads of the processor are used, if you don't know what to do, remove it, --topn 1 --strata causes ngm to write only uniquely mapping reads to the output
	#rm $i.trimmed.fq.gz

	samtools sort -@ $THREAD $i1.trimmed.fq.bam -o $i1.trimmed.fq.bam.sort.bam   # sort .bam files using samtools
	mv $i1.count.txt $WKDIR/required_files
	#rm $i1.trimmed.fq.bam
	
	#bedtools intersect -a $i.trimmed.fq.bam.sort.bam -b $rRNA -v > $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam  # removal of reads mapping to rRNA loci
	bedtools intersect -a $i1.trimmed.fq.bam.sort.bam -b $rRNA -v > $i1.trimmed.fq.bam.sort.bam.rRNAfilt.bam  # removal of reads mapping to rRNA loci
	#rm $i.trimmed.fq.bam.sort.bam

	# Labelling of duplicated reads and removal of optical duplicates
	java -jar $PICARD MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true I=$i1.trimmed.fq.bam.sort.bam.rRNAfilt.bam O=$i1.final.bam M=$WKDIR/QC/$SNAME.markdup.metrics.txt
	#rm $i.trimmed.fq.bam.sort.bam.rRNAfilt.bam
	
	#Quality control and statistics about mapped samples
	samtools flagstat $i1.final.bam >> $WKDIR/QC/$SNAME.final.flagstat_analysis.txt   # flagstat analysis

	fastqc -o $WKDIR/QC $i1.final.bam

done

multiqc -s -o $WKDIR/QC $WKDIR/QC

# Preparation of coverage files for visualization in IGV

#mkdir $WKDIR/IGV_files

#for i1 in $WKDIR/*.final.bam
#do
#	samtools index $i1
#	SNAME=$(echo $i1 | sed 's:/.*/::g')
#	if [ $PAIRED == 'yes' ]
#	then
#		bamCoverage -b $i1 -o $WKDIR/IGV_files/$SNAME.bw --normalizeUsing CPM -bs 5 -p $THREAD -e
#	else
#		bamCoverage -b $i1 -o $WKDIR/IGV_files/$SNAME.bw --normalizeUsing CPM -bs 5 -p $THREAD
#	fi
#done



#mkdir $WKDIR/count
#mkdir $WKDIR/diff_expr_analysis

#for i1 in $WKDIR/*.final.bam
#do
	#if [ $PAIRED == 'yes' ]
	#then
	#	htseq-count -f bam -s $STRANDED -r pos -t gene -i ID $i1 $FEATURES > $i1.count.txt  # read count  for each gene with htseq-count
	#else
	#	htseq-count -f bam -s $STRANDED -t gene -i ID $i1 $FEATURES > $i1.count.txt
	#fi
	#mv $i1.count.txt $WKDIR/count
#done

#for i1 in $WKDIR/count/*.count.txt
#do
#	head -n -5 $i1 > $i1.crop.txt  # clear count files for flags
#done

#cp $WKDIR/count/*.crop.txt $WKDIR/diff_expr_analysis
#cp $FILES/edgeR_analysis.R $WKDIR/diff_expr_analysis
#cp $FILES/Targets.txt $WKDIR/diff_expr_analysis
