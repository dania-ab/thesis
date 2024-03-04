# Preparation of coverage files for visualization in IGV
FILES=$(pwd)
WKDIR=$(echo $FILES | sed 's:/required_files::g')

#for i1 in $WKDIR/*.final.bam
#do
#	samtools index $i1
#	SNAME=$(echo $i1 | sed 's:/.*/::g')
	
#	bamCoverage -b $i1 -o $WKDIR/IGV_files/$SNAME.bw --normalizeUsing CPM -bs 5 -p 15 -e
#done


#mkdir $WKDIR/count
#mkdir $WKDIR/diff_expr_analysis

#for i1 in $WKDIR/*.final.bam
#do
#	htseq-count -f bam -s no -r pos -t gene -i ID $i1 /Users/u0120981/Desktop/RNA/required_files/C_albicans_SC5314_A22_current_features_haploid.gff > $i1.count.txt  # read count  for each gene with htseq-count
#	mv $i1.count.txt $WKDIR/count
#done

for i1 in $WKDIR/count/*.count.txt
do
	echo "loop $i1"
	head -n 5 $i1 > $i1.crop.txt  # clear count files for flags
done

cp $WKDIR/count/*.crop.txt $WKDIR/diff_expr_analysis
cp $FILES/edgeR_analysis.R $WKDIR/diff_expr_analysis
cp $FILES/Targets.txt $WKDIR/diff_expr_analysis
