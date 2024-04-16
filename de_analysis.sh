
FILES=$(pwd)
WKDIR=$(echo $FILES | sed 's:/required_files::g')

cp $WKDIR/count/*.crop.txt $WKDIR/diff_expr_analysis
cp $FILES/edgeR_analysis.R $WKDIR/diff_expr_analysis
cp $FILES/Targets.txt $WKDIR/diff_expr_analysis

for i1 in $WKDIR/*.final.bam
do
samtools index $i1
SNAME=$(echo $i1 | sed 's:/.*/::g')
bamCoverage -b $i1 -o $WKDIR/IGV_files/$SNAME.bw --normalizeUsing CPM -p 15 -e
done
