
FILES=$(pwd)
WKDIR=$(echo $FILES | sed 's:/required_files::g')

cp $WKDIR/count/*.crop.txt $WKDIR/diff_expr_analysis
cp $FILES/edgeR_analysis.R $WKDIR/diff_expr_analysis
cp $FILES/Targets.txt $WKDIR/diff_expr_analysis
