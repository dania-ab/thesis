for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(\.q\.gz$)')
    do
      i=$WKDIR/$SNAME
      fastqc -o $WKDIR/QC_raw $i
    done
