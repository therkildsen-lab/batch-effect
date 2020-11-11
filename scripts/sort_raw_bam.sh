OUTPUTPATH='/workdir/batch-effect/bam/'

for K in {2..168}; do
  SAMPLE=`head /workdir/batch-effect/sample_lists/sample_table_unmerged.tsv -n $K | tail -n 1 | cut -f 4`
  LANE=`head /workdir/batch-effect/sample_lists/sample_table_unmerged.tsv -n $K | tail -n 1 | cut -f 2`
  SEQID=`head /workdir/batch-effect/sample_lists/sample_table_unmerged.tsv -n $K | tail -n 1 | cut -f 3`
  DATATYPE=`head /workdir/batch-effect/sample_lists/sample_table_unmerged.tsv -n $K | tail -n 1 | cut -f 6`
  
  PREFIX=$SAMPLE'_'$SEQID'_'$LANE'_'$DATATYPE'_bt2_gadMor3'
  
  if [ $DATATYPE = se ]; then
    INPUTPATH='/workdir/cod/greenland-cod/bam/'
  else
    INPUTPATH='/workdir/batch-effect/bam/'
  fi
  #echo $INPUTPATH$PREFIX'.bam'
  #echo $OUTPUTPATH$PREFIX'_sorted.bam'
  samtools sort $INPUTPATH$PREFIX'.bam' -o $OUTPUTPATH$PREFIX'_sorted.bam' -@ 24
done
