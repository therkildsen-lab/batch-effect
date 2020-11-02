BASEDIR=/workdir/batch-effect/
INPUTDIR=/workdir/cod/greenland-cod/adapter_clipped/
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_pe.txt
SAMPLETABLE=$BASEDIR/sample_lists/sample_table_unmerged.tsv

for SAMPLEFILE in `cat $SAMPLELIST`; do
  SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
  SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
  LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
  SAMPLE_SEQ_ID=$SAMPLE_ID"_"$SEQ_ID"_"$LANE_ID
  
  ## Extract data type from the sample table
  DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
  
  ## The input and output path and file prefix
  SAMPLEADAPT=$INPUTDIR$SAMPLE_SEQ_ID
  SAMPLEQUAL=$BASEDIR"polyg_trimmed/"$SAMPLE_SEQ_ID

  /workdir/programs/fastp --trim_poly_g -L -A \
  -i $SAMPLEADAPT"_adapter_clipped_f_paired.fastq.gz" \
  -I $SAMPLEADAPT"_adapter_clipped_r_paired.fastq.gz" \
  -o $SAMPLEQUAL"_adapter_clipped_qual_filtered_f_paired.fastq.gz" \
  -O $SAMPLEQUAL"_adapter_clipped_qual_filtered_r_paired.fastq.gz" \
  -h $SAMPLEQUAL"_adapter_clipped_qual_filtered_fastp.html"
done
