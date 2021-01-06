N_CORE_MAX=10
COUNT=0
for FILE in `cat /workdir/cod/greenland-cod/sample_lists/bam_list_realigned_mincov_contamination_filtered_batch_effect.txt`; do
PREFIX=${FILE%.bam}
echo $PREFIX
samtools view -h $FILE | grep -v XS:i | samtools view -buS - | samtools sort -o $PREFIX"_uniqueonly.bam" &
COUNT=$(( COUNT + 1 ))
if [ $COUNT == $N_CORE_MAX ]; then
  wait
  COUNT=0
fi
done
