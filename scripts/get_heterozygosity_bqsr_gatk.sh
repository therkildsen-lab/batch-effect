#!/bin/bash
BAMLIST=/workdir/batch-effect/sample_lists/bam_list_realigned.txt
OUTDIR=/workdir/batch-effect/angsd/heterozygosity/
JOB_INDEX=0
JOBS=40
MINDP=2
MAXDP=10
MINQ=0
MINMAPQ=20
for LINE in `cat $BAMLIST`; do
  NAME_TEMP=`echo "${LINE%.*}"`
  NAME=`echo "${NAME_TEMP##*/}"`
	echo $NAME
  OUTBASE=$NAME'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ'_bqsr_gatk'
  
	/workdir/programs/angsd0.931/angsd/angsd \
  -i $NAME_TEMP'_bqsr.bam' \
  -anc /workdir/cod/reference_seqs/gadMor3.fasta \
  -out $OUTDIR$OUTBASE \
  -GL 1 \
  -doSaf 1 \
  -P 1 \
  -doCounts 1 \
  -setMinDepth $MINDP \
  -setMaxDepth $MAXDP \
  -minQ $MINQ \
  -minmapq $MINMAPQ &
  
  JOB_INDEX=$(( JOB_INDEX + 1 ))
	if [ $JOB_INDEX == $JOBS ]; then
		wait
		JOB_INDEX=0
	fi
done

wait

for LINE in `cat $BAMLIST`; do
  NAME_TEMP=`echo "${LINE%.*}"`
  NAME=`echo "${NAME_TEMP##*/}"`
	echo $NAME
  OUTBASE=$NAME'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ'_bqsr_gatk'
  
  /workdir/programs/angsd0.931/angsd/misc/realSFS \
  ${OUTDIR}${OUTBASE}.saf.idx \
  -P $JOBS \
  > ${OUTDIR}${OUTBASE}.ml
done