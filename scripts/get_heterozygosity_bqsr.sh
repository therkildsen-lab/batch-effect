#!/bin/bash
BATCH=$1
SAMPLESIZE=$2

OUTDIR=/workdir/batch-effect/angsd/heterozygosity/
JOB_INDEX=0
JOBS=10
MINDP=2
MAXDP=10
MINQ=0
MINMAPQ=30
KMAX=$(( SAMPLESIZE-1  ))
for K in $(seq 0 $KMAX); do
  LINE=`head /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_${BATCH}.txt -n $(( K+1 )) | tail -n 1`
  NAME_TEMP=`echo "${LINE%.*}"`
  NAME=`echo "${NAME_TEMP##*/}"`
	echo $NAME
  OUTBASE=$NAME'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ'_bqsr'
  
	/workdir/programs/angsd0.931/angsd/angsd \
  -i $LINE \
  -anc /workdir/cod/reference_seqs/gadMor3.fasta \
  -out $OUTDIR$OUTBASE \
  -GL 3 \
  -tmpdir /workdir/batch-effect/angsd_tmpdir_${BATCH}/subdir${K} \
  -doSaf 1 \
  -P 2 \
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

for K in $(seq 0 $KMAX); do
  LINE=`head /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_${BATCH}.txt -n $(( K+1 )) | tail -n 1`
  NAME_TEMP=`echo "${LINE%.*}"`
  NAME=`echo "${NAME_TEMP##*/}"`
	echo $NAME
  OUTBASE=$NAME'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ'_bqsr'
  
  /workdir/programs/angsd0.931/angsd/misc/realSFS \
  ${OUTDIR}${OUTBASE}.saf.idx \
  -P $JOBS \
  > ${OUTDIR}${OUTBASE}.ml
done