#!/bin/bash
BAMLIST=/workdir/batch-effect/sample_lists/bam_list_realigned.txt
BASEDIR=/workdir/batch-effect/
REFERENCE=/workdir/cod/reference_seqs/gadMor3.fasta
VCF=/workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mincov_filtered_mindp249_maxdp1142_minind111_minq20.vcf

JOB_INDEX=0
JOBS=20

for LINE in `cat $BAMLIST`; do
  NAME_TEMP=`echo "${LINE%.*}"`
  NAME=`echo "${NAME_TEMP##*/}"`
	echo $NAME

  /programs/gatk-4.2.0.0/gatk BaseRecalibrator \
  -I $LINE \
  -R $REFERENCE \
  --known-sites $VCF \
  -O ${BASEDIR}/bam/${NAME}_recal_data.table &
  
  JOB_INDEX=$(( JOB_INDEX + 1 ))
	if [ $JOB_INDEX == $JOBS ]; then
		wait
		JOB_INDEX=0
	fi
done

wait

JOB_INDEX=0
for LINE in `cat $BAMLIST`; do
  NAME_TEMP=`echo "${LINE%.*}"`
  NAME=`echo "${NAME_TEMP##*/}"`
	echo $NAME
  /programs/gatk-4.2.0.0/gatk ApplyBQSR \
  -R $REFERENCE \
  -I $LINE \
  --bqsr-recal-file ${BASEDIR}/bam/${NAME}_recal_data.table \
  -O ${BASEDIR}/bam/${NAME}_bqsr.bam &
  
  JOB_INDEX=$(( JOB_INDEX + 1 ))
	if [ $JOB_INDEX == $JOBS ]; then
		wait
		JOB_INDEX=0
	fi
done