#!/bin/bash
BAMLIST=/workdir/batch-effect/sample_lists/bam_list_realigned.txt
BASEDIR=/workdir/batch-effect/
REFERENCE=/workdir/cod/reference_seqs/gadMor3.fasta
VCF=/workdir/batch-effect/angsd/global_snp_list_bqsr.vcf

JOB_INDEX=0
JOBS=20

for LINE in `cat $BAMLIST`; do
  NAME_TEMP=`echo "${LINE%.*}"`
  NAME=`echo "${NAME_TEMP##*/}"`
	echo $NAME
  /programs/gatk-4.2.0.0/gatk ApplyBQSR \
  -R $REFERENCE \
  -I $LINE \
  --bqsr-recal-file ${BASEDIR}/bam/recal_data_all_samples.table \
  -O ${BASEDIR}/bam/${NAME}_bqsr.bam &
  
  JOB_INDEX=$(( JOB_INDEX + 1 ))
	if [ $JOB_INDEX == $JOBS ]; then
		wait
		JOB_INDEX=0
	fi
done