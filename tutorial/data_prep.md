Data Preparation
================

## Poly-G tail trimming

``` bash
zcat /workdir/cod/greenland-cod/adapter_clipped/QQL2011_884_55191_7_adapter_clipped_f_paired.fastq.gz | head -n 500000 | gzip > /workdir/batch-effect/tutorial/data/before_trimming_f.fastq.gz
zcat /workdir/cod/greenland-cod/adapter_clipped/QQL2011_884_55191_7_adapter_clipped_r_paired.fastq.gz | head -n 500000 | gzip > /workdir/batch-effect/tutorial/data/before_trimming_r.fastq.gz
```

## Base quality miscalibration

``` bash
## Come up with bam lists
rm /workdir/batch-effect/tutorial/data/bam_list_bq_se.txt
rm /workdir/batch-effect/tutorial/data/bam_list_bq_pe.txt

for LINE in `egrep 'IKE' /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_se.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  echo /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted_bq.bam >> /workdir/batch-effect/tutorial/data/bam_list_bq_se.txt
done

for LINE in `egrep 'IKE' /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_pe.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  echo /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted_bq.bam >> /workdir/batch-effect/tutorial/data/bam_list_bq_pe.txt
done

cp /workdir/batch-effect/tutorial/data/bam_list_bq_se.txt /workdir/batch-effect/tutorial/data/bam_list_bq.txt
cat /workdir/batch-effect/tutorial/data/bam_list_bq_pe.txt >> /workdir/batch-effect/tutorial/data/bam_list_bq.txt

## Come up with a sample table
head -n 1 /workdir/batch-effect/sample_lists/sample_table_merged.tsv >  /workdir/batch-effect/tutorial/data/sample_table_bq.tsv
egrep 'IKE' /workdir/batch-effect/sample_lists/sample_table_merged.tsv >>  /workdir/batch-effect/tutorial/data/sample_table_bq.tsv

## Subset bam files
for LINE in `egrep 'IKE' /workdir/batch-effect/sample_lists/bam_list_realigned.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  samtools view -b $LINE "LG03" > /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted_bq.bam
done
```

## Reference bias

``` bash
rm /workdir/batch-effect/tutorial/data/bam_list_rb_se.txt
rm /workdir/batch-effect/tutorial/data/bam_list_rb_pe.txt
for LINE in `egrep 'KNG|QQL|BUK|PAA' /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_se.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  echo /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted_rb.bam >> /workdir/batch-effect/tutorial/data/bam_list_rb_se.txt
done

for LINE in `egrep 'KNG|QQL|BUK|PAA' /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_pe.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  echo /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted_rb.bam >> /workdir/batch-effect/tutorial/data/bam_list_rb_pe.txt
done

cp /workdir/batch-effect/tutorial/data/bam_list_rb_se.txt /workdir/batch-effect/tutorial/data/bam_list_rb.txt
cat /workdir/batch-effect/tutorial/data/bam_list_rb_pe.txt >> /workdir/batch-effect/tutorial/data/bam_list_rb.txt

## Come up with a sample table
head -n 1 /workdir/batch-effect/sample_lists/sample_table_merged.tsv >  /workdir/batch-effect/tutorial/data/sample_table_rb.tsv
egrep 'KNG|QQL|BUK|PAA' /workdir/batch-effect/sample_lists/sample_table_merged.tsv >>  /workdir/batch-effect/tutorial/data/sample_table_rb.tsv

## Subset bam files
for LINE in `egrep 'KNG|QQL|BUK|PAA' /workdir/batch-effect/sample_lists/bam_list_realigned.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  samtools view -b $LINE "LG06:6000000-6500000" > /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted_rb.bam
done
```

## DNA degradation

``` bash
## Come up with bam lists
for LINE in `egrep 'UUM' /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_se.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  echo /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted.bam >> /workdir/batch-effect/tutorial/data/bam_list_degradation_se.txt
done

for LINE in `egrep 'UUM' /workdir/batch-effect/sample_lists/bam_list_per_pop/bam_list_realigned_pe.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  echo /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted.bam >> /workdir/batch-effect/tutorial/data/bam_list_degradation_pe.txt
done

cp /workdir/batch-effect/tutorial/data/bam_list_degradation_se.txt /workdir/batch-effect/tutorial/data/bam_list_degradation.txt
cat /workdir/batch-effect/tutorial/data/bam_list_degradation_pe.txt >> /workdir/batch-effect/tutorial/data/bam_list_degradation.txt

## Come up with a sample table
head -n 1 /workdir/batch-effect/sample_lists/sample_table_merged.tsv >  /workdir/batch-effect/tutorial/data/sample_table_degradation.tsv
egrep 'UUM' /workdir/batch-effect/sample_lists/sample_table_merged.tsv >>  /workdir/batch-effect/tutorial/data/sample_table_degradation.tsv

## Subset bam files
for LINE in `egrep 'UUM' /workdir/batch-effect/sample_lists/bam_list_realigned.txt`; do
  TEMP=${LINE%.*}
  PREFIX=${TEMP/*\/}
  samtools view -b $LINE "LG03:5000000-15000000" > /workdir/batch-effect/tutorial/data/${PREFIX}_subsetted.bam
done
```
