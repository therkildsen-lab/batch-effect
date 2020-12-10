Core data processing and analysis pipeline
================

## Load packages

``` r
library(tidyverse)
library(cowplot)
```

## Come up with some sample lists and tables

``` r
base_dir = "/workdir/batch-effect/"
# Read in the full sample table from the Greenland cod project
sample_table_full <- read_tsv("../../cod/greenland-cod/sample_lists/sample_table_merged_mincov_contamination_filtered.tsv") %>%
  bind_cols(bam_list = read_lines("../../cod/greenland-cod/sample_lists/bam_list_realigned_mincov_contamination_filtered.txt"))
# Select a subset of populations that were sequenced in both Hiseq and Nextseq platforms
# Note that there are some QQL samples that appear to be "merged", but they were merged from lane 1 and 2
sample_table_merged <- filter(sample_table_full, 
                                    population %in% c("IKE2011", "QQL2011", "ITV2011", "KNG2011", "BUK2011", "NAR2008", "UUM2010", "PAA2011", "ATP2011")) %>%
  filter(data_type != "pese")
## Raw bam file list for the pair end samples
bam_list_merged_pe <- sample_table_merged %>%
  filter(data_type=="pe") %>%
  mutate(bam_list = str_c(base_dir, "bam/", sample_seq_id, "_bt2_gadMor3_sorted.bam")) %>%
  dplyr::select(bam_list)
## Raw bam file list for the single end samples
bam_list_merged_se <- sample_table_merged %>%
  filter(data_type=="se") %>%
  mutate(bam_list = str_c(base_dir, "bam/", sample_seq_id, "_bt2_gadMor3_sorted.bam")) %>%
  dplyr::select(bam_list)
## Overlap clipped bamlist for all samples
bam_list_overlap_clipped <- sample_table_merged %>%
  mutate(bam_list = ifelse(data_type == "pe", 
                           str_c(base_dir, "bam/", sample_seq_id, "_bt2_gadMor3_sorted_dedup_overlapclipped.bam"),
                           str_c(base_dir, "bam/", sample_seq_id, "_bt2_gadMor3_sorted_dedup.bam"))) %>%
  dplyr::select(bam_list)
## Indel realigned bamlist for all samples
bam_list_realigned <- sample_table_merged %>%
  mutate(bam_list = ifelse(data_type == "pe", 
                           str_c(base_dir, "bam/", sample_seq_id, "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam"),
                           str_c(base_dir, "bam/", sample_seq_id, "_bt2_gadMor3_sorted_dedup_realigned.bam"))) %>%
  dplyr::select(bam_list, data_type)

## Get unmerged sample table for the select subset of samples
sample_table_unmerged <- read_tsv("../../cod/greenland-cod/sample_lists/sample_table.tsv") %>%
  semi_join(sample_table_merged, by=c("sample_id"="sample_id_corrected"))
## Fastq list of pe and se samples
fastq_list_pe <- filter(sample_table_unmerged, lane_number == 7)$prefix
fastq_list_se <- filter(sample_table_unmerged, lane_number != 7)$prefix

## Sample distribution
sample_table_merged %>%
  ggplot(aes(x=population_new, fill=data_type)) +
  geom_bar(color="black") +
  theme_cowplot() +
  coord_flip()
```

![](pipeline_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Write the objects created above
sample_table_merged %>%
  dplyr::select(-bam_list) %>%
  write_tsv("../sample_lists/sample_table_merged.tsv")
write_tsv(bam_list_merged_pe, "../sample_lists/bam_list_merged_pe.txt", col_names = F)
write_tsv(bam_list_merged_se, "../sample_lists/bam_list_merged_se.txt", col_names = F)
write_tsv(bam_list_overlap_clipped, "../sample_lists/bam_list_dedup_overlapclipped.txt", col_names = F)
bam_list_realigned %>%
  dplyr::select(bam_list) %>%
  write_tsv("../sample_lists/bam_list_realigned.txt", col_names = F)
bam_list_realigned %>%
  filter(data_type=="pe") %>%
  dplyr::select(bam_list) %>%
  write_tsv("../sample_lists/bam_list_per_pop/bam_list_realigned_pe.txt", col_names = F)
bam_list_realigned %>%
  filter(data_type=="se") %>%
  dplyr::select(bam_list) %>%
  write_tsv("../sample_lists/bam_list_per_pop/bam_list_realigned_se.txt", col_names = F)
write_tsv(sample_table_unmerged, "../sample_lists/sample_table_unmerged.tsv")
write_lines(fastq_list_pe, "../sample_lists/fastq_list_pe.txt")
write_lines(fastq_list_se, "../sample_lists/fastq_list_se.txt")
```

## Process the PE data

#### polyG trimming

This is run for comparison’s sake. The resulting fastq files will not be
used in downstream analyses.

``` bash
echo 'BASEDIR=/workdir/batch-effect/
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
done' > /workdir/batch-effect/scripts/trim_polyg_batch_effect.sh

nohup bash /workdir/batch-effect/scripts/trim_polyg_batch_effect.sh > /workdir/batch-effect/nohups/cut_right_batch_effect.nohups &
```

#### Sliding window trimming

``` bash
echo 'SAMPLELIST=/workdir/cod/greenland-cod/sample_lists/sample_list_pe_1.tsv 
SAMPLETABLE=/workdir/cod/greenland-cod/sample_lists/sample_table_pe.tsv
BASEDIR=/workdir/cod/greenland-cod/

for SAMPLEFILE in `cat $SAMPLELIST`; do
  SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
  SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
  LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
  SAMPLE_SEQ_ID=$SAMPLE_ID"_"$SEQ_ID"_"$LANE_ID
  
  ## Extract data type from the sample table
  DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
  
  ## The input and output path and file prefix
  SAMPLEADAPT=$BASEDIR"adapter_clipped/"$SAMPLE_SEQ_ID
  SAMPLEQUAL=$BASEDIR"qual_filtered/"$SAMPLE_SEQ_ID

  /workdir/programs/fastp --trim_poly_g -L -A --cut_right \
  -i $SAMPLEADAPT"_adapter_clipped_f_paired.fastq.gz" \
  -I $SAMPLEADAPT"_adapter_clipped_r_paired.fastq.gz" \
  -o $SAMPLEQUAL"_adapter_clipped_qual_filtered_f_paired.fastq.gz" \
  -O $SAMPLEQUAL"_adapter_clipped_qual_filtered_r_paired.fastq.gz" \
  -h $SAMPLEQUAL"_adapter_clipped_qual_filtered_fastp.html"
done' > /workdir/cod/greenland-cod/scripts/cut_right_batch_effect.sh

nohup bash /workdir/cod/greenland-cod/scripts/cut_right_batch_effect.sh > /workdir/cod/greenland-cod/nohups/cut_right_batch_effect.nohups &
mv /workdir/cod/greenland-cod/qual_filtered/* /workdir/batch_effect/qual_filtered/
```

#### Map to reference genome

This step also includes quality filtering and sorting, but I will rerun
the quality filtering and sorting step with no quality filter, both for
PE and SE data.

``` bash
nohup bash /workdir/data-processing/scripts/low_coverage_mapping.sh \
/workdir/batch-effect/sample_lists/fastq_list_pe.txt \
/workdir/batch-effect/sample_lists/sample_table_unmerged.tsv \
/workdir/batch-effect/qual_filtered/ \
/workdir/batch-effect/ \
_adapter_clipped_qual_filtered_f_paired.fastq.gz \
_adapter_clipped_qual_filtered_r_paired.fastq.gz \
very-sensitive \
/workdir/cod/reference_seqs/gadMor3.fasta \
gadMor3 \
> /workdir/batch-effect/nohups/low_coverage_mapping_pe.nohup &
```

## Process the PE and SE data together

These steps are necessary for the SE data becasue intermediate bam files
were deleted.

#### Sort the raw bam files

Save the following script as
`/workdir/batch-effect/scripts/sort_raw_bam.sh`

``` bash
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
```

``` bash
nohup bash /workdir/batch-effect/scripts/sort_raw_bam.sh > /workdir/batch-effect/nohups/sort_raw_bam.nohup &
```

#### Merge four samples that were sequenced in multiple lanes

``` bash
## Merge four samples that were sequenced in multiple lanes
samtools merge /workdir/batch-effect/bam/QQL2011_842_merged_merged_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_842_14247X225_1_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_842_14247X47_2_se_bt2_gadMor3_sorted.bam -@ 24
samtools merge /workdir/batch-effect/bam/QQL2011_844_merged_merged_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_844_14247X226_1_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_844_14247X49_2_se_bt2_gadMor3_sorted.bam -@ 24
samtools merge /workdir/batch-effect/bam/QQL2011_846_merged_merged_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_846_14247X232_1_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_846_14247X61_2_se_bt2_gadMor3_sorted.bam -@ 24
samtools merge /workdir/batch-effect/bam/QQL2011_852_merged_merged_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_852_14247X228_1_se_bt2_gadMor3_sorted.bam /workdir/batch-effect/bam/QQL2011_852_14247X40_2_se_bt2_gadMor3_sorted.bam -@ 24
```

## Deduplicate

I slightly modified the original script
(`/workdir/data-processing/scripts/deduplicate_clipoverlap.sh`) so that
the `minq20` part is no longer in the file name.

``` bash
## SE
nohup bash /workdir/batch-effect/scripts/deduplicate_clipoverlap.sh \
/workdir/batch-effect/sample_lists/bam_list_merged_se.txt \
/workdir/batch-effect/sample_lists/sample_table_merged.tsv \
/workdir/batch-effect/ \
gadMor3 \
> /workdir/batch-effect/nohups/deduplicate_se.nohup &
## PE
nohup bash /workdir/batch-effect/scripts/deduplicate_clipoverlap.sh \
/workdir/batch-effect/sample_lists/bam_list_merged_pe.txt \
/workdir/batch-effect/sample_lists/sample_table_merged.tsv \
/workdir/batch-effect/ \
gadMor3 \
> /workdir/batch-effect/nohups/deduplicate_clipoverlap_pe.nohup &
```

## Indel realignment with PE and SE data combined

``` bash
cp /workdir/batch-effect/sample_lists/bam_list_dedup_overlapclipped.txt \
/workdir/batch-effect/sample_lists/bam_list_dedup_overlapclipped.list

nohup bash /workdir/data-processing/scripts/realign_indels.sh \
/workdir/batch-effect/sample_lists/bam_list_dedup_overlapclipped.list \
/workdir/batch-effect/ \
/workdir/cod/reference_seqs/gadMor3.fasta \
gadMor3 \
> /workdir/batch-effect/nohups/realign_indels.nohup &
```

## Get genotype likelihoods

#### polyG trimmed PE samples (original)

This is with PE samples that went through polyG trimming instead of
sliding window trimming.

``` bash
## original setup
cd /workdir/cod/greenland-cod/
/workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_realigned_mincov_contamination_filtered_batch_effect.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/bam_list_realigned_mincov_contamination_filtered_batch_effect_mindp2_maxdp661_minind2_minq20 \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 32 -setMinDepth 2 -setMaxDepth 661 -minInd 2 -minQ 20 -minMaf 0.01 \
-sites angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.txt \
-rf angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.chrs \
>& nohups/get_gl_bam_list_realigned_mincov_contamination_filtered_batch_effect.log &

## more stringent filtering
cd /workdir/cod/greenland-cod/
/workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_realigned_mincov_contamination_filtered_batch_effect.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/bam_list_realigned_mincov_contamination_filtered_batch_effect_mindp2_maxdp661_minind2_stringent_filter \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 32 -setMinDepth 2 -setMaxDepth 661 -minInd 2 -minMaf 0.01 \
-minQ 33 -minMapQ 25 -uniqueOnly 1 \
-sites angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.txt \
-rf angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.chrs \
>& nohups/get_gl_bam_list_realigned_mincov_contamination_filtered_batch_effect_stringent_filter.log &
```

#### Sliding window trimmed PE samples (new)

This is with PE samples that went through sliding window trimming.

I will first try the same LD-trimmed SNP list from the Greenland cod
project

``` bash
cd /workdir/batch-effect/
nohup /workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_realigned.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/bam_list_realigned \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 16 -setMinDepth 2 -setMaxDepth 661 -minInd 2 -minQ 20 -minMaf 0.01 \
-doIBS 1 -makematrix 1 -doCov 1 \
-sites /workdir/cod/greenland-cod/angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.txt \
-rf /workdir/cod/greenland-cod/angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.chrs \
>& nohups/get_gl_bam_list_realigned.log &
```

#### Sliding window trimmed PE samples (new)

``` bash
nohup python2 /workdir/programs/pcangsd/pcangsd.py \
-beagle /workdir/batch-effect/angsd/bam_list_realigned.beagle.gz \
-minMaf 0.05 \
-threads 8 \
-o /workdir/batch-effect/angsd/bam_list_realigned_pcangsd \
> /workdir/batch-effect/nohups/run_pcangsd.nohup &
```