Poly-G tail
================

  - [Load packages](#load-packages)
  - [Read in the sample table](#read-in-the-sample-table)
  - [G content comparisons](#g-content-comparisons)
      - [Without trimming](#without-trimming)
      - [Poly-G trimming with fastp in default setting
        (`--trim_poly_g`)](#poly-g-trimming-with-fastp-in-default-setting---trim_poly_g)
      - [Sliding window trimming with fastp in default setting
        (`--cut_right`)](#sliding-window-trimming-with-fastp-in-default-setting---cut_right)
      - [Base content in three samples](#base-content-in-three-samples)
      - [Example of one read](#example-of-one-read)
      - [Combined plot](#combined-plot)
      - [Persistence of poly-G tails in bam
        files](#persistence-of-poly-g-tails-in-bam-files)
  - [Pipeline after sliding-window
    trimming](#pipeline-after-sliding-window-trimming)
      - [Map to reference genome](#map-to-reference-genome)
      - [Sort the raw bam files](#sort-the-raw-bam-files)
      - [Merge four samples that were sequenced in multiple
        lanes](#merge-four-samples-that-were-sequenced-in-multiple-lanes)
      - [Deduplicate](#deduplicate)
      - [Indel realignment](#indel-realignment)
      - [SNP calling](#snp-calling)
      - [SAF, MAF, Fst estimation, and read depth count in each batch of
        data](#saf-maf-fst-estimation-and-read-depth-count-in-each-batch-of-data)
      - [LD pruning](#ld-pruning)
      - [Get the covariance matrix with ANGSD using LD pruned SNP
        list](#get-the-covariance-matrix-with-angsd-using-ld-pruned-snp-list)

## Load packages

``` r
library(tidyverse)
library(cowplot)
```

## Read in the sample table

``` r
sample_table_unmerged_pe <- read_tsv("../sample_lists/sample_table_unmerged.tsv") %>%
  filter(lane_number==7) %>%
  mutate(prefix_new=str_c(sample_id, "_", seq_id, "_", lane_number))
```

## G content comparisons

#### Without trimming

``` bash
## Run FastQC on three random samples after adapter trimming but before polyG or sliding window trimming
fastqc /workdir/cod/greenland-cod/adapter_clipped/QQL2011_884_55191_7_adapter_clipped_f_paired.fastq.gz
fastqc /workdir/cod/greenland-cod/adapter_clipped/UUM2010_068_55111_7_adapter_clipped_f_paired.fastq.gz
fastqc /workdir/cod/greenland-cod/adapter_clipped/IKE2011_976_55124_7_adapter_clipped_f_paired.fastq.gz
```

#### Poly-G trimming with fastp in default setting (`--trim_poly_g`)

This step is the same as in the original pipeline that we used for the
Greenland project, but I had to rerun this for a subset of samples
because the intermediate files for the Greenland project were deleted.

``` bash
## Run polyg trimming
echo 'BASEDIR=/workdir/batch-effect/
INPUTDIR=/workdir/cod/greenland-cod/adapter_clipped/
SAMPLELIST=$BASEDIR/sample_lists/fastq_list_pe.txt
SAMPLETABLE=$BASEDIR/sample_lists/sample_table_unmerged.tsv

for SAMPLEFILE in `cat $SAMPLELIST`; do
  SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
  SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
  LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
  SAMPLE_SEQ_ID=$SAMPLE_ID"_"$SEQ_ID"_"$LANE_ID
  
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

## Run fastqc
nohup fastqc /workdir/batch-effect/polyg_trimmed/*fastq.gz -t 10 > /workdir/batch-effect/run_fastqc_polyg_trimmed.nohup &
```

#### Sliding window trimming with fastp in default setting (`--cut_right`)

``` bash
## Run sliding window trimming
echo 'SAMPLELIST=/workdir/cod/greenland-cod/sample_lists/sample_list_pe_1.tsv 
SAMPLETABLE=/workdir/cod/greenland-cod/sample_lists/sample_table_pe.tsv
BASEDIR=/workdir/batch-effect/
INPUTDIR=/workdir/cod/greenland-cod/adapter_clipped/

for SAMPLEFILE in `cat $SAMPLELIST`; do
  SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
  SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
  LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
  SAMPLE_SEQ_ID=$SAMPLE_ID"_"$SEQ_ID"_"$LANE_ID
  
  ## The input and output path and file prefix
  SAMPLEADAPT=$INPUTDIR$SAMPLE_SEQ_ID
  SAMPLEQUAL=$BASEDIR"cut_right/"$SAMPLE_SEQ_ID

  /workdir/programs/fastp --trim_poly_g -L -A --cut_right \
  -i $SAMPLEADAPT"_adapter_clipped_f_paired.fastq.gz" \
  -I $SAMPLEADAPT"_adapter_clipped_r_paired.fastq.gz" \
  -o $SAMPLEQUAL"_adapter_clipped_qual_filtered_f_paired.fastq.gz" \
  -O $SAMPLEQUAL"_adapter_clipped_qual_filtered_r_paired.fastq.gz" \
  -h $SAMPLEQUAL"_adapter_clipped_qual_filtered_fastp.html"
done' > /workdir/batch-effect/scripts/cut_right_batch_effect.sh

nohup bash /workdir/batch-effect/scripts/cut_right_batch_effect.sh > /workdir/batch-effect/nohups/cut_right_batch_effect.nohups &

## run fastqc with three random samples
fastqc /workdir/batch-effect/qual_filtered/QQL2011_884_55191_7_adapter_clipped_qual_filtered_f_paired.fastq.gz
fastqc /workdir/batch-effect/qual_filtered/UUM2010_068_55111_7_adapter_clipped_qual_filtered_f_paired.fastq.gz
fastqc /workdir/batch-effect/qual_filtered/IKE2011_976_55124_7_adapter_clipped_qual_filtered_f_paired.fastq.gz
```

#### Base content in three samples

``` r
random_sample <- filter(sample_table_unmerged_pe, sample_id %in% c("QQL2011_884", "UUM2010_068", "IKE2011_976"))
paths <- c("../../cod/greenland-cod/adapter_clipped/", "../polyg_trimmed/", "../cut_right/")
suffice <- c("_adapter_clipped_f_paired_fastqc", "_adapter_clipped_qual_filtered_f_paired_fastqc", "_adapter_clipped_qual_filtered_f_paired_fastqc")
types <- c("no trimming", "poly-G trimming", "sliding-window trimming")
for (j in 1:3){
  path <- paths[j]
  suffix <- suffice[j]
  type <- types[j]
  for (i in 1:3){
    sample_id <- random_sample$sample_id[i]
    sample_prefix <- random_sample$prefix_new[i]
    file_name <- str_c(path, sample_prefix, suffix)
    unzip(str_c(file_name, ".zip"), exdir = path, overwrite = FALSE)
    fastqc_data <- read_lines(file = str_c(file_name, "/fastqc_data.txt"))
    first_line <- which(str_detect(fastqc_data, ">>Per base sequence content")) + 1
    last_line <- which(str_detect(fastqc_data, ">>Per sequence GC content")) - 2
    per_base_seq_content_polyg_trimmed <- fastqc_data[first_line:last_line] %>%
      read_tsv() %>%
      rename(position=`#Base`) %>%
      pivot_longer(2:5, names_to = "base", values_to = "percentage") %>%
      mutate(sample_id=sample_id, type = type)
    if (j == 1 & i == 1) {
      per_base_seq_content_polyg_trimmed_final <- per_base_seq_content_polyg_trimmed
    } else {
      per_base_seq_content_polyg_trimmed_final <- bind_rows(per_base_seq_content_polyg_trimmed_final, per_base_seq_content_polyg_trimmed)
    }
  }
}

seq_content_p <- per_base_seq_content_polyg_trimmed_final %>%
  mutate(position = as_factor(position)) %>%
  ggplot(aes(x=position, y=percentage, color=base, group=base)) +
  geom_line(size=0.8) +
  scale_color_manual(values = c("#749dae", "#5445b1", "orange", "#cd3341")) +
  xlab("read position (in bp)") +
  facet_grid(str_c("ind ", as.numeric(as.factor(sample_id)))~type) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
seq_content_p
```

![](polyg_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### Example of one read

###### Extract the read and trim it

``` bash
## Save the read as a separate file
zcat /workdir/backup/cod/greenland_cod/fastq/8467_3270_55085_HMK3YBGX2_NAR2008_006_GCTACGCT_CTCTCTAT_R1.fastq.gz | grep ATCCCGCACCCTCCCATTTCTCTTCAACAACAACAACCTCCGCCGCCCATCCCGTGTCACACACGGGCGCGCGGGGGGGGGGGGGGGGGGTGGCGCGGGGC -A2 -B1 > /workdir/batch-effect/misc/polyg_example.fastq
## Trim polyG tail
/workdir/programs/fastp --trim_poly_g -L -A -i /workdir/batch-effect/misc/polyg_example.fastq -o /workdir/batch-effect/misc/polyg_example_trim_polyg.fastq
## Window trim
/workdir/programs/fastp --trim_poly_g -L -A --cut_right -i /workdir/batch-effect/misc/polyg_example.fastq -o /workdir/batch-effect/misc/polyg_example_polyg_cut_right.fastq
```

###### Visualization

``` r
for (file in rev(c("../misc/polyg_example.fastq", "../misc/polyg_example_trim_polyg.fastq", "../misc/polyg_example_polyg_cut_right.fastq"))){
  sequence_vector <- read_lines(file)[2] %>%
    str_split(pattern = "") %>%
    .[[1]]
  base_p <- tibble(base=sequence_vector) %>%
    {mutate(., id=seq(nrow(.)))} %>%
    ggplot(aes(x=id, y=0, fill=base)) +
    geom_tile(color="black", size=0.3) +
    #scale_fill_viridis_d(begin = 0.15) +
    scale_fill_manual(values = c("#749dae", "#5445b1", "orange", "#cd3341")) +
    theme_void() +
    theme(legend.position = "top")
  print(base_p)
}
```

![](polyg_files/figure-gfm/unnamed-chunk-8-1.svg)<!-- -->![](polyg_files/figure-gfm/unnamed-chunk-8-2.svg)<!-- -->![](polyg_files/figure-gfm/unnamed-chunk-8-3.svg)<!-- -->

``` r
base_p_final <- base_p +
  coord_cartesian(clip = 'off') +
  annotate("text", 48.5, 2.2, label="sliding-window\ntrimming", size = 3.8) +
  annotate("segment", 48.5, 1.3, xend=48.5, yend=0.6, arrow=arrow(length = unit(0.1, "npc")), size=1) +
  annotate("text", 101.5, 2.2, label="poly-G\ntrimming", size = 3.8) +
  annotate("segment", 101.5, 1.3, xend=101.5, yend=0.6, arrow=arrow(length = unit(0.1, "npc")), size=1) #+
  #annotate("text", 76, 2.5, label='bold("base")', size = 3.8, parse=TRUE) +
  #theme(legend.title=element_blank())

base_p_final
```

![](polyg_files/figure-gfm/unnamed-chunk-9-1.svg)<!-- -->

``` r
for (file in rev(c("../misc/polyg_example.fastq", "../misc/polyg_example_trim_polyg.fastq", "../misc/polyg_example_polyg_cut_right.fastq"))){
  quality_vector <- read_lines(file)[4] %>%
    charToRaw() %>%
    as.integer() %>% 
    {.-33}
  quality_p <- tibble(quality=quality_vector) %>%
    {mutate(., id=seq(nrow(.)))} %>%
    ggplot(aes(x=id, y=0, fill=quality)) +
    geom_tile(color="black", size=0.3) +
    scale_fill_continuous(high = "darkgreen", low = "lightgreen") +
    #scale_fill_viridis_c(option = "D", begin = 0.4, direction = -1) +
    theme_void() +
    labs(fill = "quality score") +
    theme(legend.position = "bottom")
  print(quality_p)
}
```

![](polyg_files/figure-gfm/unnamed-chunk-10-1.svg)<!-- -->![](polyg_files/figure-gfm/unnamed-chunk-10-2.svg)<!-- -->![](polyg_files/figure-gfm/unnamed-chunk-10-3.svg)<!-- -->

``` r
quality_p_final <- quality_p +
  coord_cartesian(clip = 'off') +
  annotate("text", 48.5, -2.2, label="sliding-window\ntrimming", size = 3.8) +
  annotate("segment", 48.5, -1.3, xend=48.5, yend=-0.6, arrow=arrow(length = unit(0.1, "npc")), size=1) +
  annotate("text", 101.5, -2.2, label="poly-G\ntrimming", size = 3.8) +
  annotate("segment", 101.5, -1.3, xend=101.5, yend=-0.6, arrow=arrow(length = unit(0.1, "npc")), size=1) #+
  #annotate("text", 76, -2.5, label='bold("base quality score")', size = 3.8, parse=TRUE) +
  #theme(legend.title=element_blank())
quality_p_final
```

![](polyg_files/figure-gfm/unnamed-chunk-11-1.svg)<!-- -->

###### Alternative visualization

``` r
top2 <- tibble(base=sequence_vector, quality=quality_vector) %>%
  {mutate(., position=seq(nrow(.)))} %>%
  ggplot(aes(x=position, y=quality)) +
  geom_line(color="grey") +
  geom_point(aes(fill=base), color="white", shape=22, size=2) +
  scale_fill_manual(values = c("#749dae", "#5445b1", "orange", "#cd3341")) +
  geom_segment(x = 48.5, xend = 48.5, y=0, yend=46, size=0.3) +
  geom_segment(x = 101.5, xend = 101.5, y=0, yend=46, size=0.3) +
  annotate("text", 48.5, 55, label="sliding-window\ntrimming", lineheight = 0.8, size=4.5) +
  annotate("text", 101.5, 55, label="poly-G\ntrimming", lineheight = 0.8, size=4.5) +
  scale_x_continuous(limits=c(0, 151), expand = c(0, 2)) +
  scale_y_continuous(limits=c(0, NA), expand = c(0, 3)) +
  coord_cartesian(clip = 'off') +
  labs(x="read position (in bp)", y="base quality") +
  theme_cowplot() +
  theme(legend.title = element_text(size=15),
        legend.text = element_text(size=12))
top2
```

![](polyg_files/figure-gfm/unnamed-chunk-12-1.svg)<!-- -->

#### Combined plot

``` r
top <- plot_grid(base_p_final, quality_p_final, rows = 2, rel_heights=c(1, 1.2)) + theme(plot.margin=unit(c(1, 1.5, 1, 1), unit="line"))
plot_grid(top, seq_content_p, labels = c('A', 'B'), label_size = 15, nrow = 2, rel_heights = c(2.5, 5))
```

![](polyg_files/figure-gfm/unnamed-chunk-13-1.svg)<!-- -->

``` r
figure_2 <- plot_grid(NULL, top2, NULL, seq_content_p, labels = c(NA, 'A', NA, 'B'), label_size = 15, ncol = 1, rel_heights = c(0.1, 2, 0.4, 4.5), label_y=c(NA, 1.08, NA, 1.02))
figure_2
```

![](polyg_files/figure-gfm/unnamed-chunk-14-1.svg)<!-- -->

``` r
ggsave("../figures/figure_2.pdf", figure_2, width=10, height=7, unit="in")
```

#### Persistence of poly-G tails in bam files

``` bash
fastqc /workdir/cod/greenland-cod/bam/QQL2011_884_55191_7_pe_bt2_gadMor3_minq20_sorted_dedup_overlapclipped_realigned.bam
fastqc /workdir/cod/greenland-cod/bam/UUM2010_068_55111_7_pe_bt2_gadMor3_minq20_sorted_dedup_overlapclipped_realigned.bam
fastqc /workdir/cod/greenland-cod/bam/IKE2011_976_55124_7_pe_bt2_gadMor3_minq20_sorted_dedup_overlapclipped_realigned.bam
mv /workdir/cod/greenland-cod/bam/*fastqc* /workdir/cod/greenland-cod/fastqc/
```

``` r
random_sample <- filter(sample_table_unmerged_pe, sample_id %in% c("QQL2011_884", "UUM2010_068", "IKE2011_976"))
path <- "../../cod/greenland-cod/fastqc/"
suffix <- "_pe_bt2_gadMor3_minq20_sorted_dedup_overlapclipped_realigned_fastqc"
type <- "poly-G trimming, mapped"
for (i in 1:3){
  sample_id <- random_sample$sample_id[i]
  sample_prefix <- random_sample$prefix_new[i]
  file_name <- str_c(path, sample_prefix, suffix)
  unzip(str_c(file_name, ".zip"), exdir = path, overwrite = FALSE)
  fastqc_data <- read_lines(file = str_c(file_name, "/fastqc_data.txt"))
  first_line <- which(str_detect(fastqc_data, ">>Per base sequence content")) + 1
  last_line <- which(str_detect(fastqc_data, ">>Per sequence GC content")) - 2
  per_base_seq_content_polyg_trimmed <- fastqc_data[first_line:last_line] %>%
    read_tsv() %>%
    rename(position=`#Base`) %>%
    pivot_longer(2:5, names_to = "base", values_to = "percentage") %>%
    mutate(sample_id=sample_id, type = type)
  if (i == 1) {
    per_base_seq_content_polyg_trimmed_final <- per_base_seq_content_polyg_trimmed
  } else {
    per_base_seq_content_polyg_trimmed_final <- bind_rows(per_base_seq_content_polyg_trimmed_final, per_base_seq_content_polyg_trimmed)
  }
}

seq_content_p <- per_base_seq_content_polyg_trimmed_final %>%
  mutate(position = as_factor(position)) %>%
  ggplot(aes(x=position, y=percentage, color=base, group=base)) +
  geom_line(size=0.8) +
  scale_color_manual(values = c("#749dae", "#5445b1", "orange", "#cd3341")) +
  facet_grid(str_c("ind ", as.numeric(as.factor(sample_id)))~type) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
seq_content_p
```

![](polyg_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Pipeline after sliding-window trimming

Sliding-window trimming is proved to be effective at removing polyG
tails. Therefore, we rerun everything with sliding-window trimmed fastq
files (the same adatper clipped fastq files are used for the HiSeq batch
where polyG is not an issue).

#### Map to reference genome

This was only run with NextSeq data because the raw bam files still
exist for the HiSeq batch.

All subsequent steps are run with NextSeq and HiSeq combined, because
intermediate files were deleted for the HiSeq batch.

Note that the reference genome (`gadMor3.fasta`) was downloaded from the
NCBI (<https://www.ncbi.nlm.nih.gov/assembly/GCF_902167405.1/>), but we
altered its chromosome names to match those of the `gadMor2` genome. See
[here](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/original_pipeline.md#first-rename-the-chromosomes-in-the-gadmor3-genome)
for we did the renaming.

``` bash
nohup bash /workdir/data-processing/scripts/low_coverage_mapping.sh \
/workdir/batch-effect/sample_lists/fastq_list_pe.txt \
/workdir/batch-effect/sample_lists/sample_table_unmerged.tsv \
/workdir/batch-effect/cut_right/ \
/workdir/batch-effect/ \
_adapter_clipped_qual_filtered_f_paired.fastq.gz \
_adapter_clipped_qual_filtered_r_paired.fastq.gz \
very-sensitive \
/workdir/cod/reference_seqs/gadMor3.fasta \
gadMor3 \
> /workdir/batch-effect/nohups/low_coverage_mapping_pe.nohup &
```

#### Sort the raw bam files

Note that a minimum mapping quality filter is no longer used in this
step.

Save the following script as
`/workdir/batch-effect/scripts/sort_raw_bam.sh`.

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

#### Deduplicate

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

#### Indel realignment

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

#### SNP calling

Average depth per site across all individuals is \~92x before polyG
trimming

``` bash
cd /workdir/batch-effect/
nohup bash /workdir/genomic-data-analysis/scripts/angsd_global_snp_calling.sh \
/workdir/batch-effect/sample_lists/bam_list_realigned.txt \
/workdir/batch-effect/ \
/workdir/cod/reference_seqs/gadMor3.fasta \
46 184 20 20 0.05 20 \
> /workdir/batch-effect/nohups/global_snp_calling_bam_list_realigned.nohup &
```

#### SAF, MAF, Fst estimation, and read depth count in each batch of data

``` bash
## MAF and SAF (minInd=20)
nohup bash /workdir/genomic-data-analysis/scripts/get_maf_per_pop.sh \
/workdir/batch-effect/ \
/workdir/batch-effect/sample_lists/sample_table_merged.tsv \
6 \
bam_list_realigned_ \
/workdir/cod/reference_seqs/gadMor3.fasta \
/workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20.txt \
20 184 20 20 \
> /workdir/batch-effect/nohups/get_maf_per_pop.nohup &
## Fst
nohup bash /workdir/genomic-data-analysis/scripts/get_fst.sh \
/workdir/batch-effect/angsd/popminind20/ \
/workdir/batch-effect/sample_lists/sample_table_merged.tsv \
6 \
_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_popminind20 \
> /workdir/batch-effect/nohups/get_fst.nohup &
## Get depth count from se samples without a mapping quality filter (minInd=2)
cd /workdir/batch-effect
nohup /workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_per_pop/bam_list_realigned_se.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/popminind2/bam_list_realigned_se_anymapq \
-doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 16 -setMinDepth 2 -minInd 2 -minQ 20 \
-sites /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20.txt \
-rf /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20.chrs \
>& nohups/get_depth_anymapq_bam_list_realigned_se.log &
```

#### LD pruning

###### Downsample the mafs file to \~1,000,000 SNPs

``` r
library(tidyverse)
## Read in the mafs file
mafs <- read_tsv("../angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20.mafs.gz")
## Total number of SNPs
nrow(mafs)
```

    ## [1] 5204764

``` r
## In a downsampled SNP list, keep 1 SNP in every 5 SNPs (this yields 1028474 SNPs)
mafs_filtered <- mutate(mafs, 
                        keep = knownEM >= 0.05,
                        row_number = ifelse(keep, cumsum(keep), -1),
                        keep = ifelse(row_number%%5!=0, F, keep),
                        keep = ifelse(!str_detect(chromo, "LG"), F, keep)) %>%
  dplyr::select(-row_number)
filter(mafs_filtered, keep==T) %>% nrow()
```

    ## [1] 1028474

``` r
## Filter out SNPs in inversions
lg = c("LG01", "LG02", "LG07", "LG12")
min_pos = c(11495229, 641666, 16830817, 546885)
max_pos = c(28834444, 4495925, 26562127, 13864687)
for (i in 1:4){
  if (i==1){
      mafs_inversion_filtered <- mutate(mafs_filtered, keep = ifelse(chromo == lg[i] & position > min_pos[i] & position < max_pos[i], FALSE, keep))
  } else {
      mafs_inversion_filtered <- mutate(mafs_inversion_filtered, keep = ifelse(chromo == lg[i] & position > min_pos[i] & position < max_pos[i], FALSE, keep))
  }
}
mafs_inversion_filtered %>%
  filter(keep==T) %>%
  ggplot(aes(x=position, y=chromo, fill=chromo)) +
  ggridges::geom_density_ridges(alpha=0.5) +
  scale_fill_viridis_d() +
  theme_cowplot() +
  theme(legend.position = "none")
```

![](polyg_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
mafs_inversion_filtered %>%
  filter(keep==T) %>%
  nrow()
```

    ## [1] 944554

###### Write the filtered SNP list

``` r
mafs_inversion_filtered %>%
  filter(keep==T) %>%
  dplyr::select(1:4) %>%
  write_tsv("../angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.txt.gz", col_names = F)
(which(mafs_inversion_filtered$keep)+1) %>% # +1 is necessary because the beagle file has a title line
  as.integer() %>%
  write_lines("../angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.pos.idx")
```

###### Downsample the beagle file

``` bash
zcat /workdir/batch-effect/angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20.beagle.gz | \
awk 'NR==FNR{ pos[$1]; next }FNR in pos' \
/workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.pos.idx - | \
gzip \
> /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.beagle.gz
```

###### Estimate LD

``` bash
nohup /workdir/programs/ngsLD/ngsLD \
--geno /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.beagle.gz \
--pos /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.txt.gz \
--n_ind 163 \
--n_sites 944554 \
--out /workdir/batch-effect/angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.ld \
--probs \
--rnd_sample 1 \
--seed 42 \
--max_kb_dist 10 \
--n_threads 30 \
> /workdir/batch-effect/nohups/run_ngsLD_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.nohup &
```

###### Remove some columns in the LD file

``` r
ld <- read_tsv("../angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.ld", col_names = F)
ld %>%
  dplyr::select(-X2, -X3, -X5, -X6) %>%
  write_tsv("../angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_reformatted.ld", col_names = F)
```

###### Run the LD pruning script

``` bash
nohup perl /workdir/programs/ngsLD/scripts/prune_graph.pl \
--in_file /workdir/batch-effect/angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_reformatted.ld \
--max_kb_dist 10 \
--min_weight 0.5 \
--out /workdir/batch-effect/angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.ld \
> /workdir/batch-effect/nohups/ld_prune_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.nohup &
```

###### Generate LD pruned SNP list

``` r
downsampled_snp_list <- read_tsv("../angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled.txt.gz", col_names = c("lg", "position", "major", "minor"))
ld_pruned_snp_list <- read_delim("../angsd/bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.ld", delim = ":", col_names = c("lg", "position")) %>%
  arrange(lg, position) %>%
  mutate(keep=TRUE) %>%
  left_join(downsampled_snp_list, ., by = c("lg", "position")) %>%
  mutate(keep=ifelse(is.na(keep), FALSE, keep))
identical(downsampled_snp_list$lg, ld_pruned_snp_list$lg)
identical(downsampled_snp_list$position, ld_pruned_snp_list$position)
## Generate a SNP list
ld_pruned_snp_list %>%
  filter(keep == TRUE) %>%
  dplyr::select(lg, position, major, minor) %>%
  write_tsv("../angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.txt", col_names = F)
## Generate a chromosome list
ld_pruned_snp_list %>%
  filter(keep == TRUE) %>%
  .$lg %>%
  unique() %>%
  write_lines("../angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.chrs")
```

``` bash
## Index the snp list
/workdir/programs/angsd0.931/angsd/angsd sites index /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.txt
```

#### Get the covariance matrix with ANGSD using LD pruned SNP list

``` bash
cd /workdir/batch-effect/

nohup /workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_realigned.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/bam_list_realigned_downsampled_unlinked \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 16 -setMinDepth 46 -setMaxDepth 184 -minInd 20 -minQ 20 -minMapQ 20 -minMaf 0.05 \
-doIBS 2 -makematrix 1 -doCov 1 \
-sites /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.txt \
-rf /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.chrs \
>& nohups/run_pca_bam_list_realigned_downsampled_unlinked.log &
```
