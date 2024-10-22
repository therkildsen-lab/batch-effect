A tutorial for the detection and mitigation of batch effects with
low-coverage whole genome sequencing data
================

  - [Introduction](#introduction)
  - [Some preparation](#some-preparation)
  - [Poly-G tail](#poly-g-tail)
      - [Poly-G tail trimming with
        fastp](#poly-g-tail-trimming-with-fastp)
      - [Sliding window trimming with
        fastp](#sliding-window-trimming-with-fastp)
      - [Run FastQC](#run-fastqc)
      - [Visualize base composition](#visualize-base-composition)
      - [Heterozygosity estimation](#heterozygosity-estimation)
      - [Visualization](#visualization)
  - [Base quality score
    miscalibration](#base-quality-score-miscalibration)
      - [Heterozygosity estimation](#heterozygosity-estimation-1)
      - [Visualization](#visualization-1)
  - [Reference bias](#reference-bias)
      - [SNP calling and PCA](#snp-calling-and-pca)
      - [Get saf file per population](#get-saf-file-per-population)
      - [Get 2dSFS and Fst](#get-2dsfs-and-fst)
      - [Get depth in HiSeq-125SE batch without mapping quality
        filter](#get-depth-in-hiseq-125se-batch-without-mapping-quality-filter)
      - [Examine Fst results](#examine-fst-results)
      - [Proportion of reads with lower mapping quality scores as a
        proxy for reference
        bias](#proportion-of-reads-with-lower-mapping-quality-scores-as-a-proxy-for-reference-bias)
      - [Fst before and after
        correction](#fst-before-and-after-correction)
      - [Examine PCA results](#examine-pca-results)
      - [Corrected PCA](#corrected-pca)
  - [DNA degradation](#dna-degradation)
      - [Heterozygosity estimation](#heterozygosity-estimation-2)
      - [Visualization](#visualization-2)
  - [Conclusion](#conclusion)

## Introduction

This is a tutorial accompanying the paper [Batch effects in population
genomic studies with low-coverage whole genome sequencing data: causes,
detection, and
mitigation](https://doi.org/10.22541/au.162791857.78788821/v2) by R.
Nicolas Lou and Nina Overgaard Therkildsen. In this tutorial, we provide
a small subset of the data used in our paper and offer hands-on
exercises for the detection and mitigation of batch effects in
low-coverage whole genome sequencing (lcWGS) data.

We will focus on four main causes of batch effects in this tutorial:

1.  Difference in sequencing chemistry (two vs. four channel) leading to
    the presence / absence of poly-G tails

2.  Separate sequencing runs leading to different levels of
    miscalibration in base quality scores

3.  Difference in read type and read length leading to different levels
    of reference bias / alignment error

4.  Difference in levels of DNA degradation

After taking this tutorial, you will learn about the different causes of
batch effects, and will be able to perform some simple bioinformatic
methods to detect and mitigate their impacts.

This tutorial is under active development and will be finalized before
the publication of our paper.

## Some preparation

To start, clone this GitHub repo to your Linux server (e.g. `git clone
https://github.com/therkildsen-lab/batch-effect.git`). Even though the
sizes of our individual files are not very large, this may still take
some time.

Change your working directory to the tutorial folder within this GitHub
repo (e.g. `cd /workdir/batch-effect/tutorial`) and inspect the content
of this folder (e.g. `ls`). As you can see, there is a folder called
`data` that stores all the data needed for this tutorial, and some
markdown files.

Then, define the following variables on your Linux server. Make sure
that you **modify these paths** to reflect the locations in your file
system. If some of the programs listed below are not installed on your
server, install them first.

``` bash
## Location of the tutorial folder in the batch-effect Github repo
BASEDIR=/workdir/batch-effect/tutorial
## Path to samtools (https://github.com/samtools/)
SAMTOOLS=samtools
## Path to fastp (https://github.com/OpenGene/fastp)
FASTP=/workdir/programs/fastp_0.19.7/fastp
## Path to FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
FASTQC=fastqc
## Path to AGNSD (https://github.com/ANGSD/angsd)
ANGSD=/programs/angsd0.930/angsd/angsd 
## Path to the realSFS module in AGNSD
REALSFS=/workdir/programs/angsd0.931/angsd/misc/realSFS
## Path to a subsetted gadMor3 reference genome
REFERENCE=$BASEDIR/data/gadMor3_subsetted.fasta
```

Index the subsetted reference genome.

``` bash
$SAMTOOLS faidx $REFERENCE
```

You will need to create a `results` folder within the `tutorial` folder
to store your output files. To do so, run

``` bash
makdir $BASEDIR/results
```

Also, you will need to use R in this tutorial for data wrangling and
visualization. We recommend you use RStudio server, but you can also
launch another session on your Linux server and use R there. Install the
following packages with `install.packages()` if you don’t have them
installed already, and load them using `library()` as shown below.

``` r
library(tidyverse)
library(cowplot)
library(ggstatsplot)
library(statsExpressions)
```

## Poly-G tail

Across Illumina platforms, an important change in recent years is the
shift from a four-channel system (used e.g., on HiSeq instruments) where
each DNA base is detected with a different fluorescent dye, to a
two-channel chemistry, that uses the combinations of two different dyes.
With the two-channel system (implemented on newer platforms like NextSeq
and NovaSeq), G is called when there is little to no fluorescence
signal. Accordingly, the absence of a signal can result from a true G
base in the DNA template, but any low-intensity fluorescence signal
(regardless of the true base) may also lead to a G call, which becomes
problematic.

Here, we will first inspect a pair of fastq files
(`before_trimming_f.fastq.gz` and `before_trimming_r.fastq.gz`)
generated by a NextSeq sequencer. Specifically, we will perform two
different trimming operations, and will compare the base composition per
position in the raw and trimmed fastq files.

Then, we will estimate heterozygosity from bam files to assess the
impact that incomplete poly-G tail removal has on downstream analyses.

#### Poly-G tail trimming with fastp

Poly-G tail trimming looks for long stretches of G base at the end of
sequencing reads and removes them from the fastq file.

``` bash
$FASTP --trim_poly_g -L -A \
  -i $BASEDIR/data/before_trimming_f.fastq.gz \
  -I $BASEDIR/data/before_trimming_r.fastq.gz \
  -o $BASEDIR/results/poly_g_trimming_f.fastq.gz \
  -O $BASEDIR/results/poly_g_trimming_r.fastq.gz \
  -h $BASEDIR/results/poly_g_trimming.html \
  -j $BASEDIR/results/poly_g_trimming.json
```

#### Sliding window trimming with fastp

Sliding window trimming moves a window from the beginning to the end of
a read. Once the average base quality score falls below a certain
threshold, everything after the window will be removed.

``` bash
$FASTP --trim_poly_g -L -A --cut_right \
  -i $BASEDIR/data/before_trimming_f.fastq.gz \
  -I $BASEDIR/data/before_trimming_r.fastq.gz \
  -o $BASEDIR/results/sliding_window_trimming_f.fastq.gz \
  -O $BASEDIR/results/sliding_window_trimming_r.fastq.gz \
  -h $BASEDIR/results/sliding_window_trimming.html \
  -j $BASEDIR/results/sliding_window_trimming.json
```

#### Run FastQC

We use fastQC to compute the average base composition of each position
along the sequencing reads. We will only use the forward reads here.

``` bash
$FASTQC $BASEDIR/data/before_trimming_f.fastq.gz -o $BASEDIR/results/
$FASTQC $BASEDIR/results/poly_g_trimming_f.fastq.gz -o $BASEDIR/results/
$FASTQC $BASEDIR/results/sliding_window_trimming_f.fastq.gz -o $BASEDIR/results/
```

#### Visualize base composition

Run the following with R to compile and plot the base composition
generation by FastQC.

``` r
## Define the paths to the results/ folder as your basedir
basedir="/workdir/batch-effect/tutorial/results/"
## Loop through the FastQC results from the raw fastq file, poly-G trimmed fastq file, and sliding-window trimmed fastq file
file_list=c("before_trimming_f_fastqc", "poly_g_trimming_f_fastqc", "sliding_window_trimming_f_fastqc")
type_list <- c("no trimming", "poly-G trimming", "sliding-window trimming")
for (i in 1:3){
  file <- file_list[i]
  type <- type_list[i]
  unzip(str_c(basedir, file, ".zip"), exdir = basedir)
  fastqc_data <- read_lines(file = str_c(basedir, file, "/fastqc_data.txt"))
  first_line <- which(str_detect(fastqc_data, ">>Per base sequence content")) + 1
  last_line <- which(str_detect(fastqc_data, ">>Per sequence GC content")) - 2
  per_base_seq_content_polyg_trimmed <- fastqc_data[first_line:last_line] %>%
    read_tsv() %>%
    rename(position=`#Base`) %>%
    pivot_longer(2:5, names_to = "base", values_to = "percentage") %>%
    mutate(type = type)
  if (i == 1) {
    per_base_seq_content_polyg_trimmed_final <- per_base_seq_content_polyg_trimmed
  } else {
    per_base_seq_content_polyg_trimmed_final <- bind_rows(per_base_seq_content_polyg_trimmed_final, per_base_seq_content_polyg_trimmed)
  }
}
## Plot
seq_content_p <- per_base_seq_content_polyg_trimmed_final %>%
  mutate(position = as_factor(position)) %>%
  ggplot(aes(x=position, y=percentage, color=base, group=base)) +
  geom_line(size=0.8) +
  scale_color_manual(values = c("#749dae", "#5445b1", "orange", "#cd3341")) +
  xlab("read position (in bp)") +
  facet_wrap(~type, nrow = 3) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
seq_content_p
```

![](tutorial_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

The raw fastq file has a strong enrichment of G bases toward the end of
sequencing reads.

#### Heterozygosity estimation

The differentially trimmed fastq files need to be aligned to a reference
genome and go through several quality control steps. For the sake of
time, we will skip these steps in this tutorial. If you are interested
in them, please see our full
[pipeline](https://github.com/therkildsen-lab/batch-effect) for our
batch-effect paper as well as a
[tutorial](https://github.com/nt246/lcwgs-guide-tutorial) for the
processing and analysis of lcWGS data.

Here, the quality-controlled alignment files (in bam format) are
provided to you directly. They are from “pop 5” in our paper, and only
contain a 10Mb section of a chromosome. Samples that were sequenced on a
NextSeq platform were either trimmed using the poly-G trimming option in
fastp, or through sliding window trimming. Samples from the HiSeq
platform were not trimmed. The names of these bam files are stored in
the `bam_list_poly_g.txt` file under `data/`.

First, we loop through `bam_list_poly_g.txt` to estimate heterozygosity
from these bam files using ANGSD.

``` bash
MINDP=2
MAXDP=10
MINQ=20
MINMAPQ=30

for LINE in `cat $BASEDIR/data/bam_list_poly_g.txt`; do
    NAME_TEMP=`echo "${LINE%.*}"`
    NAME=`echo "${NAME_TEMP##*/}"`
      echo $NAME
    OUTBASE=$NAME'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ
    
    ## Get saf file
    $ANGSD \
    -i $LINE \
    -anc $REFERENCE \
    -out $BASEDIR/results/$OUTBASE \
    -doSaf 1 \
    -GL 1 \
    -P 8 \
    -doCounts 1 \
    -setMinDepth $MINDP \
    -setMaxDepth $MAXDP \
    -minQ $MINQ \
    -minmapq $MINMAPQ 

    ## Estimate sfs
    /workdir/programs/angsd0.931/angsd/misc/realSFS \
    $BASEDIR/results/$OUTBASE'.saf.idx' \
    -tole 0.0000001 \
    -P 8 \
    -seed 42 \
    > $BASEDIR/results/$OUTBASE'.ml'
done
```

#### Visualization

In R, we read in the estimated heterozygosity and compare the estimated
heterozygosity in the same samples that have gone through different
types of quality trimming.

Note that each NextSeq sample in the sample table has a bam file that
has gone through poly-G trimming, and a bam file that has gone through
sliding-window trimming. HiSeq samples were not trimmed, so the two bam
files for each HiSeq sample are essentially the same.

``` r
## Define the tutorial/ folder as your basedir
basedir="/workdir/batch-effect/tutorial/"
## Read in a table of sample information
sample_table <- read_tsv(str_c(basedir, "data/sample_table_poly_g.tsv"))
## Read in the estimated heterozygosity for each sample in the sample table
for (i in seq_len(nrow(sample_table))){
  data_type <- sample_table$data_type[i]
  if (data_type=="se"){
    prefix <- str_c(sample_table$sample_seq_id[i], "_bt2_gadMor3_sorted_dedup_realigned_subsetted_")
  } else {
    prefix <- str_c(sample_table$sample_seq_id[i], "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_subsetted_")
  }
  for (trimming in c("poly_g_trimming", "sliding_window_trimming")) {
    ml <- read_delim(str_c(basedir, "results/", prefix, trimming, "_mindp2_maxdp10_minq20_minmapq30.ml"), col_names = FALSE, delim = " ")
    het <- ml$X2/(ml$X1 + ml$X2 + ml$X3)
    line <- tibble(het=het, sample_id=sample_table$sample_id_corrected[i], data_type=sample_table$data_type[i], trimming=trimming)
    if (i==1 & trimming=="poly_g_trimming") {
      het_final <- line
    } else{
      het_final <- bind_rows(het_final, line)
    }
  }
}
## Plot
set.seed(42)
het_final %>%
  mutate(batch=ifelse(data_type=="pe", "NextSeq-150PE", "HiSeq-125SE")) %>%
  mutate(trimming=ifelse(trimming=="poly_g_trimming", "Poly-G trimming", "Sliding window trimming")) %>%
  ggplot(aes(x="", y=het)) +
  geom_boxplot(outlier.alpha = 0, color="black", size=0.2, width=0.2) +
  geom_jitter(aes(color=batch), height = 0, width = 0.1, size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  ylab(expression(paste("heterozygosity (in ", 10^-3, ")"))) +
  xlab(" ") +
  coord_flip() +
  facet_grid(trimming~"pop 5") +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
```

![](tutorial_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

As shown above, remnant poly-G tails severely inflate heterozygosity
estimates in NextSeq samples if only poly-G trimming is performed. In
contrast, sliding window trimming can substantially correct for such
batch effects. This suggests that poly-G tails are a major cause of
batch effects in heterozygosity estimates in our data.

## Base quality score miscalibration

In an ideal scenario, a base quality score should accurately reflect the
probability of the base call being correct. In practice, however, these
scores are often incorrectly calibrated (Callahan et al., 2016; Ni &
Stoneking, 2016), which can lead to batch effects if the levels of such
biases differ across sequencing runs.

Here, we will compare the heterozygosity estimation using a relaxed
vs. stringent base quality filter to detect base quality score
miscalibration in samples from the same Atlantic cod population that are
split into two sequencing batches. We will start from sliding-window
trimmed bam files in “pop 5”; these bam files contain an entire
chromosome (LG06).

#### Heterozygosity estimation

We use ANGSD to estimate the heterozygosity of individual bam files in
`bam_list_bq.txt` using two different minimum base quality filters (20
vs. 33).

``` bash
MINDP=2
MAXDP=10
MINMAPQ=30

for MINQ in {20,33}; do
  for LINE in `cat $BASEDIR/data/bam_list_bq.txt`; do
      NAME_TEMP=`echo "${LINE%.*}"`
      NAME=`echo "${NAME_TEMP##*/}"`
      echo $NAME
      OUTBASE=$NAME'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ
      
      ## Get saf file
      $ANGSD \
      -i $LINE \
      -anc $REFERENCE \
      -out $BASEDIR/results/$OUTBASE \
      -doSaf 1 \
      -GL 1 \
      -P 8 \
      -doCounts 1 \
      -setMinDepth $MINDP \
      -setMaxDepth $MAXDP \
      -minQ $MINQ \
      -minmapq $MINMAPQ 

      ## Estimate sfs
      /workdir/programs/angsd0.931/angsd/misc/realSFS \
      $BASEDIR/results/$OUTBASE'.saf.idx' \
      -tole 0.0000001 \
      -P 8 \
      -seed 42 \
      > $BASEDIR/results/$OUTBASE'.ml'
  done
done
```

#### Visualization

We perform visualization and statistical analysis using R.

``` r
## Define the path to the tutorial/ folder as basedir
basedir="/workdir/batch-effect/tutorial/"
## Sample table for this set of samples
sample_table <- read_tsv(str_c(basedir, "data/sample_table_bq.tsv"))
for (i in seq_len(nrow(sample_table))){
  data_type <- sample_table$data_type[i]
  if (data_type=="se"){
    prefix <- str_c(sample_table$sample_seq_id[i], "_bt2_gadMor3_sorted_dedup_realigned_subsetted_bq_mindp2_maxdp10_minq")
  } else {
    prefix <- str_c(sample_table$sample_seq_id[i], "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_subsetted_bq_mindp2_maxdp10_minq")
  }
  for (minq in c(20, 33)) {
    ml <- read_delim(str_c(basedir, "results/", prefix, minq, "_minmapq30.ml"), col_names = FALSE, delim = " ")
    het <- ml$X2/(ml$X1 + ml$X2 + ml$X3)
    line <- tibble(het=het, sample_id=sample_table$sample_id_corrected[i], data_type=sample_table$data_type[i], minq=minq)
    if (i==1 & minq==20) {
      het_final <- line
    } else{
      het_final <- bind_rows(het_final, line)
    }
  }
}
## Plot
set.seed(42)
het_final %>%
  mutate(data_type=ifelse(data_type=="pe", "NextSeq-150PE (more miscalibration in base quality scores)", "HiSeq-125SE (less miscalibration in base quality scores)")) %>%
  mutate(filter=ifelse(minq==20, "relaxed base quality filter (minQ = 20)", "stringent base quality filter (minQ = 33)")) %>%
  mutate(data_type=factor(data_type, levels = c("NextSeq-150PE (more miscalibration in base quality scores)", "HiSeq-125SE (less miscalibration in base quality scores)"))) %>%
  ggstatsplot::grouped_ggwithinstats(
  x = filter,
  y = het,
  #type = "np", # non-parametric statistics
  point.path.args=list(alpha=0.2),
  xlab = element_blank(),
  ylab = "estimated heterozygosity",
  grouping.var = data_type,
  ggtheme = theme_ggstatsplot(),
  bf.message = FALSE,
  ggplot.component = list(scale_color_viridis_d(begin = 0.4, end = 0.7, option = "A"),
                          theme(panel.grid = element_blank(),
                                axis.line = element_line()))
  )
```

![](tutorial_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

As shown above, NextSeq-150PE samples get significantly lower
heterozygosity estimates after the filter is applied, and the
HiSeq-125PE samples get significantly higher heterozygosity estimates
with the filter.

This suggests that base qualities are likely to be overestimated in the
NextSeq-150PE batch, and underestimated in the HiSeq-125PE batch. This
could cause inflation of heterozygosity in the NextSeq-150PE batch
(i.e. sequencing errors are interpreted as true mutations). Therefore,
using a more stringent base quality threshold for both batches can be an
effective mitigation strategy.

## Reference bias

Compared to longer paired-end reads, shorter single-end reads carrying
bases that are different from the reference are less likely to be
aligned to the reference genome (either correctly or incorrectly) with
high confidence (e.g., due to the presence of an analogous sequence on
the genome), and therefore tend to receive low mapping quality scores.
As a result, when a stringent mapping quality threshold is imposed,
shorter single-end reads are more likely to be affected.

Here, we perform PCA and Fst estimation on a 0.5Mb section of a
chromosome using samples from four Atlantic cod populations (pop 2, 3,
4, 6). Some of the samples were sequenced with paired-end 150bp reads,
and others were sequenced with single-end 125bp reads.

#### SNP calling and PCA

First, we perform SNP calling and PCA using ANGSD.

``` bash
$ANGSD -b $BASEDIR/data/bam_list_rb.txt \
-anc $REFERENCE \
-out $BASEDIR/results/bam_list_rb \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 2 -makematrix 1 -doCov 1 -P 8 \
-SNP_pval 1e-6 -setMinDepth 10 -setMaxDepth 92 -minInd 2 -minQ 20 -minMaf 0.05 -minMapQ 20 

zcat $BASEDIR/results/bam_list_rb.mafs.gz | cut -f 1,2,3,4 | tail -n +2 > $BASEDIR/results/bam_list_rb.snp_list.txt
$ANGSD sites index $BASEDIR/results/bam_list_rb.snp_list.txt
```

#### Get saf file per population

We get the sample allele frequency likelihoods (saf files) using ANGSD
in each of the two batches (i.e., NextSeq-150PE, HiSeq-125SE).

``` bash
## NextSeq-150PE
$ANGSD \
-b $BASEDIR/data/bam_list_rb_pe.txt \
-anc $REFERENCE \
-out $BASEDIR/results/bam_list_rb_pe \
-doSaf 1 \
-GL 1 \
-P 8 \
-doCounts 1 \
-dumpCounts 1 \
-setMinDepth 5 \
-minQ 20 \
-minmapq 20 \
-sites $BASEDIR/results/bam_list_rb.snp_list.txt

## HiSeq-125SE
$ANGSD \
-b $BASEDIR/data/bam_list_rb_se.txt \
-anc $REFERENCE \
-out $BASEDIR/results/bam_list_rb_se \
-doSaf 1 \
-GL 1 \
-P 8 \
-doCounts 1 \
-dumpCounts 1 \
-setMinDepth 5 \
-minQ 20 \
-minmapq 20 \
-sites $BASEDIR/results/bam_list_rb.snp_list.txt
```

#### Get 2dSFS and Fst

We compute Fst between the two batches.

``` bash
$REALSFS \
$BASEDIR/results/bam_list_rb_pe.saf.idx \
$BASEDIR/results/bam_list_rb_se.saf.idx \
-P 8 \
-seed 42 \
-maxiter 500 \
> $BASEDIR/results/bam_list_rb_pe_bam_list_rb_se.2dSFS

$REALSFS fst index \
$BASEDIR/results/bam_list_rb_pe.saf.idx \
$BASEDIR/results/bam_list_rb_se.saf.idx \
-sfs $BASEDIR/results/bam_list_rb_pe_bam_list_rb_se.2dSFS \
-fstout $BASEDIR/results/bam_list_rb_pe_bam_list_rb_se.alpha_beta

$REALSFS fst print \
$BASEDIR/results/bam_list_rb_pe_bam_list_rb_se.alpha_beta.fst.idx \
> $BASEDIR/results/bam_list_rb_pe_bam_list_rb_se.alpha_beta.txt
```

#### Get depth in HiSeq-125SE batch without mapping quality filter

We count the sequencing depth at each SNP in the HiSeq-125SE batch
without a mapping quality filter.

``` bash
$ANGSD \
-b $BASEDIR/data/bam_list_rb_se.txt \
-anc $REFERENCE \
-out $BASEDIR/results/bam_list_rb_se_minmapq0 \
-P 8 \
-doCounts 1 \
-dumpCounts 1 \
-setMinDepth 2 \
-minQ 20 \
-minmapq 0 \
-sites $BASEDIR/results/bam_list_rb.snp_list.txt
```

#### Examine Fst results

Plot Fst in R.

``` r
basedir="/workdir/batch-effect/tutorial/"
fst <- read_tsv(str_c(basedir, "results/bam_list_rb_pe_bam_list_rb_se.alpha_beta.txt"), col_names = FALSE) %>%
  transmute(pos=X2, fst=X3/X4)
fst %>%
  ggplot(aes(x=pos, y=fst)) +
  geom_point(size=0.1) +
  theme_cowplot()
```

![](tutorial_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

At some SNPs, Fst is very high between samples from the same populations
but two different batches, which is not expected. Let’s investigate what
is causing such high Fst values.

#### Proportion of reads with lower mapping quality scores as a proxy for reference bias

We calculate the proportion of reads with mapping quality scores lower
than 20 in the HiSeq-125SE batch for each SNP using R, and check whether
this relates to the Fst outliers.

``` r
## Depth in HiSeq-125SE without mapping quality filter
depth_minmapq0 <- read_tsv(str_c(basedir, "results/bam_list_rb_se_minmapq0.pos.gz")) %>%
  rename(minmapq0=totDepth)
## Depth in HiSeq-125SE with minimum mapping quality = 20
depth_minmapq20 <- read_tsv(str_c(basedir, "results/bam_list_rb_se.pos.gz"))%>%
  rename(minmapq20=totDepth)
## Calculate depth ratio for each SNP
depth_joined <- left_join(depth_minmapq20, depth_minmapq0) %>%
  mutate(depth_ratio=1-minmapq20/minmapq0) %>%
  dplyr::select(pos, depth_ratio)
## Join fst with depth ratio
fst_depth_ratio <- fst %>%
  left_join(depth_joined) %>%
  mutate(type=ifelse(fst>0.2, "Fst outliers", "all other SNPs"),
         type=fct_relevel(type, c("Fst outliers", "all other SNPs")))
## Fst vs. depth ratio
fst_depth_ratio %>%
  ggplot(aes(x=depth_ratio, y=fst)) +
  geom_point(size=0.2) +
  theme_cowplot()
```

![](tutorial_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
## Get test stats
fst_depth_ratio %>%
  ggbetweenstats(y=depth_ratio, x=type, output="subtitle", bf.message = FALSE)
```

    ## paste(italic("t")["Welch"], "(", "18.08", ") = ", "11.17", ", ", 
    ##     italic("p"), " = ", "1.51e-09", ", ", widehat(italic("g"))["Hedges"], 
    ##     " = ", "2.75", ", CI"["95%"], " [", "1.73", ", ", "3.76", 
    ##     "], ", italic("n")["obs"], " = ", "5,202")

``` r
fst_depth_ratio_stats <- fst_depth_ratio %>%
  centrality_description(type, depth_ratio)
## Distrution of depth ratio in Fst outliers (those with Fst > 0.2) vs. all other SNPs
fst_depth_ratio %>%
  ggplot(aes(x=depth_ratio)) +
  geom_density(mapping = aes(fill=type), alpha=0.3) +
  geom_vline(data=fst_depth_ratio_stats, aes(xintercept = depth_ratio)) +
  geom_label(data=fst_depth_ratio_stats, aes(label=expression), y=4, parse=TRUE) +
  scale_fill_viridis_d() +
  labs(x="proportion of reads with mapping quality lower than 20",
       y="density",
       subtitle=expression(paste(italic("t")["Welch"], "(", "18.08", ") = ", "-11.17", ", ", 
    italic("p"), " = ", "1.51e-09", ", ", widehat(italic("g"))["Hedges"], 
    " = ", "-2.75", ", CI"["95%"], " [", "-3.76", ", ", "-1.73", 
    "], ", italic("n")["obs"], " = ", "5,202"))) +
  theme_ggstatsplot() +
  theme(panel.grid = element_blank(),
        axis.line = element_line())
```

![](tutorial_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

Fst outliers are enriched in regions where a high proportion of
single-end 125bp reads get low mapping quality scores, a signal of
heightened level reference bias in one batch of our data. Therefore, we
can exclude such SNPs (e.g., those with \> 10% reads having mapping
quality scores lower than 20 in the single-end 125bp batch) from the
analysis.

#### Fst before and after correction

Compare Fst before and after correction in R.

``` r
fst_depth_ratio %>%
  filter(depth_ratio < 0.1) %>%
  mutate(type="after") %>%
  bind_rows(mutate(fst_depth_ratio, type="before")) %>%
  mutate(type=fct_relevel(type, c("before", "after"))) %>%
  ggplot(aes(x=pos, y=fst)) +
  geom_point(size=0.1) +
  facet_wrap(~type, ncol=1) +
  theme_cowplot()
```

![](tutorial_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

On this segment of the genome, our mitigation strategy is particularly
successful: almost all of the outlier SNPs are removed.

#### Examine PCA results

Perform eigendecomposition and data visualization in R.

``` r
## Define the tutorial/ folder as your basedir
basedir="/workdir/batch-effect/tutorial/"
## This is how the different populations are named in our paper
rename_pop <- tibble(population = c("ITV2011", "KNG2011", "QQL2011", "BUK2011", "IKE2011", "PAA2011", "ATP2011", "NAR2008", "UUM2010"),
                     population_new =c("pop 1", "pop 2", "pop 3", "pop 4", "pop 5", "pop 6", "pop 7", "pop 8", "pop 9"))
## Join new population names with the sample table
sample_table <- read_tsv(str_c(basedir, "data/sample_table_rb.tsv")) %>%
  dplyr::select(-population_new) %>%
  left_join(rename_pop) %>%
  mutate(batch=ifelse(data_type=="se", "HiSeq-125SE", "NextSeq-150PE"))
## Order of samples in the covariance matrix is the same as in the bam list, so we extract the order of samples from the bam list
bam_list <- read_tsv(str_c(basedir, "data/bam_list_rb.txt"), col_names = FALSE) %>%
  transmute(sample_id=str_sub(X1, 37, 47))
## Read in the covariance matrix
genome_cov <- read_tsv(str_c(basedir, "results/bam_list_rb.covMat"), col_names = F)[1:nrow(sample_table),1:nrow(sample_table)] %>% as.matrix
## Perform eigendecomposition on the covariance matrix, and join it with the sample table
pca_table <- eigen(genome_cov)$vectors %>% 
  as_tibble() %>%
  transmute(PC1=V1, PC2=V2) %>%
  bind_cols(bam_list, .) %>%
  left_join(sample_table, ., by=c("sample_id_corrected"="sample_id"))
## Plot PCA result
pca_table %>%
  left_join(rename_pop) %>%
  ggplot(aes(x=PC1, y=PC2, color=batch)) +
  facet_wrap(~population_new) +
  geom_point() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.5))
```

![](tutorial_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

Because this subset of our dataset is very small, individual level PCA
has limited power, and therefore there is not a signal of batch effects
in our PCA results. However, we still provide the code to mitigate batch
effects when they are observed in the PCA results.

#### Corrected PCA

First, let’s come up with a new SNP list excluding SNPs that are
affected by reference bias.

``` r
## Come up with a new SNP list
read_tsv(str_c(basedir, "results/bam_list_rb.snp_list.txt"), col_names = FALSE) %>%
  semi_join(filter(fst_depth_ratio, depth_ratio < 0.1), by=c("X2"="pos")) %>%
  write_tsv(str_c(basedir, "results/bam_list_rb.depth_ratio_filtered_snp_list.txt"), col_names = FALSE)
```

Now, let’s run PCA at this new SNP list.

``` bash
## Index the new SNP list with ANGSD
$ANGSD sites index $BASEDIR/results/bam_list_rb.depth_ratio_filtered_snp_list.txt

## Run PCA with ANGSD, restricting the analysis on our new SNP list
$ANGSD -b $BASEDIR/data/bam_list_rb.txt \
-anc $REFERENCE \
-out $BASEDIR/results/bam_list_rb_depth_ratio_filtered \
-GL 1 -doMajorMinor 1 -doCounts 1 -doIBS 2 -makematrix 1 -doCov 1 -P 8 \
-minQ 20 -minMapQ 20 \
-sites $BASEDIR/results/bam_list_rb.depth_ratio_filtered_snp_list.txt
```

Use R to perform eigendecomposition and visualization.

``` r
basedir="/workdir/batch-effect/tutorial/"
rename_pop <- tibble(population = c("ITV2011", "KNG2011", "QQL2011", "BUK2011", "IKE2011", "PAA2011", "ATP2011", "NAR2008", "UUM2010"),
                     population_new =c("pop 1", "pop 2", "pop 3", "pop 4", "pop 5", "pop 6", "pop 7", "pop 8", "pop 9"))
sample_table <- read_tsv(str_c(basedir, "data/sample_table_rb.tsv")) %>%
  dplyr::select(-population_new) %>%
  left_join(rename_pop) %>%
  mutate(batch=ifelse(data_type=="se", "HiSeq-125SE", "NextSeq-150PE"))
bam_list <- read_tsv(str_c(basedir, "data/bam_list_rb.txt"), col_names = FALSE) %>%
  transmute(sample_id=str_sub(X1, 37, 47))
genome_cov_depth_ratio_filtered <- read_tsv(str_c(basedir, "results/bam_list_rb_depth_ratio_filtered.covMat"), col_names = F)[1:nrow(sample_table),1:nrow(sample_table)] %>% as.matrix
pca_table_depth_ratio_filtered <- eigen(genome_cov_depth_ratio_filtered)$vectors %>% 
  as_tibble() %>%
  transmute(PC1=V1, PC2=V2) %>%
  bind_cols(bam_list, .) %>%
  left_join(sample_table, ., by=c("sample_id_corrected"="sample_id"))
pca_table_depth_ratio_filtered %>%
  left_join(rename_pop) %>%
  ggplot(aes(x=PC1, y=PC2, color=batch)) +
  facet_wrap(~population_new) +
  geom_point() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.5))
```

![](tutorial_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

Again, because of the size of our subsetted data, the PCA result is not
substantially altered. For PCA comparison before and after batch effect
correction with our full dataset, see Figure 1B and Figure S5 of our
paper.

## DNA degradation

A major consequence of DNA degradation is deamination of cytosines
(i.e., transition of C bases into U bases), causing enrichment of C-to-T
and G-to-A substitutions in more degraded batches of data. Similar to
base quality score miscalibration, these errors will also inflate
diversity estimates, as degradation patterns will be regarded as true
variants, and can cause batch effects if two batches of data are
differentially degraded.

Here, we have samples from the same population (pop 9), some of which
are heavily degraded (determined by gDNA gel results), and some are
relatively well preserved. We sequenced these samples in two different
batches. To detect whether DNA degradation may inflate heterozygosity in
the more degraded batch of samples, we estimate heterozygosity with and
without transitions. These bam files contain a 10Mb segment of one
chromosome.

#### Heterozygosity estimation

We use ANGSD to estimate heterozygosity including and excluding
transitions.

``` bash
MINDP=2
MAXDP=10
MINMAPQ=30
MINQ=33

## NOTRANS=0: transitions are included in the analysis. NOTRANS=1: transitions are excluded from the analysis.
for NOTRANS in {0,1}; do
  for LINE in `cat $BASEDIR/data/bam_list_degradation.txt`; do
      NAME_TEMP=`echo "${LINE%.*}"`
      NAME=`echo "${NAME_TEMP##*/}"`
      echo $NAME
      OUTBASE=$NAME'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ'_notrans'$NOTRANS
      
      ## Get saf file
      $ANGSD \
      -i $LINE \
      -anc $REFERENCE \
      -out $BASEDIR/results/$OUTBASE \
      -doSaf 1 \
      -GL 1 \
      -P 8 \
      -doCounts 1 \
      -setMinDepth $MINDP \
      -setMaxDepth $MAXDP \
      -minQ $MINQ \
      -minmapq $MINMAPQ \
      -noTrans $NOTRANS

      ## Estimate sfs
      /workdir/programs/angsd0.931/angsd/misc/realSFS \
      $BASEDIR/results/$OUTBASE'.saf.idx' \
      -tole 0.0000001 \
      -P 8 \
      -seed 42 \
      > $BASEDIR/results/$OUTBASE'.ml'
  done
done
```

#### Visualization

Using R, we compare the change in heterozygosity after transitions are
excluded.

``` r
## Define path to the tutorial/ folder as basedir
basedir="/workdir/batch-effect/tutorial/"
## Read in the sample table for this analysis
sample_table <- read_tsv(str_c(basedir, "data/sample_table_degradation.tsv"))
## Loop through all samples and read in heterozygosity estimated with and without transitions
for (i in seq_len(nrow(sample_table))){
  data_type <- sample_table$data_type[i]
  if (data_type=="se"){
    prefix <- str_c(sample_table$sample_seq_id[i], "_bt2_gadMor3_sorted_dedup_realigned_subsetted_degradation_mindp2_maxdp10_minq33_minmapq30_notrans")
  } else {
    prefix <- str_c(sample_table$sample_seq_id[i], "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_subsetted_degradation_mindp2_maxdp10_minq33_minmapq30_notrans")
  }
  for (notrans in c(0, 1)) {
    ml <- read_delim(str_c(basedir, "results/", prefix, notrans, ".ml"), col_names = FALSE, delim = " ")
    het <- ml$X2/(ml$X1 + ml$X2 + ml$X3)
    line <- tibble(het=het, sample_id=sample_table$sample_id_corrected[i], data_type=sample_table$data_type[i], notrans=notrans)
    if (i==1 & notrans==0) {
      het_final <- line
    } else{
      het_final <- bind_rows(het_final, line)
    }
  }
}
## Plot
set.seed(42)
het_final %>%
  pivot_wider(names_from = notrans, values_from = het) %>%
  mutate(delta=`1`-`0`) %>%
  mutate(data_type=ifelse(data_type=="pe", "NextSeq-150PE (well-preserved)", "HiSeq-125SE (degraded)")) %>%
  mutate(data_type=factor(data_type, levels = c("NextSeq-150PE (well-preserved)", "HiSeq-125SE (degraded)"))) %>%
  ggstatsplot::ggbetweenstats(x = data_type, 
                              y = delta,  
                              type = "p", 
                              pairwise.comparisons = TRUE,
                              ggsignif.args = list(textsize = 3),
                              ggplot.component = list(coord_cartesian(ylim=c(-0.005, 0.0001)), 
                                                      scale_color_manual(values = c("#5DC863FF", "#3B528BFF")), 
                                                      theme(panel.grid = element_blank(),
                                                            axis.line = element_line()))) +
  ylab("change in heterozygosity estimates \nafter excluding transitions") +
  xlab("sample type") +
  geom_hline(yintercept = 0, linetype=2, color="red")
```

![](tutorial_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

Exclusion of transitions will always cause a decrease in estimated
heterozygosity. However, in this case, we see that the more degraded
samples are more negatively affected by this exclusion, suggesting that
the deamination of cytocines can inflate heterozygosity estimates if the
transitions are included in the more degraded samples. Therefore, to
compare batches of data with different DNA degradation levels, it is
better to exclude the transitions (while bearing in mind that the
absolute values of heterozygosity would be downward biased).

## Conclusion

In this tutorial, we used some simple bioinformatic approaches (e.g.,
visualizing the base composition of sequencing reads, comparing
heterozygosity estimates using different base quality filters and
with/without transitions, locating SNPs that have a high proportion of
low-mapping-quality reads mapped to them) to detect and mitigate four
different types of batch effects. We observed that batch effects can be
common when different lcWGS datasets are combined, but they can be
substantially reduced with simple bioinformatic methods. We conclude
that combining pre-existing lcWGS datasets is still a powerful approach
as long as batch effects are carefully accounted for.
