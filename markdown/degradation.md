DNA degradation
================

  - [Estimate heterozygositing without
    transitions](#estimate-heterozygositing-without-transitions)
      - [Run ANGSD](#run-angsd)
      - [Genome-wide, relaxed vs. stringent mapping quality filter and
        including vs. excluding
        transitions](#genome-wide-relaxed-vs-stringent-mapping-quality-filter-and-including-vs-excluding-transitions)
      - [Effect of excluding transitions in degraded vs. well-preserved
        samples](#effect-of-excluding-transitions-in-degraded-vs-well-preserved-samples)
  - [A closer look at private
    alleles](#a-closer-look-at-private-alleles)
      - [Extract private alleles and examine proportion of base
        substitutions](#extract-private-alleles-and-examine-proportion-of-base-substitutions)
  - [Exclude certain SNPs to mitigate the batch effects in PCA
    results](#exclude-certain-snps-to-mitigate-the-batch-effects-in-pca-results)
      - [Come up with a new SNP list](#come-up-with-a-new-snp-list)
      - [Run ANGSD](#run-angsd-1)
      - [PCA result with original SNP
        list](#pca-result-with-original-snp-list)
      - [Supplementary Figure: ascertainment bias when using SNP lists
        that exclude either PE or SE private
        alleles](#supplementary-figure-ascertainment-bias-when-using-snp-lists-that-exclude-either-pe-or-se-private-alleles)
      - [PCA result with a SNP list that exclude both PE and SE
        SNPs](#pca-result-with-a-snp-list-that-exclude-both-pe-and-se-snps)
  - [Effectiveness of different mitigation
    strategies](#effectiveness-of-different-mitigation-strategies)
      - [Effectiveness of base quality filtering and transition
        exclusion on the heterozygosity estimation in one population
        (UUM2010)](#effectiveness-of-base-quality-filtering-and-transition-exclusion-on-the-heterozygosity-estimation-in-one-population-uum2010)
      - [All three populations affected by
        degradation](#all-three-populations-affected-by-degradation)
      - [Effectiveness on PCA result](#effectiveness-on-pca-result)
  - [Assemble Figure 5](#assemble-figure-5)

``` r
library(tidyverse)
library(ggstatsplot)
library(cowplot)
library(ggsignif)
source("/workdir/genomic-data-analysis/scripts/individual_pca_functions.R")
```

Here, we examine the effect of DNA degradation by excluding transitions
when estimating heterozygosity. In addition, we study the composition of
private alleles in both batches of data and try to exclude these private
alleles (as well as regions affected by reference bias) from PCAs to
help alleviate batch effects.

## Estimate heterozygositing without transitions

#### Run ANGSD

``` bash
## Relaxed quality filter of 20
nohup bash /workdir/cod/greenland-cod/scripts/get_heterozygosity_notrans.sh \
/workdir/batch-effect/ \
/workdir/batch-effect/sample_lists/bam_list_realigned.txt \
/workdir/cod/reference_seqs/gadMor3.fasta \
2 \
10 \
20 \
30 \
> /workdir/batch-effect/nohups/get_heterozygosity_notrans.nohup &
## More stringent quality filter of 33
nohup nice -n 19 bash /workdir/cod/greenland-cod/scripts/get_heterozygosity_notrans.sh \
/workdir/batch-effect/ \
/workdir/batch-effect/sample_lists/bam_list_realigned.txt \
/workdir/cod/reference_seqs/gadMor3.fasta \
2 \
10 \
33 \
30 \
> /workdir/batch-effect/nohups/get_heterozygosity_notrans_stringent.nohup &

## Rerun a few outlier individuals 
## The optimization algorithm may have stuck at a local optimum for these individuals so I'll give them another try

## Relaxed quality filter of 20, excluding transitions
echo "/workdir/batch-effect/bam/ITV2011_672_55181_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/PAA2011_708_55139_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/BUK2011_029_55108_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/KNG2011_413_55159_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/UUM2010_038_55119_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/KNG2011_377_55155_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam" \
> /workdir/batch-effect/sample_lists/bam_list_realigned_rerun_heterozygosity_3.txt
nohup bash /workdir/cod/greenland-cod/scripts/get_heterozygosity_notrans.sh \
/workdir/batch-effect/ \
/workdir/batch-effect/sample_lists/bam_list_realigned_rerun_heterozygosity_3.txt \
/workdir/cod/reference_seqs/gadMor3.fasta \
2 \
10 \
20 \
30 \
> /workdir/batch-effect/nohups/rerun_heterozygosity_3.nohup &
## Stringent quality filter of 33, excluding transitions
echo "/workdir/batch-effect/bam/NAR2008_002_55156_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/ATP2011_118_14247X184_6_se_bt2_gadMor3_sorted_dedup_realigned.bam
/workdir/batch-effect/bam/IKE2011_978_55137_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/QQL2011_860_55166_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/BUK2011_045_55130_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/ITV2011_714_55194_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/UUM2010_038_55119_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam
/workdir/batch-effect/bam/QQL2011_886_55196_7_pe_bt2_gadMor3_sorted_dedup_overlapclipped_realigned.bam"\
> /workdir/batch-effect/sample_lists/bam_list_realigned_rerun_heterozygosity_4.txt
nohup nice -n 19 bash /workdir/cod/greenland-cod/scripts/get_heterozygosity_notrans.sh \
/workdir/batch-effect/ \
/workdir/batch-effect/sample_lists/bam_list_realigned_rerun_heterozygosity_4.txt \
/workdir/cod/reference_seqs/gadMor3.fasta \
2 \
10 \
33 \
30 \
> /workdir/batch-effect/nohups/rerun_heterozygosity_4.nohup &
```

#### Genome-wide, relaxed vs. stringent mapping quality filter and including vs. excluding transitions

``` r
sample_table <- read_tsv("../sample_lists/sample_table_merged.tsv")
for (i in 1:nrow(sample_table)){
  sample_seq_id <- sample_table$sample_seq_id[i]
  sample_id <- sample_table$sample_id_corrected[i]
  population <- sample_table$population[i]
  data_type <- sample_table$data_type[i]
  if (str_detect(data_type,"pe")){
    path <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_stringent <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10_minq33_minmapq30")
  } else {
    path <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_stringent <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_realigned_mindp2_maxdp10_minq33_minmapq30")
  }
  het_relaxed <- read_delim(str_c(path, ".ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, tran="Including transitions", filter="relaxed")
  het_relaxed_notrans <- read_delim(str_c(path, "_notrans.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, tran="Excluding transitions", filter="relaxed")
  het_stringent <- read_delim(str_c(path_stringent, ".ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, tran="Including transitions", filter="stringent")
  het_stringent_notrans <- read_delim(str_c(path_stringent, "_notrans.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, tran="Excluding transitions", filter="stringent")
  het_combined <- bind_rows(het_relaxed, het_relaxed_notrans, het_stringent, het_stringent_notrans)
  if(i==1){
    het_final <- het_combined
  } else {
    het_final <- bind_rows(het_final, het_combined)
  }
}
het_per_ind <- het_final %>%
  unite(col = type, tran, filter, sep = " ") %>%
  dplyr::select(sample_id, population, data_type, type, het) %>%
  pivot_wider(names_from = type, values_from = het)
## Calculate average heterozygosity while filtering out inversions
set.seed(42)
het_final %>%
  ggplot(aes(x=population, y=het)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color=data_type), height = 0, size=0.8) +
  facet_grid(filter~tran) +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8))
```

![](degradation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
het_final %>%
  ggplot(aes(x=n_sites, y=het, color=data_type)) +
  geom_point(height = 0, size=1) +
  geom_smooth(se = F, color="black", aes(group=data_type)) +
  facet_grid(filter~tran) +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8))
```

    ## Warning: Ignoring unknown parameters: height

![](degradation_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

#### Effect of excluding transitions in degraded vs. well-preserved samples

``` r
set.seed(42)
delta_het <- het_final %>%
  dplyr::select(het, sample_id, population, data_type, tran, filter) %>%
  pivot_wider(names_from = tran,  values_from = het) %>%
  mutate(delta = `Excluding transitions`-`Including transitions`,
         degradation = ifelse(population %in% c("UUM2010", "NAR2008", "ATP2011") & data_type=="se", "more degraded", "less degraded"),
         type=str_c(data_type, degradation, sep = "\n"))
delta_het %>%
  ggplot(aes(x=type, y=delta)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0) +
  facet_wrap(~filter) +
  theme_cowplot()
```

![](degradation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
set.seed(42)
delta_het %>%
  filter(filter=="stringent") %>%
  ggplot(aes(x=type, y=delta)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0) +
  #annotate(geom = "text", label="p-value=0.014", x=1.5, y=0.002) +
  scale_x_discrete(labels=c("NextSeq-150PE\n(well-preserved)\n", "HiSeq-125SE\n(well-preserved)\n", "HiSeq-125SE\n(degraded)\n")) +
  ylab("change in heterozygosity\nafter excluding transitions") +
  theme_cowplot() +
  theme(axis.title.x = element_blank())
```

![](degradation_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
## Only including the three populations where samples are split based on degradation level
delta_het %>%
  filter(filter=="stringent", population %in% c("UUM2010", "NAR2008", "ATP2011")) %>%
  ggplot(aes(x=type, y=delta)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0) +
  #annotate(geom = "text", label="p-value=0.014", x=1.5, y=0.002) +
  scale_x_discrete(labels=c("NextSeq-150PE\n(well-preserved)\n", "HiSeq-125SE\n(degraded)\n")) +
  ylab("change in heterozygosity\nafter excluding transitions") +
  ylim(c(NA, 0)) +
  theme_cowplot() +
  theme(axis.title.x = element_blank())
```

![](degradation_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
set.seed(42)
p_a <- delta_het %>%
  filter(filter=="stringent") %>%
  ggstatsplot::ggbetweenstats(x = type, 
                              y = delta,  
                              type = "p", 
                              p.adjust.method = "holm",
                              pairwise.comparisons = TRUE,
                              ggsignif.args = list(textsize = 3),
                              ggplot.component = list(coord_cartesian(ylim=c(-0.004, 0.0001)),
                                                      theme(panel.grid = element_blank(),
                                                            axis.line = element_line()))) +
  #geom_signif(comparisons = list(c("pe\nless degraded", "se\nmore degraded"), c("se\nless degraded", "se\nmore degraded")), y_position = c(-0.001, -0.0012)) +
  scale_x_discrete(labels=c("NextSeq-150PE\n(well-preserved)\n", "HiSeq-125SE\n(well-preserved)\n", "HiSeq-125SE\n(degraded)\n")) +
  ylab("change in heterozygosity estimates \nafter excluding transitions") +
  xlab("sample type") +
  geom_hline(yintercept = 0, linetype=2, color="red")
```

    ## Warning: Ignoring unknown parameters: segment.linetype

``` r
print(p_a)
```

![](degradation_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

This shows that DNA damage has an stronger effect on SE samples (which
are more degraded).

## A closer look at private alleles

#### Extract private alleles and examine proportion of base substitutions

``` r
maf_se <- read_tsv("../angsd/popminind20/se_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_popminind20.mafs.gz") %>%
  transmute(lg = chromo, position = position, major=major, minor = minor, se_maf = knownEM, se_nind=nInd)
maf_pe <- read_tsv("../angsd/popminind20/pe_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_popminind20.mafs.gz")%>%
  transmute(lg = chromo, position = position, major=major, minor = minor, pe_maf = knownEM, pe_nind=nInd)
maf_joined <- inner_join(maf_se, maf_pe) %>%
  mutate(delta = abs(se_maf- pe_maf))
private_alleles <- bind_rows((filter(maf_joined, pe_maf<0.01 | pe_maf>0.99) %>% filter(se_maf>0.1 & se_maf<0.9) %>% transmute(major=major, minor=minor, batch = "HiSeq-125SE\n(a subset of samples\nare degraded)")),
          (filter(maf_joined, se_maf<0.01 | se_maf>0.99) %>% filter(pe_maf>0.1 & pe_maf<0.9) %>% transmute(major=major, minor=minor, batch = "NextSeq-150PE\n(well-preserved)"))) %>%
  mutate(base_substitution = str_c(major, "-to-", minor)) %>%
  mutate(base_substitution = case_when(
    base_substitution %in% c("A-to-C", "T-to-G") ~ "A-to-C\nT-to-G", 
    base_substitution %in% c("A-to-G", "T-to-C") ~ "A-to-G\nT-to-C", 
    base_substitution %in% c("A-to-T", "T-to-A") ~ "A-to-T\nT-to-A", 
    base_substitution %in% c("C-to-A", "G-to-T") ~ "C-to-A\nG-to-T", 
    base_substitution %in% c("C-to-G", "G-to-C") ~ "C-to-G\nG-to-C", 
    base_substitution %in% c("C-to-T", "G-to-A") ~ "C-to-T\nG-to-A"
  )) 
p_b <- private_alleles %>%
  group_by(base_substitution, batch) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(batch) %>%
  mutate(frequency = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=base_substitution, y=frequency, fill=batch, group=batch)) +
  #geom_line() + 
  #geom_point() +
  geom_col(position = "dodge", color="black")+
  scale_fill_viridis_d(begin = 0.25, end="0.75") +
  ylim(c(NA, 0.3)) +
  labs(x="type of base substitution", y="frequency") +
  theme_ggstatsplot() +
  theme(legend.position = "top",
        panel.grid =element_blank(),
        axis.line = element_line())
print(p_b)
```

![](degradation_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
p_b_2 <- private_alleles %>%
  arrange(base_substitution) %>%
  mutate(base_substitution = factor(base_substitution)) %>%
  ggbarstats(
  x = base_substitution,
  y = batch,
  xlab = "Batch",
  legend.title = "Base substitution",
  palette = "Set2"
  )
print(p_b_2)
```

![](degradation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

C-to-T and G-to-A transitions are enriched in samples that are more
degraded

## Exclude certain SNPs to mitigate the batch effects in PCA results

Here, we try to exclude private alleles and regions affected by
reference bias in PCA to address batch effect caused by DNA degradation
and reference bias

#### Come up with a new SNP list

``` r
maf_pe <- read_tsv("../angsd/popminind20/pe_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked_popminind20.mafs.gz") %>%
  transmute(lg = chromo, position = position, major=major, minor = minor, pe_maf = knownEM, pe_nind=nInd)

maf_se <- read_tsv("../angsd/popminind20/se_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked_popminind20.mafs.gz") %>%
  transmute(lg = chromo, position = position, major=major, minor = minor, se_maf = knownEM, se_nind=nInd)

maf_joined <- inner_join(maf_pe, maf_se) %>% 
  mutate(delta=abs(se_maf-pe_maf))

maf_joined_excluding_private <- filter(maf_joined, !((pe_maf<0.01 | pe_maf>0.99)&(se_maf>0.1 & se_maf<0.9))) %>%
  filter(!((se_maf<0.01 | se_maf>0.99)&(pe_maf>0.1 & pe_maf<0.9)))

maf_excluding_pe <- maf_se %>%
  filter(!(se_maf<0.1 | se_maf>0.9))

maf_excluding_se <-  maf_pe %>%
  filter(!(pe_maf<0.1 | pe_maf>0.9))

#anymapq_depth <- read_tsv("../angsd/popminind2/bam_list_realigned_se_anymapq.pos.gz") %>%
#  rename(lg=chr, position=pos, total_depth_anymapq=totDepth)
#mapq20_depth <- read_tsv("../angsd/popminind20/se_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked_popminind20.pos.gz") %>%
#  rename(lg=chr, position=pos, total_depth_mapq20=totDepth)
#depth <- inner_join(anymapq_depth, mapq20_depth) %>%
#  mutate(depth_ratio=total_depth_mapq20/total_depth_anymapq)
#
#maf_joined_excluding_private_filtering_depth <- maf_joined_excluding_private %>%
#  left_join(depth) %>%
#  filter(depth_ratio > 0.9)
  
original_snp_list <- read_tsv("/workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.txt", col_names = c("lg", "position", "major", "minor")) 
  
new_snp_list <- semi_join(original_snp_list, maf_joined_excluding_private)
write_tsv(new_snp_list, "../angsd/global_snp_list_private_snps.txt", col_names = F)

se_snp_list <- semi_join(original_snp_list, maf_excluding_pe)
write_tsv(se_snp_list, "../angsd/global_snp_list_se_snps.txt", col_names = F)

pe_snp_list <- semi_join(original_snp_list, maf_excluding_se)
write_tsv(pe_snp_list, "../angsd/global_snp_list_pe_snps.txt", col_names = F)
```

#### Run ANGSD

``` bash
## Private SNPs, depth ratio filtered
cd /workdir/batch-effect/
/workdir/programs/angsd0.931/angsd/angsd sites index /workdir/batch-effect/angsd/global_snp_list_private_snps.txt
nohup /workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_realigned.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/bam_list_realigned_private_snps \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 16 -setMinDepth 2 -setMaxDepth 661 -minInd 2 -minQ 20 -minMapQ 20 -minMaf 0.05 \
-doIBS 2 -makematrix 1 -doCov 1 \
-sites /workdir/batch-effect/angsd/global_snp_list_private_snps.txt \
-rf /workdir/batch-effect/angsd/global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_downsampled_unlinked.chrs \
>& nohups/get_gl_bam_list_realigned_private_snps.log &
## SE SNPs
cd /workdir/batch-effect/
/workdir/programs/angsd0.931/angsd/angsd sites index /workdir/batch-effect/angsd/global_snp_list_se_snps.txt
nohup /workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_realigned.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/bam_list_realigned_se_snps \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 16 -setMinDepth 2 -setMaxDepth 661 -minInd 2 -minQ 20 -minMapQ 20 -minMaf 0.05 \
-doIBS 2 -makematrix 1 -doCov 1 \
-sites /workdir/batch-effect/angsd/global_snp_list_se_snps.txt \
-rf /workdir/cod/greenland-cod/angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.chrs \
>& nohups/get_gl_bam_list_realigned_se_snps.log &
## PE SNPs
cd /workdir/batch-effect/
/workdir/programs/angsd0.931/angsd/angsd sites index /workdir/batch-effect/angsd/global_snp_list_pe_snps.txt
nohup /workdir/programs/angsd0.931/angsd/angsd \
-b sample_lists/bam_list_realigned.txt \
-anc /workdir/cod/reference_seqs/gadMor3.fasta \
-out angsd/bam_list_realigned_pe_snps \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
-P 16 -setMinDepth 2 -setMaxDepth 661 -minInd 2 -minQ 20 -minMapQ 20 -minMaf 0.05 \
-doIBS 2 -makematrix 1 -doCov 1 \
-sites /workdir/batch-effect/angsd/global_snp_list_pe_snps.txt \
-rf /workdir/cod/greenland-cod/angsd/global_snp_list_bam_list_realigned_mincov_contamination_filtered_mindp151_maxdp661_minind102_minq20_downsampled_unlinked.chrs \
>& nohups/get_gl_bam_list_realigned_pe_snps.log &
```

#### PCA result with original SNP list

``` r
sample_table <- read_tsv("../sample_lists/sample_table_merged.tsv")
genome_cov <- read_tsv("../angsd/bam_list_realigned_downsampled_unlinked.covMat", col_names = F)[1:163,]
PCA(genome_cov, sample_table$sample_id_corrected, sample_table$data_type, 1, 2, show.ellipse = F, show.line = T)
```

    ## Warning: Computation failed in `stat_conline()`:
    ## there is no package called 'miscTools'

![](degradation_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
pca_table_data_type <- pca_table[,1:6] %>% rename(data_type=population)
PCA(genome_cov, sample_table$sample_id_corrected, sample_table$population, 1, 2, show.ellipse = F, show.line = T)
```

    ## Warning: Computation failed in `stat_conline()`:
    ## there is no package called 'miscTools'

![](degradation_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
pca_table_population <- pca_table[,1:6]

genome_dist <- read_tsv("../angsd/bam_list_realigned_downsampled_unlinked.ibsMat", col_names = F)[1:163,]
PCoA(genome_dist, sample_table$sample_id_corrected, sample_table$data_type, 10, 1, 2, show.ellipse = F, show.line = T)
```

    ## Warning: Computation failed in `stat_conline()`:
    ## there is no package called 'miscTools'

![](degradation_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
PCoA(genome_dist, sample_table$sample_id_corrected, sample_table$population, 10, 1, 2, show.ellipse = F, show.line = T)
```

    ## Warning: Computation failed in `stat_conline()`:
    ## there is no package called 'miscTools'

![](degradation_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
pca_table_data_type_summary <- group_by(pca_table_data_type, data_type) %>%
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2), PC3_mean=mean(PC3), PC4_mean=mean(PC4))
pca_table_population_summary <- group_by(pca_table, population) %>%
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2), PC3_mean=mean(PC3), PC4_mean=mean(PC4))
pca_table_combined <- left_join(pca_table_data_type, pca_table_population)

pca_table_combined %>%
  left_join(pca_table_data_type_summary) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, color=population)) +
  geom_segment(aes(x=PC1, y=PC2, xend=PC1_mean, yend=PC2_mean, color=population), size = 0.1) +
  geom_label(aes(x=PC1_mean, y=PC2_mean, label=data_type)) +
  ylim(-0.15, NA) +
  theme_cowplot()
```

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_segment).

![](degradation_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
pca_table_combined %>%
  left_join(pca_table_population_summary) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, color=data_type)) +
  geom_segment(aes(x=PC1, y=PC2, xend=PC1_mean, yend=PC2_mean, color=data_type), size = 0.1) +
  geom_label(aes(x=PC1_mean, y=PC2_mean, label=population)) +
  ylim(-0.15, NA) +
  theme_cowplot()
```

    ## Warning: Removed 2 rows containing missing values (geom_point).
    
    ## Warning: Removed 2 rows containing missing values (geom_segment).

![](degradation_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->

``` r
pca_table_combined %>%
  left_join(pca_table_data_type_summary) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, color=population)) +
  geom_segment(aes(x=PC1, y=PC2, xend=PC1_mean, yend=PC2_mean, color=population), size = 0.1) +
  #geom_label(aes(x=PC1_mean, y=PC2_mean, label=data_type)) +
  ylim(-0.15, NA) +
  facet_wrap(~data_type) +
  theme_cowplot()
```

    ## Warning: Removed 2 rows containing missing values (geom_point).
    
    ## Warning: Removed 2 rows containing missing values (geom_segment).

![](degradation_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->

``` r
pca_table_combined %>%
  left_join(pca_table_population_summary) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, color=data_type)) +
  geom_segment(aes(x=PC1, y=PC2, xend=PC1_mean, yend=PC2_mean, color=data_type), size = 0.1) +
  #geom_label(aes(x=PC1_mean, y=PC2_mean, label=population)) +
  ylim(-0.15, NA) +
  facet_wrap(~population) +
  theme_cowplot()
```

    ## Warning: Removed 2 rows containing missing values (geom_point).
    
    ## Warning: Removed 2 rows containing missing values (geom_segment).

![](degradation_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->

#### Supplementary Figure: ascertainment bias when using SNP lists that exclude either PE or SE private alleles

``` r
rename_pop <- tibble(population = c("ITV2011", "KNG2011", "QQL2011", "all pops"),
                     population_new =c("pop 1", "pop 2", "pop 3", "all pops"))
pca_combined <- bind_rows(bind_cols(pca_pe, type="NextSeq-150PE SNPs"), 
                          bind_cols(pca_se, type="HiSeq-125SE SNPs")) %>%
  mutate(batch=ifelse(data_type=="se", "HiSeq-125SE", "NextSeq-150PE"))
pca_combined_select_pops <- filter(pca_combined, population %in% c("KNG2011", "QQL2011", "ITV2011"))
pca_plot <- pca_combined_select_pops %>%
  bind_rows(mutate(pca_combined, population = "all pops")) %>%
  left_join(rename_pop) %>%
  mutate(population_new=fct_relevel(population_new, c("pop 1", "pop 2", "pop 3", "all pops"))) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_combined, color="grey", size=0.5) +
  geom_point(aes(color=batch), size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(population_new~type) +
  ylim(c(NA, 0.25)) +
  xlim(c(-0.15, NA)) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.5),
        legend.position = c(0.78, 0.94),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.x = element_text(face = "bold", size=20),
        legend.key = element_rect(fill = "white", colour = "black"))
pca_plot
```

    ## Warning: Removed 16 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

![](degradation_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

#### PCA result with a SNP list that exclude both PE and SE SNPs

``` r
genome_cov <- read_tsv("../angsd/bam_list_realigned_private_snps.covMat", col_names = F)[1:163,]
PCA(genome_cov, sample_table$sample_id_corrected, sample_table$data_type, 1, 2, show.ellipse = F, show.line = T)
```

    ## Warning: Computation failed in `stat_conline()`:
    ## there is no package called 'miscTools'

![](degradation_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
pca_table %>%
  dplyr::select(1:4) %>%
  mutate(data_type=population, population=str_sub(individual, 1, 7)) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, color=data_type)) +
  ylim(NA, 0.15) +
  facet_wrap(~population) +
  theme_cowplot()
```

![](degradation_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

## Effectiveness of different mitigation strategies

#### Effectiveness of base quality filtering and transition exclusion on the heterozygosity estimation in one population (UUM2010)

``` r
## Within KNG2011 only
set.seed(42)
het_final %>%
  filter(str_detect(sample_id, "UUM2010")) %>%
  ggplot(aes(x=data_type, y=het)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(height = 0) +
  facet_grid(filter~tran, scales = "free") +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8))
```

![](degradation_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
## PE vs SE before filtering and excluding transitions
t.test(filter(het_final, str_detect(sample_id, "UUM2010"), filter=="relaxed", tran == "Including transitions", data_type=="se")$het,
       filter(het_final, str_detect(sample_id, "UUM2010"), filter=="relaxed", tran == "Including transitions", data_type=="pe")$het)$p.value
```

    ## [1] 0.00346448

``` r
## PE vs SE after filtering and excluding transitions
t.test(filter(het_final, str_detect(sample_id, "UUM2010"), filter=="stringent", tran == "Excluding transitions", data_type=="se")$het,
       filter(het_final, str_detect(sample_id, "UUM2010"), filter=="stringent", tran == "Excluding transitions", data_type=="pe")$het)$p.value
```

    ## [1] 0.001337545

This shows that after excluding transitions and using stringent depth
and quality filters, batch effects on heterozygosity estimate are
somewhat reduced.

#### All three populations affected by degradation

``` r
rename_pop <- tibble(population = c("ATP2011", "NAR2008", "UUM2010"),
                     population_new =c("pop 4", "pop 5", "pop 6"))
set.seed(42)
p_c <- het_final %>%
  left_join(rename_pop) %>%
  filter(!is.na(population_new), filter=="stringent") %>%
  mutate(batch=ifelse(data_type=="pe", "well-preserved", "degraded")) %>%
  mutate(type=ifelse(tran=="Including transitions", "Before", "After")) %>%
  mutate(type=fct_relevel(type, c("Before", "After"))) %>%
  ggplot(aes(x=population, y=het)) +
  geom_boxplot(outlier.alpha = 0, color="black", size=0.2, width=0.2) +
  geom_jitter(aes(color=batch), height = 0, width = 0.1, size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  ylim(c(0,NA)) +
  ylab("heterozygosity") +
  facet_grid(population_new~type, scales = "free_y") +
  xlab(" ") +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        legend.position = c(0.75, 0.93),
        legend.text = element_text(size=9),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.x = element_text(face = "bold", size=20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"))
print(p_c)
```

![](degradation_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

#### Effectiveness on PCA result

``` r
genome_cov <- read_tsv("../angsd/bam_list_realigned_downsampled_unlinked.covMat", col_names = F) %>%
  as.matrix()
PCA(genome_cov, sample_table$sample_id_corrected, sample_table$data_type, 1, 2, show.ellipse = F, show.line = F)
```

``` r
pca_before <- pca_table[,1:4] %>%
  rename(data_type=population) %>%
  left_join(transmute(sample_table, individual=sample_id_corrected, population=population))

genome_cov <- read_tsv("../angsd/bam_list_realigned_private_snps.covMat", col_names = F) %>%
  as.matrix()
PCA(genome_cov, sample_table$sample_id_corrected, sample_table$data_type, 1, 2, show.ellipse = F, show.line = F)
```

``` r
pca_after <- pca_table[,1:4] %>%
  rename(data_type=population) %>%
  left_join(transmute(sample_table, individual=sample_id_corrected, population=population))
pca_combined <- bind_rows(bind_cols(pca_before, type="Before"), 
                          bind_cols(pca_after, type="After")) %>%
  mutate(type=fct_relevel(type, c("Before", "After"))) %>%
  mutate(batch=ifelse(data_type=="se", "degraded", "well-preserved"))
pca_combined_select_pops <- filter(pca_combined, population %in% c("ATP2011", "NAR2008", "UUM2010"))
```

``` r
p_d <- pca_combined_select_pops %>%
  left_join(rename_pop) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_combined, color="grey", size=0.5) +
  geom_point(aes(color=batch), size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(population_new~type) +
  ylim(-0.1, NA) +
  xlim(-0.15, NA) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.5),
        legend.position = c(0.75, 0.93),
        legend.text = element_text(size=9),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.x = element_text(face = "bold", size=20),
        legend.key = element_rect(fill = "white", colour = "black"))
p_d
```

    ## Warning: Removed 12 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

![](degradation_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## Assemble Figure 5

``` r
bottom <- cowplot::plot_grid(p_b, p_a,nrow = 1, labels = c("C", "D"), scale = 0.9)
top <- cowplot::plot_grid(p_c, p_d, nrow = 1, labels = c("A", "B"), scale = 0.9)
```

    ## Warning: Removed 12 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

``` r
figure <- cowplot::plot_grid(top, bottom, nrow = 2, rel_heights = c(5.3, 7))
print(figure)
```

![](degradation_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->