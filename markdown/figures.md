Figures
================

``` r
library(tidyverse)
```

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'pillar'

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'tibble'

    ## Warning: replacing previous import 'lifecycle::last_warnings' by
    ## 'rlang::last_warnings' when loading 'hms'

``` r
library(cowplot)
library(RcppCNPy)
source("/workdir/genomic-data-analysis/scripts/individual_pca_functions.R")
```

## Figure 1

#### Heterozygosity

Here, I am comparing heterozygosity estimated from:

1.  polyG trimmed, relaxed read quality filtering of 20 (before)

2.  sliding window trimmed, stringent read quality filtering of 33
    (after)

Neither have transitions filtered, have a stringent maximum depth filter
of 10, and a stringent mapping quality filter of 30.

``` r
sample_table <- read_tsv("../sample_lists/sample_table_merged.tsv")
rename_pop <- tibble(population = c("ITV2011", "KNG2011", "QQL2011", "BUK2011", "IKE2011", "PAA2011", "ATP2011", "NAR2008", "UUM2010"),
                     population_new =c("pop 1", "pop 2", "pop 3", "pop 4", "pop 5", "pop 6", "pop 7", "pop 8", "pop 9"))
for (i in 1:nrow(sample_table)){
  sample_seq_id <- sample_table$sample_seq_id[i]
  sample_id <- sample_table$sample_id_corrected[i]
  population <- sample_table$population[i]
  data_type <- sample_table$data_type[i]
  if (str_detect(data_type,"pe")){
    path_original <- str_c("../../cod/greenland-cod/angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_minq20_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_new <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10_minq33_minmapq30")
  } else {
    path_original <- str_c("../../cod/greenland-cod/angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_minq20_sorted_dedup_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_new <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_realigned_mindp2_maxdp10_minq33_minmapq30")
  }
  het_original <- read_delim(str_c(path_original, ".ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="Before")
  het_stringent <- read_delim(str_c(path_new, ".ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="After")
  het_combined <- bind_rows(het_original, het_stringent)
  if(i==1){
    het_final <- het_combined
  } else {
    het_final <- bind_rows(het_final, het_combined)
  }
}
het_gg <- het_final %>%
  left_join(rename_pop) %>%
  mutate(type=fct_relevel(type, c("Before", "After"))) %>%
  mutate(batch=ifelse(data_type=="pe", "NextSeq-150PE", "HiSeq-125SE")) 
set.seed(42)
het_plot <- het_gg %>%
  filter(population_new %in% str_c("pop ", 1:6)) %>%
  ggplot(aes(x="", y=het*10^3)) +
  geom_boxplot(outlier.alpha = 0, color="black", size=0.2, width=0.2) +
  #geom_jitter(data=(het_gg %>% dplyr::select(-population_new) %>% filter(! population %in% c("KNG2011", "QQL2011", "ITV2011"))) , color="grey", height = 0, width = 0.3, size=1) +
  geom_jitter(aes(color=batch), height = 0, width = 0.1, size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  ylab(expression(paste("heterozygosity (in ", 10^-3, ")"))) +
  facet_grid(type~population_new, scales = "free_y") +
  xlab(" ") +
  #scale_y_continuous(limits = c(0.002, 0.008), breaks = 0.002*(1:4)) +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        legend.position = c(0.9, 0.09),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.y = element_text(face = "bold", size=20),
        strip.text.x = element_text(size=18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
het_plot
```

![](figures_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

#### PCA

Here, I am comparing PCA result from

1.  the original LD pruned SNP list

2.  a subsetted SNP list with SNPs with depth ratio \< 0.9 filtered out

Both have gone through sliding window trimming, have base quality and
mapping quality filters of 20, minimum individual filter of 20, and etc.
Both PCA are performed through ANGSD’s `-doCov 1`.

Also, note that two outlier individuals from UUM2010 are removed from
these plots for clearer results.

``` r
pca_combined <- bind_rows(bind_cols(pca_before, type="Before"), 
                          bind_cols(pca_after, type="After")) %>%
  mutate(type=fct_relevel(type, c("Before", "After"))) %>%
  mutate(batch=ifelse(data_type=="se", "HiSeq-125SE", "NextSeq-150PE")) %>%
  filter(! individual %in% c("UUM2010_036", "UUM2010_038"))
pca_plot <- pca_combined %>%
  left_join(rename_pop) %>%
  filter(population_new %in% str_c("pop ", 1:6)) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_combined, color="grey", size=0.5) +
  geom_point(aes(color=batch), size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(type~population_new) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.8),
        legend.position = c(0.9, 0.09),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.y = element_text(face = "bold", size=20),
        strip.text.x = element_text(size=18),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
pca_plot
```

![](figures_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

#### Fst

Here, I am plotting Fst estimated from sliding window trimmed data with
relaxed mapping quality and read quality filter (both 20). There is a
minimum number of individuals filter of 20 for both batches of data. I
compare the Fst result before and after applying a depth ratio filter of
0.9.

``` r
maf_se <- read_tsv("../angsd/popminind20/se_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_popminind20.mafs.gz") %>%
  transmute(lg = chromo, position = position, major=major, minor = minor, se_maf = knownEM, se_nind=nInd)
maf_pe <- read_tsv("../angsd/popminind20/pe_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_popminind20.mafs.gz")%>%
  transmute(lg = chromo, position = position, major=major, minor = minor, pe_maf = knownEM, pe_nind=nInd)
fst <- read_tsv("../angsd/popminind20/pe_se_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_popminind20.alpha_beta.txt", col_names = F) %>%
  mutate(X5=X3/X4) %>%
  transmute(lg=X1, position = X2, alpha=X3, beta=X4, fst = X5)
anymapq_depth <- read_tsv("../angsd/popminind2/bam_list_realigned_se_anymapq.pos.gz") %>%
  rename(lg=chr, position=pos, total_depth_anymapq=totDepth)
mapq20_depth <- read_tsv("../angsd/popminind20/se_global_snp_list_bam_list_realigned_mindp46_maxdp184_minind20_minq20_popminind20.pos.gz") %>%
  rename(lg=chr, position=pos, total_depth_mapq20=totDepth)
maf_joined <- inner_join(maf_se, maf_pe) %>%
  left_join(fst) %>%
  filter(str_detect(lg, "LG"))
depth <- inner_join(anymapq_depth, mapq20_depth) %>%
  mutate(depth_ratio=total_depth_mapq20/total_depth_anymapq)
fst_combined <- maf_joined %>%
  left_join(depth) %>%
  filter(depth_ratio > 0.9) %>%
  mutate(type="After") %>%
  bind_rows(maf_joined %>% left_join(depth) %>% mutate(type="Before")) %>%
  mutate(type=fct_relevel(type, c("Before", "After")))
## per SNP
fst_plot <- fst_combined %>%
  ggplot(aes(x=position/10^6, y=fst, color=lg)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = rep(c("black", "darkgrey"), 12)) +
  scale_x_continuous(breaks = c(0,15)) +
  facet_grid(type~lg, scales = "free_x", space = "free_x") +
  xlab("position (in Mb)") +
  ylab(expression(F[ST]~between~batches)) +
  theme_cowplot() +
  theme(panel.spacing = unit(0.0, "lines")) +
  theme(legend.position = "none",
        strip.text.y = element_text(face = "bold", size=20))
fst_plot
```

![](figures_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

#### Combine the plots

``` r
figure_1 <- plot_grid(het_plot, pca_plot, fst_plot, labels = c('A', 'B', 'C'), label_size = 20, nrow = 3, rel_heights = c(4, 3.8, 4))
figure_1
```

![](figures_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave("../figures/figure_1.png", figure_1, width = 15, height=11.8, unit="in", dpi="retina")
```

## Cover image

#### Option 1

``` r
panel_a <- het_gg %>%
  filter(population_new %in% str_c("pop ", 2)) %>%
  ggplot(aes(x="", y=het*10^3)) +
  geom_boxplot(outlier.alpha = 0, color="black", size=0.2, width=0.2) +
  geom_jitter(aes(color=batch), height = 0, width = 0.1, size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  ylab(expression(paste("heterozygosity (in ", 10^-3, ")"))) +
  facet_grid(type~population_new, scales = "free_y") +
  xlab(" ") +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        legend.position = c(0.52, 0.93),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.y = element_text(face = "bold", size=20),
        strip.text.x = element_text(size=18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
panel_b <- pca_combined %>%
  left_join(rename_pop) %>%
  filter(population_new %in% str_c("pop ", 2)) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_combined, color="grey", size=0.5) +
  geom_point(aes(color=batch), size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(type~population_new) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.8),
        legend.position = c(0.52, 0.93),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.y = element_text(face = "bold", size=20),
        strip.text.x = element_text(size=18),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
panel_c <- fst_combined %>%
  filter(lg=="LG06") %>%
  ggplot(aes(x=position/10^6, y=fst, color=lg)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = rep(c("black", "darkgrey"), 12)) +
  facet_grid(type~lg, scales = "free_x", space = "free_x") +
  xlab("position (in Mb)") +
  ylab(expression(F[ST]~between~batches)) +
  ylim(c(0, 0.85)) +
  theme_cowplot() +
  theme(panel.spacing = unit(0.5, "lines")) +
  theme(panel.border = element_rect(colour="black",size=0.8),
        legend.position = "none",
        strip.text.x = element_text(size=18),
        strip.text.y = element_text(face = "bold", size=20))
plot_grid(panel_b, nrow = 2, rel_heights = c(1, 0.06)) %>%
  plot_grid(panel_a, ., labels = c('A', 'B'), label_size = 20, nrow = 1, rel_widths = c(1, 1)) %>%
  plot_grid(panel_c, labels=c(NA, 'C'), label_size = 20, nrow = 2, rel_heights = c(1.5, 1))
```

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](figures_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

#### Option B

``` r
panel_a <- het_gg %>%
  filter(population_new %in% str_c("pop ", 1:3)) %>%
  ggplot(aes(x="", y=het*10^3)) +
  geom_boxplot(outlier.alpha = 0, color="black", size=0.2, width=0.2) +
  geom_jitter(aes(color=batch), height = 0, width = 0.1, size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  ylab(expression(paste("heterozygosity (in ", 10^-3, ")"))) +
  facet_grid(type~population_new, scales = "free_y") +
  xlab(" ") +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        legend.position = c(0.78, 0.09),
        legend.key.size = unit(0.8, 'lines'),
        strip.text.y = element_text(face = "bold", size=20),
        strip.text.x = element_text(size=18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
panel_b <- pca_combined %>%
  left_join(rename_pop) %>%
  filter(population_new %in% str_c("pop ", 1:3)) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_combined, color="grey", size=0.5) +
  geom_point(aes(color=batch), size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(type~population_new) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.8),
        legend.position = c(0.78, 0.42),
        legend.key.size = unit(0.8, 'lines'),
        strip.text.y = element_text(face = "bold", size=20),
        strip.text.x = element_text(size=18),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())

cover_image <- plot_grid(panel_a, panel_b, labels = c('A', 'B'), label_size = 20, nrow = 2, rel_heights = c(4.2, 3.8))
cover_image
```

![](figures_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave("../misc/cover_image.pdf", plot=cover_image, width = 8, height = 8, units = "in")
```

## Figure S3

This is a sequential visualization of estimated heterozygosity after
each issue is taken care of

``` r
sample_table <- read_tsv("../sample_lists/sample_table_merged.tsv")
rename_pop <- tibble(population = c("ITV2011", "KNG2011", "QQL2011", "BUK2011", "IKE2011", "PAA2011", "ATP2011", "NAR2008", "UUM2010"),
                     population_new =c("pop 1", "pop 2", "pop 3", "pop 4", "pop 5", "pop 6", "pop 7", "pop 8", "pop 9"))
for (i in 1:nrow(sample_table)){
  sample_seq_id <- sample_table$sample_seq_id[i]
  sample_id <- sample_table$sample_id_corrected[i]
  population <- sample_table$population[i]
  data_type <- sample_table$data_type[i]
  if (str_detect(data_type,"pe")){
    path_original <- str_c("../../cod/greenland-cod/angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_minq20_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_new <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10")
  } else {
    path_original <- str_c("../../cod/greenland-cod/angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_minq20_sorted_dedup_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_new <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_realigned_mindp2_maxdp10")
  }
  het_original <- read_delim(str_c(path_original, ".ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="Before mitigation")
  
  het_relaxed <- read_delim(str_c(path_new, "_minq20_minmapq30.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="After sliding-window trimming")
  
  het_stringent <- read_delim(str_c(path_new, "_minq33_minmapq30.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="After stringent quality filtering")

  het_stringent_notrans <- read_delim(str_c(path_new, "_minq33_minmapq30_notrans.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="After filtering transitions")
  
  het_combined <- bind_rows(het_original, het_relaxed, het_stringent, het_stringent_notrans)
  if(i==1){
    het_final <- het_combined
  } else {
    het_final <- bind_rows(het_final, het_combined)
  }
}
het_gg <- het_final %>%
  left_join(rename_pop) %>%
  mutate(type=fct_relevel(type, c("Before mitigation", "After sliding-window trimming", "After stringent quality filtering", "After filtering transitions"))) %>%
  mutate(batch=ifelse(data_type=="pe", "NextSeq-150PE", "HiSeq-125SE")) 
set.seed(42)
het_plot <- het_gg %>%
  ggplot(aes(x="", y=het*10^3)) +
  geom_boxplot(outlier.alpha = 0, color="black", size=0.2, width=0.2) +
  #geom_jitter(data=(het_gg %>% dplyr::select(-population_new) %>% filter(! population %in% c("KNG2011", "QQL2011", "ITV2011"))) , color="grey", height = 0, width = 0.3, size=1) +
  geom_jitter(aes(color=batch), height = 0, width = 0.1, size=1.5) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  ylab(expression(paste("heterozygosity (in ", 10^-3, ")"))) +
  facet_grid(population_new~type, scales = "free_y") +
  xlab(" ") +
  #scale_y_continuous(limits = c(0.002, 0.008), breaks = 0.002*(1:4)) +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        legend.position = c(0.84, 0.98),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
het_plot
```

![](figures_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Figure S4

This compares filtering based methods with base quality score
recalibration.

``` r
sample_table <- read_tsv("../sample_lists/sample_table_merged.tsv")
rename_pop <- tibble(population = c("ITV2011", "KNG2011", "QQL2011", "BUK2011", "IKE2011", "PAA2011", "ATP2011", "NAR2008", "UUM2010"),
                     population_new =c("pop 1", "pop 2", "pop 3", "pop 4", "pop 5", "pop 6", "pop 7", "pop 8", "pop 9"))
for (i in 1:nrow(sample_table)){
  sample_seq_id <- sample_table$sample_seq_id[i]
  sample_id <- sample_table$sample_id_corrected[i]
  population <- sample_table$population[i]
  data_type <- sample_table$data_type[i]
  if (str_detect(data_type,"pe")){
    path_original <- str_c("../../cod/greenland-cod/angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_minq20_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_new <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10")
  } else {
    path_original <- str_c("../../cod/greenland-cod/angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_minq20_sorted_dedup_realigned_mindp2_maxdp10_minq20_minmapq30")
    path_new <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_realigned_mindp2_maxdp10")
  }
  het_relaxed <- read_delim(str_c(path_new, "_minq20_minmapq30.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="Relaxed base quality filter")
  
  het_stringent <- read_delim(str_c(path_new, "_minq33_minmapq30.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="Stringent base quality filter")

  het_bqsr_angsd <- read_delim(str_c(path_new, "_minq0_minmapq30_bqsr.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="BQSR with ANGSD")

  het_bqsr_gatk <- read_delim(str_c(path_new, "_minq0_minmapq30_bqsr_gatk.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="BQSR with GATK")

  het_combined <- bind_rows(het_relaxed, het_stringent, het_bqsr_angsd, het_bqsr_gatk)
  if(i==1){
    het_final <- het_combined
  } else {
    het_final <- bind_rows(het_final, het_combined)
  }
}
bqsr_gg <- het_final %>%
  left_join(rename_pop) %>%
  mutate(type=fct_relevel(type, c("Relaxed base quality filter", "Stringent base quality filter", "BQSR with ANGSD", "BQSR with GATK"))) %>%
  mutate(batch=ifelse(data_type=="pe", "NextSeq-150PE", "HiSeq-125SE")) 
set.seed(42)
bqsr_plot <- bqsr_gg %>%
  filter(!population_new %in% c("pop 7", "pop 8", "pop 9")) %>%
  ggplot(aes(x=population_new, y=het*10^3)) +
  geom_boxplot(outlier.alpha = 0, color="black") +
  geom_jitter(aes(color=batch), height = 0, size=1.5) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  scale_x_discrete(limits = rev) +
  ylab(expression(paste("heterozygosity (in ", 10^-3, ")"))) +
  facet_grid(type~., scales = "free_y") +
  xlab(" ") +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        legend.position = c(0.03, 0.98),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
bqsr_plot
```

![](figures_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Figure S5

``` r
pca_combined <- bind_rows(bind_cols(pca_before, type="Full LD-pruned\nSNP list"), 
                          bind_cols(pca_pe, type="NextSeq-150PE\nSNPs only"), 
                          bind_cols(pca_se, type="HiSeq-125SE\nSNPs only"), 
                          bind_cols(pca_private_filtered, type="All private\nSNPs removed"), 
                          bind_cols(pca_depth_ratio_filtered, type="Depth ratio filtered\nSNP list"), 
                          bind_cols(pca_depth_ratio_filtered_notrans, type="Depth ratio & transition\nfiltered SNP list")) %>%
  mutate(type=fct_relevel(type, c("Full LD-pruned\nSNP list", "NextSeq-150PE\nSNPs only", "HiSeq-125SE\nSNPs only", "All private\nSNPs removed", "Depth ratio filtered\nSNP list", "Depth ratio & transition\nfiltered SNP list"))) %>%
  mutate(batch=ifelse(data_type=="se", "HiSeq-125SE", "NextSeq-150PE")) %>%
  filter(! individual %in% c("UUM2010_036", "UUM2010_038"))
pca_plot <- bind_rows(pca_combined, mutate(pca_combined, population="all pops")) %>%
  left_join(rename_pop) %>%
  mutate(population_new=ifelse(is.na(population_new), "all pops", population_new)) %>% 
  mutate(population_new=fct_relevel(population_new, c(str_c("pop ", 1:9), "all pops"))) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_combined, color="grey", size=0.3) +
  geom_point(aes(color=batch), size=1.5) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(population_new~type) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.8),
        legend.position = c(0.89, 0.99),
        legend.key.size = unit(0.5, 'lines'),
        legend.text = element_text(size=8),
        strip.text.y = element_text(size=10),
        strip.text.x = element_text(size=10),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_blank())
pca_plot
```

![](figures_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
