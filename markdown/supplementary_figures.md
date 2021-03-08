Supplementary Figures
================

``` r
library(tidyverse)
library(cowplot)
library(RcppCNPy)
source("/workdir/genomic-data-analysis/scripts/individual_pca_functions.R")
```

## Supplementary Figure: heterozygosity with polyG trimming only

``` r
sample_table <- read_tsv("../sample_lists/sample_table_merged.tsv")
rename_pop <- tibble(population = c("ITV2011", "KNG2011", "QQL2011", "all pops"),
                     population_new =c("pop 1", "pop 2", "pop 3", "all pops"))
for (i in 1:nrow(sample_table)){
  sample_seq_id <- sample_table$sample_seq_id[i]
  sample_id <- sample_table$sample_id_corrected[i]
  population <- sample_table$population[i]
  data_type <- sample_table$data_type[i]
  if (str_detect(data_type,"pe")){
    path <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_overlapclipped_realigned_mindp2_maxdp10_minq20_minmapq30")
  } else {
    path <- str_c("../angsd/heterozygosity/", sample_seq_id,  "_bt2_gadMor3_sorted_dedup_realigned_mindp2_maxdp10_minq20_minmapq30")
  }
  het_relaxed <- read_delim(str_c(path, "_notrans.ml"), col_names = F, delim = " ") %>% 
    transmute(n_sites=(X1+X2+X3), n_snp=X2, het=n_snp/n_sites) %>%
    mutate(sample_id=sample_id, population=population, data_type=data_type, type="Before")
  if(i==1){
    het_final <- het_relaxed
  } else {
    het_final <- bind_rows(het_final, het_relaxed)
  }
}
het_gg <- het_final %>%
  left_join(rename_pop) %>%
  mutate(batch=ifelse(data_type=="pe", "NextSeq-150PE", "HiSeq-125SE")) 
set.seed(42)
het_plot <- het_gg %>%
  filter(population %in% c("KNG2011", "QQL2011", "ITV2011")) %>%
  ggplot(aes(x=population, y=het)) +
  geom_boxplot(outlier.alpha = 0, color="black", size=0.2, width=0.2) +
  geom_jitter(aes(color=batch), height = 0, width = 0.1, size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  ylab("heterozygosity") +
  facet_grid(population_new~., scales = "free_y") +
  xlab(" ") +
  scale_y_continuous(limits = c(0, 0.005), breaks = 0.001*(0:5), labels = c("0", "0.001", "0.002", "0.003", "0.004", "0.005")) +
  coord_flip() +
  theme_cowplot() +
  theme(panel.background=element_rect(colour="black", size=0.8),
        legend.position = c(0.55, 0.94),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.x = element_text(face = "bold", size=20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.key = element_rect(fill = "white", colour = "black"))
het_plot
```

![](supplementary_figures_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Supplementary Figure: PCA with PCAngsd

``` r
pca_angsd_select_pops <- filter(pca_pcangsd, population %in% c("KNG2011", "QQL2011", "ITV2011"))
pca_plot <- pca_angsd_select_pops %>%
  left_join(rename_pop) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_pcangsd, color="grey", size=0.5) +
  geom_point(aes(color=batch), size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(population_new~.) +
  ylim(NA, 0.2) +
  xlim(-0.15, NA) +
  theme_cowplot() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour="black",size=0.5),
        legend.position = c(0.55, 0.94),
        legend.key.size = unit(0.5, 'lines'),
        strip.text.x = element_text(face = "bold", size=20),
        legend.key = element_rect(fill = "white", colour = "black"))
pca_plot
```

    ## Warning: Removed 6 rows containing missing values (geom_point).

![](supplementary_figures_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Supplementary Figure: PCA with ascertainment bias

``` r
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

![](supplementary_figures_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->