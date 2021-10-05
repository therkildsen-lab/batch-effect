Sequencing depth
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
library(RcppCNPy)
library(cowplot)
library(scales)
library(ggrepel)
library(MASS)
source("/workdir/genomic-data-analysis/scripts/individual_pca_functions.R")
```

Here, we examine the effect of variation in sequencing depth by running
PCA on simulated data. We also run PCAngsd on the real cod data to
demonstrate its potential vulnerability to batch effects caused by
variable sequencing depth.

The simulation was done as part of the [lcWGS guide
project](https://github.com/therkildsen-lab/lcwgs-simulation)
(specifically see
[link](https://github.com/therkildsen-lab/lcwgs-simulation/blob/master/markdowns/simulation_workflow_spatial_pop_sim.md#uneven-coverage)),
and the simulation output are store in the directory
`/workdir/lcwgs-simulation` on the `nt246` server.

I will just include the visualization script below.

## Analysis with simulated data

#### PCA with PCAngsd

``` r
i=1
for (sample_size in c(5,10,20,40,80)){
  pop_label <- read_lines(paste0("/workdir/lcwgs-simulation/spatial_pop_sim/rep_1/sample_lists/bam_list_",sample_size,"_uneven_coverage.txt")) %>%
    str_extract('p[1-9]')
  coverage <- ifelse(1:(sample_size*9) %% 2 == 1,0.125, 4)
  ## Read covariance matrix
  cov_matrix <- npyLoad(paste0("/workdir/lcwgs-simulation/spatial_pop_sim/rep_1/angsd/pcagnsd_bam_list_",sample_size,"_uneven_coverage.cov.npy")) %>%
    as.matrix()
  ## Perform eigen decomposition
  e <- eigen(cov_matrix)
  e_value<-e$values
  x_variance<-e_value[1]/sum(e_value)*100
  y_variance<-e_value[2]/sum(e_value)*100
  e_vector <- as.data.frame(e$vectors)[,1:5]
  pca_table <- bind_cols(pop_label=pop_label, e_vector) %>%
    transmute(population=pop_label, PC1=rescale(V1, c(-1, 1)), PC2=rescale(V2, c(-1, 1)), PC3=rescale(V3, c(-1, 1)), PC4=rescale(V4, c(-1, 1)), PC5=rescale(V5, c(-1, 1)), coverage=coverage, sample_size=sample_size)
  ## Bind PCA tables and DAPC tables for all sample size and coverage combinations
  if (i==1){
    pca_table_final <- pca_table
  } else {
    pca_table_final <- bind_rows(pca_table_final,pca_table)
  }
  i=i+1
}
ggplot(pca_table_final,aes(x=PC1, y=PC2, color=population, shape=as.character(coverage))) +
  geom_point() +
  facet_grid(.~sample_size, scales="free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
```

![](depth_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
p1 <- pca_table_final %>%
  mutate(method="PCAngsd") %>%
  ggplot(aes(x=PC1, y=PC2, color=as.character(coverage))) +
  geom_point() +
  facet_grid(method~sample_size, scales="free") +
  scale_color_viridis_d(name = "coverage", option = "D", begin = 0.75, end = 0.25) +
  theme_cowplot() +
  theme(panel.border = element_rect(size = 1, color = "black"),
        text = element_text(size=11), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        legend.position = "top",
        legend.justification = "center")
p1
```

![](depth_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

### PCA with covMat

``` r
i=1
for (sample_size in c(5,10,20,40,80)){
  pop_label <- read_lines(paste0("/workdir/lcwgs-simulation/spatial_pop_sim/rep_1/sample_lists/bam_list_",sample_size,"_uneven_coverage.txt")) %>%
    str_extract('p[1-9]')
  coverage <- ifelse(1:(sample_size*9) %% 2 == 1,0.125, 4)
  ## Read covariance matrix
  cov_matrix <- read_tsv(paste0("/workdir/lcwgs-simulation/spatial_pop_sim/rep_1/angsd/bam_list_",sample_size,"_uneven_coverage.covMat"), col_names = F) %>%
    as.matrix() %>%
    .[,-(sample_size*9+1)]
  cov_matrix[is.na(cov_matrix)]<- median(cov_matrix, na.rm = T)
  ## Perform eigen decomposition
  e <- eigen(cov_matrix)
  e_value<-e$values
  x_variance<-e_value[1]/sum(e_value)*100
  y_variance<-e_value[2]/sum(e_value)*100
  e_vector <- as.data.frame(e$vectors)[,1:5]
  pca_table <- bind_cols(pop_label=pop_label, e_vector) %>%
    transmute(population=pop_label, PC1=rescale(V1, c(-1, 1)), PC2=rescale(V2, c(-1, 1)), PC3=rescale(V3, c(-1, 1)), PC4=rescale(V3, c(-1, 1)), PC5=rescale(V5, c(-1, 1)), sample_size=sample_size, coverage=coverage)
  ## Bind PCA tables and DAPC tables for all sample size and coverage combinations
  if (i==1){
    pca_table_final <- pca_table
  } else {
    pca_table_final <- bind_rows(pca_table_final,pca_table)
  }
  i=i+1
}
ggplot(pca_table_final,aes(x=PC1, y=PC2, color=population, shape=as.character(coverage))) +
  geom_point() +
  facet_grid(.~sample_size, scales="free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
```

![](depth_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
p2 <- pca_table_final %>%
  mutate(method="ANGSD (-doCov 1)") %>%
  ggplot(aes(x=PC1, y=PC2, color=as.character(coverage))) +
  geom_point() +
  scale_color_viridis_d(name = "coverage", option = "D", begin = 0.75, end = 0.25) +
  facet_grid(method~sample_size, scales="free") +
  theme_cowplot() +
  theme(panel.border = element_rect(size = 1, color = "black"),
        text = element_text(size=11), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        legend.position = "none")
p2
```

![](depth_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

### PCoA with ibsMat

``` r
i=1
for (sample_size in c(5,10,20,40,80)){
  pop_label <- read_lines(paste0("/workdir/lcwgs-simulation/spatial_pop_sim/rep_1/sample_lists/bam_list_",sample_size,"_uneven_coverage.txt")) %>%
    str_extract('p[1-9]')
  coverage <- ifelse(1:(sample_size*9) %% 2 == 1,0.125, 4)
  ## Read covariance matrix
  dist_matrix <- read_tsv(paste0("/workdir/lcwgs-simulation/spatial_pop_sim/rep_1/angsd/bam_list_",sample_size,"_uneven_coverage.ibsMat"), col_names = F) %>%
    as.matrix() %>%
    .[,-(sample_size*9+1)]
  dist_matrix[is.na(dist_matrix)] <- mean(dist_matrix, na.rm = T)
  ## Perform MDS
  mds <- cmdscale(as.dist(dist_matrix), k=5) %>%
    as.data.frame() 
  mds_table <- bind_cols(pop_label=pop_label, mds) %>%
    transmute(population=pop_label, PCo1=rescale(V1, c(-1, 1)), PCo2=rescale(V2, c(-1, 1)), PCo3=rescale(V3, c(-1, 1)), PCo4=rescale(V4, c(-1, 1)), PCo5=rescale(V5, c(-1, 1)), coverage=coverage, sample_size=sample_size)
  eigen_value <- cmdscale(as.dist(dist_matrix), k=5, eig = T)$eig
  var_explained <- round(eigen_value/sum(eigen_value)*100, 2)
  ## Bind PCoA tables and DAPC tables for all sample size and coverage combinations
  if (i==1){
    mds_table_final <- mds_table
  } else {
    mds_table_final <- bind_rows(mds_table_final,mds_table)
  }
  i=i+1
}
ggplot(mds_table_final,aes(x=PCo1, y=PCo2, color=population)) +
  geom_point() +
  facet_grid(.~sample_size, scales="free") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
```

![](depth_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
p3 <- mds_table_final %>%
  mutate(method="ANGSD (-doIBS 2)") %>%
  ggplot(aes(x=PCo1, y=PCo2, color=as.character(coverage))) +
  geom_point() +
  facet_grid(method~sample_size, scales="free") +
  scale_color_viridis_d(name = "coverage", option = "D", begin = 0.75, end = 0.25) +
  theme_cowplot() +
  theme(panel.border = element_rect(size = 1, color = "black"),
        text = element_text(size=11), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        legend.position = "none")
p3
```

![](depth_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

### Assemble plots for batch effect paper

``` r
cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.2, 1, 1), labels=c("A", "B", "C"), label_y=c(0.8, 1, 1))
```

![](depth_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Test PCAngsd with real cod data using the filtered SNP list

``` bash
nohup python2 /workdir/programs/pcangsd/pcangsd.py \
-beagle /workdir/batch-effect/angsd/bam_list_realigned_private_snps.beagle.gz \
-minMaf 0.05 \
-threads 8 \
-o /workdir/batch-effect/angsd/bam_list_realigned_private_snps_pcangsd \
> /workdir/batch-effect/nohups/run_pcangsd_private_snps.nohup &
```

#### Supplementary Figure: PCA with PCAngsd

``` r
pca_angsd_select_pops <- filter(pca_pcangsd, population %in% c("KNG2011", "QQL2011", "ITV2011"))
pca_plot <- pca_angsd_select_pops %>%
  left_join(rename_pop) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(data=pca_pcangsd, color="grey", size=0.5) +
  geom_point(aes(color=batch), size=2) +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  facet_grid(population_new~.) +
  ylim(NA, 0.15) +
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

![](depth_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->
