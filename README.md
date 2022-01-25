[![DOI](https://zenodo.org/badge/308741089.svg)](https://zenodo.org/badge/latestdoi/308741089)

## Batch effects in population genomic studies with low-coverage whole genome sequencing data: causes, detection, and mitigation

This is a GitHub repo accompanying the paper [Batch effects in population genomic studies with low-coverage whole genome sequencing data: causes, detection, and mitigation]( https://doi.org/10.1111/1755-0998.13559) by R. Nicolas Lou and Nina Overgaard Therkildsen. 

In this repo, we provide our entire bioinformatic pipeline for data processing and analysis used in our paper. In addition, we put together a tutorial based on a small subset of the data used in our paper, which give researchers an opportunity to go through some of the key bioinformatic methods that we proposed to detect and mitigate batch effects in lcWGS data. We also include some instructions on how the entirety of our raw data can be downloaded from the NCBI SRA.

<br>

**Table of Content**: 

<br>

* [A tutorial with subsetted data and hands-on exercises](https://github.com/therkildsen-lab/batch-effect/blob/main/tutorial/tutorial.md)

<br>

* Bioinformatic pipeline for the paper

    * [The original "batch-effect-naive" pipeline](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/original_pipeline.md)  

    * [Presence/absence of poly-G tails](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/polyg.md)  

    * [Difference in levels of base quality score miscalibration](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/base_quality.md)
    
        * [Base quality score recalibration](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/bqsr.md)

    * [Difference in levels of reference bias](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/reference_bias.md)

    * [Difference in levels of DNA degradation](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/degradation.md)

    * [Difference in sequencing depth](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/depth.md)

    * [Compilation of some multi-panel figures](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/figures.md)

<br>

* [Instructions for downloading the raw fastq files from the NCBI SRA](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/data_download.md)

<br>

* Miscellaneous scripts

  * [Uploading the raw data to the NCBI SRA](https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/data_upload.md)

  * [Data preparation for the tutorial](https://github.com/therkildsen-lab/batch-effect/blob/main/tutorial/data_prep.md)
