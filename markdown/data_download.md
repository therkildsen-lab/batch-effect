Instructions for downloading raw data from the NCBI SRA
================

Our raw sequencing data have been uploaded to NCBI SRA (Accession Number
SRR16894036-SRR16894202). We encourage other researcher to obtain this
publically available data either for reproducing our results or for
their own research. Here, we provide an instruction on how our raw
sequencing data files can be downloaded from the SRA.

First, please first clone this GitHub repo to your Linux server. For
example,

``` bash
git clone https://github.com/therkildsen-lab/batch-effect.git
```

Then, define the location of this GitHub repo as `BASEDIR`, and create a
subdirectory named `raw_fastq` within `BASEDIR`. For example,

``` bash
BASEDIR=/workdir/batch-effect
mkdir $BASEDIR/raw_fastq
```

Lastly, run the following loop to download the raw fastq files from the
NCBI SRA. You can also save it as a bash script and use `nohup` to run
it in the background. Gzipped fastq files will be downloaded to the
`raw_fastq` folder.

``` bash
for LINE in `cat $BASEDIR/sample_lists/SRR_Acc_List.txt`; do
  fasterq-dump --outdir $BASEDIR/raw_fastq --split-files $LINE
  gzip ${LINE}*.fastq
done
```

Metadata for these fastq files are stored in the `SraRunTable.txt` file
within the `sample_lists` folder in this GitHub repo. This metadata
table (as well as the accession list used in the previous step) was
downloaded from the SRA Run Selector:
<https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP345309&o=acc_s%3Aa>

Since SRA has altered the original fastq file names, to rerun our data
analysis, use `sample_table_unmerged_sra.tsv`, `fastq_list_sra.txt`,
`fastq_list_sra_pe.txt`, `fastq_list_sra_se.txt`, instead of the
original sample tables and fastq list until the merging step
(i.e. adapter trimming, quality filtering, read alignment and sorting).
The following script is how these new sample table and lists are
generated (you don’t need to run them).

``` r
library(tidyverse)
metadata <- read_csv("../sample_lists/SraRunTable.txt") %>%
  dplyr::select(Run, `Library Name`)
sample_table_old <- read_tsv("../sample_lists/sample_table_unmerged.tsv") %>%
  mutate(`Library Name`=str_c(sample_id, seq_id, lane_number, data_type, sep = "_"))
sample_table_new <- left_join(sample_table_old, metadata) %>%
  mutate(prefix=`Run`) %>%
  dplyr::select(-`Library Name`, -`Run`)
write_tsv(sample_table_new, "../sample_lists/sample_table_unmerged_sra.tsv")
filter(sample_table_new, data_type=="pe") %>%
  .$prefix %>%
  write_lines("../sample_lists/fastq_list_pe_sra.txt")
filter(sample_table_new, data_type=="se") %>%
  .$prefix %>%
  write_lines("../sample_lists/fastq_list_se_sra.txt")
sample_table_new$prefix %>%
  write_lines("../sample_lists/fastq_list_sra.txt")
```
