Data upload
================

``` r
library(tidyverse)
library(googlesheets4)
```

## Put together the biosample attributes

``` r
metadata <- read_sheet(sheet="collection_data", ss = "https://docs.google.com/spreadsheets/d/1I3Pj5VUWfP2eUwwZj1ggAjqgbRmrT0MrLtbQ1lOGaPE/")

sample_table <- read_tsv("../sample_lists/sample_table_merged.tsv") %>%
  dplyr::select(-population_new) %>%
  left_join(metadata, by=c("sample_id_corrected"="sample_id", "population"="population")) %>%
  mutate(date=str_c(day, month, year, sep="-"), 
         date=parse_date(date, "%d-%m-%Y")) %>%
  unnest(lat, keep_empty = TRUE) %>%
  unnest(long, keep_empty = TRUE)


sample_table_biosample <- sample_table %>%
  transmute(`*sample_name` = sample_id_corrected, 
            `*organism` = "Gadus morhua", 
            strain = "not applicable", 
            age = age, 
            `*sex` = case_when(sex=="F"~"female",
                               sex=="M"~"male"), 
            tissue = "fin clips", 
            collection_date = date ,
            geo_loc_name = str_c(country, `ShipPlace/Area`, sep=": "), 
            lat_lon = str_c(round(lat,2), " N ", round(-long,2), " W"),
            store_cond = conservation_method,
            length = length, 
            weight = weight, 
            dressed_weight = dressed_weight, 
            liver_weight = liver_weight,
            gonad_weight = gonad_weight, 
            maturity = maturity,
            ship = ship, 
            trip = trip, 
            station = station)

write_tsv(sample_table_biosample, "../sample_lists/biosample_attributes.tsv", na = "not collected")
```

## Put together the SRA metadata sheet

``` r
metadata <- read_sheet(sheet="library_prep", ss = "https://docs.google.com/spreadsheets/d/1-Koz1X_kIO5JK4hG6uTganu-OvhffijPZma4749BbHg/edit?usp=drive_web&ouid=114221884752263011810") %>%
  filter(library_sequenced==TRUE)

sample_table <- read_tsv("../sample_lists/sample_table_unmerged.tsv") %>%
  left_join(metadata, by="sample_id")

sra_metadata <- sample_table %>%
  transmute(sample_name=sample_id,
            library_ID=str_c(sample_id, seq_id, lane_number, data_type, sep = "_"), 
            title="low-coverage WGS of Gadus morhua: Greenland", 
            library_strategy="WGS",
            library_source="GENOMIC",
            library_selection="RANDOM",
            library_layout=ifelse(data_type=="se", "single", "paired"),
            platform="ILLUMINA",
            instrument_model=ifelse(data_type=="se", "Illumina HiSeq 2500", "NextSeq 500"),
            design_description="DNA was extracted from fin clips or gill tissue with the Qiagen DNeasy Blood & Tissue Kit and libraries were prepared with the protocol described in Therkildsen & Palumbi (2017).",
            filetype="fastq",
            filename=ifelse(data_type=="se", str_c(prefix, ".txt.gz"), str_c(prefix, "_R1.fastq.gz")),
            filename2=ifelse(data_type=="se", "", str_c(prefix, "_R2.fastq.gz")))
write_tsv(sra_metadata, "../sample_lists/sra_metadata.tsv")
```

## Script to transfer files

``` r
sample_table <- read_tsv("../sample_lists/sample_table_unmerged.tsv") 

put <- NULL
for (i in seq_len(nrow(sample_table))){
  if(sample_table$data_type[i]=="pe"){
    put <- c(put, 
              str_c("put ", sample_table$prefix[i], "_R1.fastq.gz"),
              str_c("put ", sample_table$prefix[i], "_R2.fastq.gz"))
  } else {
    put <- c(put, 
             str_c("put ", sample_table$prefix[i], ".txt.gz"))
  }
}
write_lines(put, "../sample_lists/sample_upload.txt")
```
