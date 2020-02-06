library(breakaway)
library(DivNet)
library(tidyverse)
library(phyloseq)
library(qiime2R)

setwd("~/Documents/GitHub/PBIN/data")
list.files()

count_tab <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)

# Clean data: filter ASV less than 10 occurences
# create new col with ASV sum
new_count_tab <- mutate(count_tab, rowSums(count_tab))
# filter only ASVs with more than 10
new_count_tab <- filter(new_count_tab, rowSums(count_tab)>100)
# remove rowSums col
new_count_tab <- select(new_count_tab, -"rowSums(count_tab)")

#upload metadata table that was manually merged with asv table
# meta_asv <- read.csv("meta_asv.csv")
meta <- read.csv("metadata3.csv", row.names = 1, header=T)

#reorder cols
col_order <- c("SampleID", "Date", "Months", "Years", "Period", "Site", "bloom2", "Total_Phosphorus_ug",
               "Phosph_Range", "Nitrog_Range", "Total_Nitrogen_mg", "Temperature_Water_Celsius", 
               "Dissolved_P" , "Dissolved_N", "Cumulative_precipitation_t1_t7_mm", "Profondeur_Secchi_cm" ,
               "Mean_temperature_t0_t7"  , "Microcystin_ug_L", "Description" , "echantillonnage")
meta <- meta[, col_order]

#ASV count table to phyloseq table
count_phy <- otu_table(new_count_tab, taxa_are_rows=T)
sample_info <- sample_data(meta)
viral_physeq <- phyloseq(count_phy, sample_info)

virTree<-read_tree("viral_tree")
#add tree to phyloseq
viral_physeq <- phyloseq(count_phy, sample_info, virTree)


```
count to phyloseq table
calculate Shannon for all samples, if it's below, then that's the cutoff for bloom no bloom
```
