library(breakaway)
packageVersion("breakaway")
library(phyloseq)
packageVersion("phyloseq")
library(DivNet)
library(tidyverse)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
meta <- read.csv("metadata3.csv", row.names=1, header=T)

#ensure col names of ASV table is same as row names of covariate info
head(colnames(ASV_count) == rownames(meta))

freq_tab <- make_frequency_count_table(ASV_count)
freq_tab %>% head(10)
freq_tab %>% tail(10)
#second col indicates the number of species observed j time, where first col contains j
#ex. 10 species observed 1223

#look to see how many species were observed in the dataset
freq_tab %>% sample_richness

est_richness <- breakaway(freq_tab)
est_richness
