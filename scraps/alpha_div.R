library(breakaway)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
meta <- read.csv("metadata3.csv", row.names=1, header=T)

#ensure col names of ASV table is same as row names of covariate info
head(colnames(ASV_count) == rownames(meta))

#create frequency table
freq_tab_list <- build_frequency_count_tables(ASV_count)


