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

#create frequency table
freq_tab_list <- build_frequency_count_tables(ASV_count)
head(freq_tab_list)

#estimating species richness
breakaway(freq_tab_list[[1]]) #run on first freq count tab
  #no plot is provided, so returned model WLRM
  #possible reasons: too many singletons, not a long enough tail, false diversity, not enough data
#test if failure was sensitive to singleton count
breakaway_nof1(freq_tab_list[[1]][-1])

#objective bayes
  #get distribution of estimates
set.seed(1234)
objective_bayes_negbin(freq_tab_list[[1]], plot=F)
objective_bayes_poisson(freq_tab_list[[1]])$results
objective_bayes_geometric(freq_tab_list[[1]])$results
objective_bayes_mixedgeo(freq_tab_list[[1]])$results



