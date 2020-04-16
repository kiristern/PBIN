library(breakaway)
packageVersion("breakaway")
library(phyloseq)
packageVersion("phyloseq")
library(DivNet)
library(tidyverse)
library(qiime2R)
library(mvpart)
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyverse)
library(stats)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
meta <- read.csv("metadata3.csv", row.names=1, header=T)

meta$Years <- as.factor(meta$Years)

#ASV count table to phyloseq table
count_phy <- otu_table(ASV_count, taxa_are_rows=T)
sample_info <- sample_data(meta)
viral_physeq <- phyloseq(count_phy, sample_info)

#upload viral tree
virTree<-read_tree("viral_tree")

#add tree to phyloseq object
viral_physeq <- phyloseq(count_phy, sample_info, virTree)

#view ASV count (top 10)
sort(taxa_sums(viral_physeq), decreasing = T)[1:10]

#check data
print(viral_physeq)




######## prep data for MRT and RDA ########

setwd("~/Documents/GitHub/PBIN")

#load in data
abund_table <- read.table("data/ASVs_counts_copy.tsv", header = T, row.names = 1, check.names = F)
#transform asv density as a proportion of the sum of all densities
abund_table3 <-decostand(abund_table, method="hellinger")

#load metadata
enviro_var <- read.csv("data/metadata3.csv", row.names=1, header=T)
#standardize environmental data
enviro_var[,c(7:12)]<-decostand(enviro_var[,c(7:12)], method="standardize")

#rename cols
enviro_var <- enviro_var %>%
  rename(
    Cumul_precip = "Cumulative_precipitation_t1_t7_mm",
    Avg_temp = "Mean_temperature_t0_t7",
    Tot_P = "Total_Phosphorus_ug",
    Tot_N = "Total_Nitrogen_mg"
  )
# when tidyverse decides not to load and don't want to restart R:
# enviro_var <- rename(enviro_var, c("Cumulative_precipitation_t1_t7_mm"="Cumul_precip",
#                                    "Mean_temperature_t0_t7" = "Avg_temp",
#                                    "Total_Nitrogen_mg" = "Tot_N",
#                                    "Total_Phosphorus_ug" = "Tot_P"))

#add cyano data to meta
cyano_var <- read.table("data/cyano/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
cyano_var <- cyano_var %>% select(1:135)
#remove X in front of date
colnames(cyano_var) <- substring(colnames(cyano_var), 2)
cyano_var_sum <- colSums(cyano_var)
cyano_var_sum <- as.data.frame(cyano_var_sum)
# write.csv(cyano_var_sum, "cyano_var_sum.csv")
env_cy <- read.csv("data/metadata_w_cyano.csv")
colnames(env_cy)
#standardize environmental data
env_cy[,c(8:14)]<-decostand(env_cy[,c(8:14)], method="standardize")

#rename cols
env_cy <- env_cy %>%
  rename(
    Cumul_precip = "Cumulative_precipitation_t1_t7_mm",
    Avg_temp = "Mean_temperature_t0_t7",
    Tot_P = "Total_Phosphorus_ug",
    Tot_N = "Total_Nitrogen_mg"
  )

#Remove rows with too many NAs
summary(env_cy)

#remove date col
env_keep <- env_cy[,-c(1,2)]

#check how many rows there are without any NAs
complete.cases(env_keep)

complete_env_keep <- env_keep[complete.cases(env_keep), ]
complete_env_keep %>% dplyr::glimpse()
summary(complete_env_keep)

abundance <- t(data.matrix(abund_table3))

#look at the species' distribution frequencies
(viral_ab <- table(unlist(abundance)))
# barplot(viral_ab, las=1, xlab = "Abundance class", ylab="Frequency")

#see how many absences
sum(abundance==0)
#look at the proportion of zeros in community data
sum(abundance==0)/(nrow(abundance)*ncol(abundance))

#comparing removed env rows with abundance samples
abund_name <- row.names(abundance)
env_row_name <- row.names(complete_env_keep)

#check which rows are not the same
abund_name %in% env_row_name
#which rows are the same
#intersect(abund_name, env_row_name)

#specific samples that are not the same
(row_remove <- setdiff(abund_name, env_row_name))
#count how many are different
length(setdiff(abund_name, env_row_name))

#remove rows (samples) that aren't in env_var from abundance
abundance_removed <- abundance[!(row.names(abundance) %in% row_remove), ]







