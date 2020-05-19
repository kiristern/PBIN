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
library(psych)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
#meta <- read.csv("metadata3.csv", row.names=1, header=T)
meta <- read.csv("meta_cmd.csv", row.names = 1, header = T)


#meta$Years <- as.factor(meta$Years)
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

#remove taxa not seen more than 3 times in at least 5% of the samples. This protects against ASV with small mean & trivially large coef of var
filt_virseq <- filter_taxa(viral_physeq, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
filt_vir <- as.data.frame(filt_virseq %>% otu_table())

# #standardize abundances to the median sequencing depth
# total <- median(sample_sums(filt_vir_seq))
# standf <- function(x, t=total) round(t * (x / sum(x)))
# st_filt_virseq <- transform_sample_counts(filt_vir_seq, standf)
# #filter the taxa using a cutoff of 3.0 for the coef of var
# virseq <- filter_taxa(st_filt_virseq, function(x) sd(x)/mean(x) > 3.0, TRUE)

#transform to relative abundance BETTER TO USE DECOSTAND METHOD=HELLINGER???
# rel_virseq  <- transform_sample_counts(filt_vir_seq, function(x) x / sum(x) )
#filter such that only ASV with mean greater than 10^5 are kept
#filt_rel_virseq <- filter_taxa(rel_virseq, function(x) mean(x) > 1e-5, TRUE)

#View phyloseq objects
# filt_vir_seq %>% otu_table()
# ntaxa(virseq)
# nsamples(virseq)
# sample_names(virseq)[1:5]
# sample_variables(virseq)
# otu_table(virseq)[1:5, 1:5]




######## prep data for MRT and RDA ########

setwd("~/Documents/GitHub/PBIN")

#transform asv density as a proportion of the sum of all densities
vir_abund_helli <- decostand(t(filt_vir), method="hellinger")

#load metadata
# enviro_var <- meta
# #standardize environmental data
# enviro_var[,c(7:13)]<-decostand(enviro_var[,c(7:13)], method="standardize")
# 
# #rename cols
# enviro_var <- enviro_var %>%
#   rename(
#     Cumul_precip = "Cumulative_precipitation_t1_t7_mm",
#     Avg_temp = "Mean_temperature_t0_t7",
#     Tot_P = "Total_Phosphorus_ug",
#     Tot_N = "Total_Nitrogen_mg"
#   )
# when tidyverse decides not to load and don't want to restart R:
# enviro_var <- rename(enviro_var, c("Cumulative_precipitation_t1_t7_mm"="Cumul_precip",
#                                    "Mean_temperature_t0_t7" = "Avg_temp",
#                                    "Total_Nitrogen_mg" = "Tot_N",
#                                    "Total_Phosphorus_ug" = "Tot_P"))

#### add cyano data to meta #####
# cyano_var <- read.table("data/cyano/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
# select_cyano_var <- cyano_var %>% select(1:135)
# 
# #transform cyano asv density as a proportion of the sum of all densities
# cyano_abund_helli <-decostand(select_cyano_var, method="hellinger")
# 
# #remove X in front of date
# colnames(select_cyano_var) <- substring(colnames(select_cyano_var), 2)
# cyano_var_sum <- colSums(select_cyano_var)
# cyano_var_sum <- as.data.frame(cyano_var_sum)
# # write.csv(cyano_var_sum, "cyano_var_sum.csv")

#### read in cyano with env var table ####
env_cy <- meta
colnames(env_cy)
#standardize environmental data
env_cy[,c(1:3, 10:15)]<-decostand(env_cy[,c(1:3, 10:15)], method="standardize")

#rename cols
env_cy <- env_cy %>%
  rename(
    Cumul_precip = "Cumulative_precipitation_t1_t7_mm",
    Avg_temp = "Mean_temperature_t0_t7",
    Tot_P = "Total_Phosphorus_ug",
    Tot_N = "Total_Nitrogen_mg"
  )

###Remove rows with too many NAs
#look at data
summary(env_cy)

#remove date col
#env_keep <- env_cy[,-1]
#remove certain env variables to reduce number of NAs
env_keep <- env_cy %>% select("Months",
                              "Years",
                              "Site",
                              "Period",
                              "bloom2",
                              "Tot_P",
                              "Tot_N",
                              "Dissolved_P", 
                              "Dissolved_N",
                              "Cumul_precip",
                              "Avg_temp",
                              "Dolicho.Abundance",
                              "Micro.Abundance",
                              "Cyano.Abundance")


summary(env_keep)
# colnames(env_keep)

#check how many rows there are without any NAs
complete.cases(env_keep)

complete_env_keep <- env_keep[complete.cases(env_keep), ]
complete_env_keep %>% dplyr::glimpse() 
summary(complete_env_keep)


#### Remove viral asvs that are not present (due to removal of NA from env vars)
vir_abund_helli
#look at the species' distribution frequencies
#viral_ab <- table(unlist(vir_abundance))
#barplot(viral_ab, las=1, xlab = "Abundance class", ylab="Frequency")

#see how many absences
sum(vir_abund_helli==0)
#look at the proportion of zeros in community data
sum(vir_abund_helli==0)/(nrow(vir_abund_helli)*ncol(vir_abund_helli))

#comparing removed env rows with cyano abundance samples
abund_name <- row.names(vir_abund_helli)
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
vir_abun_removed <- vir_abund_helli[!(row.names(vir_abund_helli) %in% row_remove), ]
nrow(vir_abun_removed)






