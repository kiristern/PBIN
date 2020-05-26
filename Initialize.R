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
library(reshape2)
library(genefilter)
library(TSA)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
meta_cyano <- read.csv("metadata_w_cmd.csv", row.names = 1, header = T)



#meta$Years <- as.factor(meta$Years)
meta_cyano$Years <- as.factor(meta_cyano$Years)


#ASV count table to phyloseq table
count_phy <- otu_table(ASV_count, taxa_are_rows=T)
sample_info <- sample_data(meta_cyano)
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
write.table(filt_vir, "filt_vir.txt", row.names = TRUE, quote=FALSE)

#filtering based on the across-sample mean; hence, if an OTU is abundant in some samples, but very rare in the rest of the samples, then it will be categorized as rare
mean_rare <- filter_taxa(viral_physeq, function(x) mean(x) < 300, TRUE)

##filter based on per-sample (not across-sample) abundance 
#define criterion for filtering; kOverA evaluates to TRUE if at least k (1) of the arguments are >A (34)
flist_a <- filterfun(kOverA(1, 300))
a <- filter_taxa(viral_physeq, flist_a)
rare_ps <- prune_taxa(!a, viral_physeq)
raredf <- as.data.frame(rare_ps %>% otu_table())


flist    <- filterfun(kOverA(1, 300))
ent.logi <- filter_taxa(viral_physeq, flist)
ent.trim <- filter_taxa(viral_physeq, flist, TRUE)
identical(ent.trim, prune_taxa(ent.logi, viral_physeq)) 
identical(sum(ent.logi), ntaxa(ent.trim))
filter_taxa(viral_physeq, flist, TRUE)


######## prep data for MRT and RDA ########
#transform asv density as a proportion of the sum of all densities
vir_abund_helli <-decostand(t(filt_vir), method="hellinger")



#### read in cyano with env var table ####
env_cy <- meta_cyano
colnames(env_cy)

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
                              "Cyano.Abundance",
                              "Micro.Abundance",
                              "Dolicho.Abundance")
colnames(env_keep)
#standardize environmental data
env_keep[,c(6:11)]<-decostand(env_keep[,c(6:11)], method="standardize")
#env_keep[,c(6:14)]<-decostand(env_keep[,c(6:14)], method="standardize")

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









