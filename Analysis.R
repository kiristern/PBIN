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

#####Abundance#####
df <- read.csv("FINAL_ASV_TAXA_ABUN.csv")

df %>% group_by(samples, ASV_ID) %>% summarise(Relative_Abundance = sum(abundance)) %>%
  ggplot(aes(x = samples, y = Relative_Abundance, fill = ASV_ID)) + 
  geom_bar(stat = "identity", show.legend = FALSE)

df_unk <- read.csv("rel_ab_top20_unknown_FINAL.csv")
df_unk %>% group_by(samples, ASV_ID) %>% summarise(Relative_Abundance = sum(abundance)) %>%
  ggplot(aes(x = samples, y = Relative_Abundance, fill = ASV_ID)) + 
  geom_bar(stat = "identity", show.legend = TRUE)

##ordination
##https://joey711.github.io/phyloseq/plot_ordination-examples.html

#just plotting asv
GP.ord <- ordinate(viral_physeq, "NMDS", "bray")
p1 = plot_ordination(viral_physeq, GP.ord)
print(p1)


##using breakaway to explore alpha diversity
##https://github.com/adw96/stamps2018/blob/master/estimation/diversity-lab.R

viral_physeq
viral_physeq %>% sample_data

#how to extract only values corresponding to period "spring" "summer" "fall"
# season <- viral_physeq %>% 
#   subset_samples(Period %in% c("Spring", "Summer", "Fall"))

#look at observed richness plot colour by month
observed <- sample_richness(viral_physeq)
summary(observed)
plot(observed, viral_physeq, color="Months")

#observed richness 
#Sequencing depth by Year
# data.frame("observed_richness" = (observed %>% summary)$estimate,
#            "Depth" = phyloseq::sample_sums(viral_physeq), # sequence depth
#            "type" = viral_physeq %>% sample_data %>% get_variable("Years")) %>%
#   ggplot(aes(x = Depth, y = observed_richness, color = type)) +
#   geom_point()
# 
# #Year by month
# data.frame("observed_richness" = (observed %>% summary)$estimate,
#            "Years" = viral_physeq %>% sample_data %>% get_variable("Years"),
#            "type" = viral_physeq %>% sample_data %>% get_variable("Months")) %>%
#   ggplot(aes(x = Years, y = observed_richness, color = type)) +
#   geom_point()
# 
# #Month by year
# data.frame("observed_richness" = (observed %>% summary)$estimate,
#            "Months" = viral_physeq %>% sample_data %>% get_variable("Months"),
#            "type" = viral_physeq %>% sample_data %>% get_variable("Years")) %>%
#   ggplot(aes(x = Months, y = observed_richness, color = type)) +
#   geom_point()

#estimate number of missing species using a species richness estimate
ba <- breakaway(viral_physeq)
ba
plot(ba, viral_physeq, color="Years")

#get sample names
# x <- viral_physeq %>% sample_data %>% rownames()

#table of sample names, richness, and error 
# alpha_err <- data.frame(x,
                        # summary(ba)$estimate,
                        # summary(ba)$error,
                        # make_design_matrix(viral_physeq, "Years"))

# richness <- summary(ba)$estimate
# std_dev <- summary(ba)$error
# 
# plot(x, richness,
#      ylim = range(c(richness-std_dev, richness+std_dev)),
#      xlim = 
#      pch=19, xlab = "Sample", ylab = "Richness +/- Std Dev",
#      main = "Alpha div with std.dev error bars"
# )
# #hack: draw arrow heads as lines
# arrows(x, richness-std_dev, x, richness+std_dev, length=0.05, angle = 90, code=3)

#for variables Years and bloom2
ba_alpha = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
           "Year" = viral_physeq %>% sample_data %>% get_variable("Years"))

ba_plot <-  ggplot(ba_alpha, aes(x = Bloom, y = ba_observed_richness))+
  geom_point()

#geom_crossbar()
ba_plot + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                       geom="crossbar", width=0.5)

#geom_errorbar()
ba_plot + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red")


#look at just one sample
tr <- viral_physeq %>% subset_samples(Date == "2006-5-30")
tr
#look at the structure of this dataset
freq_count <- tr %>% otu_table %>% make_frequency_count_table
#this is the frequency count table for this dataset
freq_count %>% head(10) #14 ASVs observed only twice  
freq_count %>% tail() #1 ASV observed 293 times

#fit breakaway to the sample
ba_tr <- breakaway(freq_count)
ba_tr #this is an alpha diversity estimate -- a special class for alpha div. estimates

#check which model breakaway picked
ba_tr$model

#take estimates and turn into df 
# summary(ba) %>%
#   add_column("SampleNames" = viral_physeq %>% otu_table %>% sample_names)

#chose a different species richness estimate
chao_est <- viral_physeq %>%
  chao1 %>%
  plot(viral_physeq, color="Months")



#betta() works like a regression model but accounts for the uncertainty in estimating diversity
# test hypothesis that different types of water systems have the same microbial div (not in this case)
bt <- betta(summary(ba)$estimate,
            summary(ba)$error,
            make_design_matrix(viral_physeq, "Site"))
bt$table
#betta() estimates that the mean ASV-level diversity in ... is 368?
#estimates the the diversity in March is significantly lower (~222) than all others but that they're all below 0??
#betta accounts for the error bars in diversity when doing hypothesis testing

#Shannon index: downweights the importance of rare taxa (reflects low abundance of a taxon)
#DivNet adjusts for different sequencing depths so don't need to rarefy (throw away data)

#Run in parallel
dv_viral_ps <- divnet(viral_physeq, ncores = 4)



