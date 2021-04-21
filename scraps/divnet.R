library(DivNet)
library(phyloseq)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
meta <- read.csv("metadata3.csv", row.names=1, header=T)
#ASV count table to phyloseq table
count_phy <- otu_table(ASV_count, taxa_are_rows=T)
sample_info <- sample_data(meta)
viral_physeq <- phyloseq(count_phy, sample_info)
#upload viral tree
virTree<-read_tree("viral_tree")
#add tree to phyloseq object
viral_physeq <- phyloseq(count_phy, sample_info, virTree)


viral_physeq %>%
  sample_richness %>%
  plot

set.seed(20200318)
divnet_phylum <- viral_physeq %>%
  divnet(ncores = 4, tuning = "careful")

divnet_phylum %>% names

divnet_phylum$shannon %>% head


divnet_phylum_char <- lee_phylum %>%
  divnet(X = "char", ncores = 4)

library(ggplot2)
divnet_phylum_char$shannon %>% 
  plot(lee_phylum, color = "char") +
  xlab("Basalt characteristic") +
  ylab("Shannon diversity estimate\n(phylum level)") +
  coord_cartesian(ylim = c(0,2))

simplifyBeta(divnet_phylum_char, lee_phylum, "bray-curtis", "char")

# You can plot this easily
simplifyBeta(divnet_phylum_char, lee_phylum, "bray-curtis", "char") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est,
             col = interaction(Covar1, Covar2))) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

merge_samples(lee_phylum, "char") %>%
  sample_shannon

plugin <- tax_glom(Lee, taxrank="Phylum") %>%
  estimate_richness(measures = "Shannon") %$% Shannon
char <- Lee %>% sample_data %$% char
t.test(plugin[char == "altered"],
       plugin[char == "glassy"])


estimates <- divnet_phylum_char$shannon %>% summary %$% estimate
ses <- sqrt(divnet_phylum_char$`shannon-variance`)
X <- breakaway::make_design_matrix(lee_phylum, "char")
betta(estimates, ses, X)$table

betta(estimates, ses, X)$table["predictorsglassy",3]
round(t.test(plugin[char == "altered"], plugin[char == "glassy"])$p.value, 2)
round(betta(estimates, ses, X)$table["predictorsglassy", 1], 3)
round(mean(plugin[char == "glassy"])-mean(plugin[char == "altered"]), 3)

betta(estimates, ses, X)$global[2]



library(magrittr)
library(phyloseq)
library(breakaway)
library(DivNet) # last checked with version 0.3.5
sessionInfo()

data(Lee)
Lee

tax_glom(Lee, taxrank = "Phylum") %>%
  sample_richness %>%
  plot

lee_phylum <- tax_glom(Lee, taxrank = "Phylum")

set.seed(20200318)
divnet_phylum <- lee_phylum %>%
  divnet(ncores = 4, tuning = "careful")

divnet_phylum %>% names

divnet_phylum$shannon %>% head


divnet_phylum_char <- lee_phylum %>%
  divnet(X = "char", ncores = 4)

library(ggplot2)
divnet_phylum_char$shannon %>% 
  plot(lee_phylum, color = "char") +
  xlab("Basalt characteristic") +
  ylab("Shannon diversity estimate\n(phylum level)") +
  coord_cartesian(ylim = c(0,2))

simplifyBeta(divnet_phylum_char, lee_phylum, "bray-curtis", "char")

# You can plot this easily
simplifyBeta(divnet_phylum_char, lee_phylum, "bray-curtis", "char") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est,
             col = interaction(Covar1, Covar2))) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

merge_samples(lee_phylum, "char") %>%
  sample_shannon

plugin <- tax_glom(Lee, taxrank="Phylum") %>%
  estimate_richness(measures = "Shannon") %$% Shannon
char <- Lee %>% sample_data %$% char
t.test(plugin[char == "altered"],
       plugin[char == "glassy"])


estimates <- divnet_phylum_char$shannon %>% summary %$% estimate
ses <- sqrt(divnet_phylum_char$`shannon-variance`)
X <- breakaway::make_design_matrix(lee_phylum, "char")
betta(estimates, ses, X)$table

betta(estimates, ses, X)$table["predictorsglassy",3]
round(t.test(plugin[char == "altered"], plugin[char == "glassy"])$p.value, 2)
round(betta(estimates, ses, X)$table["predictorsglassy", 1], 3)
round(mean(plugin[char == "glassy"])-mean(plugin[char == "altered"]), 3)

betta(estimates, ses, X)$global[2]
