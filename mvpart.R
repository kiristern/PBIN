#https://bookdown.org/forestgeoguest/mpart/mvpart.html

library(mvpart)
library(tidyverse)
library(dplyr)
library(vegan)

# meta_table <- meta#[complete.cases(meta),]
abund_table <- read.table("data/ASVs_counts_copy.tsv", header = T, row.names = 1, check.names = F)
enviro_var <- read.csv("data/metadata3.csv", row.names=1, header=T)
#change to data format
enviro_var[1] <- as.Date(enviro_var$Date)

summary(enviro_var)
env_var <- enviro_var %>% select(1:6, 12:13)
env_var <- env_var[complete.cases(env_var), ]
env_var %>% dplyr::glimpse()
summary(env_var)

abundance <- t(data.matrix(abund_table))

abund_name <- row.names(abundance)
env_row_name <- row.names(env_var)

#check which rows are not the same
abund_name %in% env_row_name
#which rows are the same
#intersect(abund_name, env_row_name)
#specific samples that are not the same
row_remove <- setdiff(abund_name, env_row_name)
#count how many are different
#length(setdiff(abund_name, env_row_name))

#remove rows (samples) that aren't in env_var from abundance
abundance <- abundance[!(row.names(abundance) %in% row_remove), ]

#transform repsonse variables
#The transformation consists of expressing each fish density as a proportion of the sum of all densities 
#in the analytical unit and taking the square root of the resulting value (Legendre and Gallagher 2001).
#The square-root portion of the transformation decreases the importance of the most abundant species.
abun_norm <- decostand(abundance, "hellinger")

#check where there is NAs
env_var %>% 
  filter_all(any_vars(is.na(.))) 

# mvpart_formula <- abundance ~ Date + Months + Years +  Site + Period + bloom2 + 
#   Cumulative_precipitation_t1_t7_mm + Mean_temperature_t0_t7 

mvpart_formula <- abun_norm ~ Period + Site + bloom2 + Months + Years + Cumulative_precipitation_t1_t7_mm +
  Mean_temperature_t0_t7 

tree <- mvpart(mvpart_formula, env_var,
               legend=FALSE, margin=0.01, cp=0, xv="pick",
               xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

str(tree)
summary(tree)

mvpart_run1 <- mvpart(
  form = mvpart_formula, 
  data = env_var,
  all.leaves = TRUE,  # annotate all nodes
  rsq = TRUE,  # give "rsq" plot
  pca = TRUE,  # plot PCA of group means and add species and site information
  wgt.ave.pca = TRUE  # plot weighted averages acorss sites for species
)

rpart.pca(tree, interact = T, wgt.ave = T)
