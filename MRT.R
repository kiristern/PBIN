#https://bookdown.org/forestgeoguest/mpart/mvpart.html

library(mvpart)
library(tidyverse)
library(dplyr)
library(vegan)
library(stats)

setwd("~/Documents/GitHub/PBIN")

# meta_table <- meta#[complete.cases(meta),]
abund_table <- read.table("data/ASVs_counts_copy.tsv", header = T, row.names = 1, check.names = F)
enviro_var <- read.csv("data/metadata3.csv", row.names=1, header=T)
#change to data format
enviro_var[1] <- as.Date(enviro_var$Date)

enviro_var <- enviro_var %>% 
  rename(
    Cumul_precip = Cumulative_precipitation_t1_t7_mm,
    Avg_temp = Mean_temperature_t0_t7
  )

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

# #check where there is NAs
# env_var %>% 
#   filter_all(any_vars(is.na(.))) 

mvpart_formula <- abun_norm ~ Months + Years + Period + Site + bloom2 + Cumul_precip + Avg_temp

tree <- mvpart(abun_norm ~ Period + Months + Years, env_var,
               legend=T, cp=0, xv="pick",
               xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)
plot(tree)
text(tree)

tree_month <- mvpart(abun_norm ~ Months, env_var,
               legend=T, cp=0, xv="pick",
               xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

(month_R2 <- RsquareAdj(tree_month)$adj.r.squared)

tree_year <- mvpart(abun_norm ~ Years, env_var,
                     legend=T, cp=0, xv="pick",
                     xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

tree_period <- mvpart(abun_norm ~ Period, env_var,
                     legend=T, cp=0, xv="pick",
                     xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

tree_site <- mvpart(abun_norm ~ Site, env_var,
                     legend=T, cp=0, xv="pick",
                     xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

tree_bloom <- mvpart(abun_norm ~ bloom2, env_var,
                     legend=T, cp=0, xv="pick",
                     xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

tree_precip <- mvpart(abun_norm ~ Cumul_precip, env_var,
                     legend=T, cp=0, xv="pick",
                     xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

tree_temp <- mvpart(abun_norm ~ Avg_temp, env_var,
                     legend=T, cp=0, xv="pick",
                     xval=nrow(abundance), xvmult=100, which=4, big.pts=T, bars=F)

#view tree details
printcp(tree)
str(tree)
summary(tree)

printcp(tree_month)

# obtain the path to the leaf nodes
tree$frame

leafnodeRows <- grepl("leaf", tree$frame$var)
nodevals <- as.numeric(rownames(tree$frame)[leafnodeRows])
rules <- path.rpart(tree, nodevals)

rulesdf <- do.call(
  "rbind", 
  lapply(rules, function(x) paste(x, collapse = " -AND- "))
)
rulesdf <- data.frame(
  nodeNumber = rownames(rulesdf), 
  rule = rulesdf[, 1], 
  stringsAsFactors = FALSE
)
rulesdf

#PCA
rpart.pca(tree, interact=T, wgt.ave = T, add.tree=T, speclabs = F)

