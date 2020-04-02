#https://bookdown.org/forestgeoguest/mpart/mvpart.html

library(mvpart)
library(tidyverse)
library(dplyr)
library(vegan)
library(stats)

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

#Remove rows with too many NAs
summary(enviro_var)

#remove date col
env_keep <- enviro_var[,-1]

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


mvpart_formula <- abundance_removed ~ Months + Years + Site + Period + bloom2 + Tot_P + Tot_N + Dissolved_P + 
                              Dissolved_N + Cumul_precip + Avg_temp

mrt <- mvpart(as.matrix(abundance_removed) ~., complete_env_keep,  legend=FALSE, margin=0.01, cp=0, xv="pick",
              xval=nrow(abundance_removed), xvmult=100, which=4)
#The graph shows the relative error RE (in green) and the cross-validated relative error CVRE (in blue) of trees of increasing size. 
#The red dot indicates the solution with the smallest CVRE, and the orange dot shows the smallest tree within one standard error of CVRE. 
#It has been suggested that instead of choosing the solution minimizing CVRE, it would be more parsimonious to opt for the smallest tree 
#for which the CVRE is within one standard error of the tree with the lowest CVRE (Breiman et al. 1984). 
#The green bars at the top indicate the number of times each size was chosen during the cross-validation process.

#The residual error (the reciprocal of the R2 of the model, in this case 5.9%%), 
#the cross-validated error, and the standard error. 
#This tree has only two leaves separated by one node. 
#This node splits the data into two groups at the threshold Dissolved_N value of -0.1542.
#Each leaf is characterized by a small barplot showing the abundances of the species, its number of 
#sites and its relative error.

#compare trees
# Using the CVRE criterion (10-group solution)
mrt.cvre <- mvpart(as.matrix(abundance_removed)~., complete_env_keep, 
                         legend=FALSE, margin=0.01, cp=0,xv="pick", 
                         xval=nrow(abundance_removed), xvmult=100,which=4)

# Choosing ourself the best number of partitions
mrt.4 <- mvpart(as.matrix(abundance_removed)~., complete_env_keep, 
                      legend=FALSE, margin=0.01, cp=0, xv="pick", 
                      xval=nrow(abundance_removed), xvmult=100,which=4)







tree <- mvpart(abundance_removed ~ Dissolved_N + Cumul_precip + Tot_N + Tot_P, complete_env_keep,
               legend=T, cp=0, xv="pick",
               xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)
plot(tree)
text(tree)


tree_month <- mvpart(as.matrix(abundance_removed)~ Months, complete_env_keep,
               legend=T, margin=0.01, cp=0, xv="pick",
               xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

(month_R2 <- RsquareAdj(tree_month)$adj.r.squared)

tree_year <- mvpart(as.matrix(abundance_removed)~ Years, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_site <- mvpart(as.matrix(abundance_removed)~ Site, complete_env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_period <- mvpart(as.matrix(abundance_removed)~ Period, complete_env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_bloom <- mvpart(as.matrix(abundance_removed)~ bloom2, complete_env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_totP <- mvpart(as.matrix(abundance_removed)~ Tot_P, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_totN <- mvpart(as.matrix(abundance_removed)~ Tot_N, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_DisP <- mvpart(as.matrix(abundance_removed)~ Dissolved_P, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_DisN <- mvpart(as.matrix(abundance_removed)~ Dissolved_N, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_precip <- mvpart(as.matrix(abundance_removed)~ Cumul_precip, complete_env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_temp <- mvpart(as.matrix(abundance_removed)~ Avg_temp, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(abundance_removed), xvmult=100, which=4, big.pts=T, bars=F)

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

