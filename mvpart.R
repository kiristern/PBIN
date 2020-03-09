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

#RDA
#Standardize the environmental variables table: the environmental data could be in different units, then a standardization is required to perform RDA
# env_var2 <- env_var
# std_env_var <- decostand(env_var2, method="standardize")

# rda_formula <- abun_norm ~ Period + Site + bloom2 + Months + Years + Cumul_precip + Avg_temp + Date
rda_formula <- function(rda_variable){
  rda_var <- rda(abun_norm ~ rda_variable, data=env_var)
  score_time_var <- scores(rda_var, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
  summary(rda_var, display=NULL)
}

rda_formula(Years)

rda_time_var <- rda(abun_norm ~ Period + Months + Years, data=env_var)
score_time_var <- scores(rda_time_var, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
summary(rda_time_var, display=NULL)

rda_period <- rda(abun_norm ~ Period, data=env_var)
score_period <- scores(rda_period, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
summary(rda_period, display=NULL)

rda_month <- rda(abun_norm ~ Months, data=env_var)
score_month <- scores(rda_month, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
summary(rda_month, display=NULL)

rda_year <- rda(abun_norm ~ Years, data=env_var)
score_year <- scores(rda_year, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
summary(rda_year, display=NULL)

rda_precip <- rda(abun_norm ~ Cumul_precip, data=env_var)
score_precip <- scores(rda_precip, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
summary(rda_precip, display=NULL)

rda_temp <- rda(abun_norm ~ Avg_temp, data=env_var)
score_temp <- scores(rda_temp, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
summary(rda_temp, display=NULL)

rda_site <- rda(abun_norm ~ Site, data=env_var)
score_site <- scores(rda_site, display=c("sp", "wa", "lc", "bp", "cn"), scaling=2)
summary(rda_site, display=NULL)

#adjusted R squared
(R2adj <- RsquareAdj(rda_period)$adj.r.squared) #the strength of the relationship between X and Y corrected for the number of X variables is
RsquareAdj(rda_year)$adj.r.squared
RsquareAdj(rda_month)$adj.r.squared
RsquareAdj(rda_precip)$adj.r.squared
RsquareAdj(rda_temp)$adj.r.squared
RsquareAdj(rda_site)$adj.r.squared


# ordiR2step(rda(abun_norm~1, data=env_var), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)

#extract site data
df_sites<-data.frame(score_time_var$sites)
colnames(df_sites)<-c("x", "y")
#plot sites
plot_time_var <-ggplot() %>%
  + geom_point(data=df_sites,aes(x,y)) %>%
  + labs(x="RDA1(x%)", y="RDA2(x%)") #x values could be taken from summary(rda_tab)
#add arrows
multiplier <- vegan:::ordiArrowMul(score_time_var$biplot)
df_arrows<- score$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)
plot_time_var <-plot_time_var %>% 
  + geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.6)
plot_time_var
#add labels
plot_time_var2 <- plot_time_var %>%
  + geom_text(data=as.data.frame(df_arrows*1.3),aes(x, y, label = rownames(df_arrows),size=6,fontface="bold"),
              color="black",alpha=0.6, size=3)
plot_time_var2
# Draw species
df_scr_sp<- as.data.frame(score$species)
colnames(df_scr_sp)<-c("x","y")
df_scr_sp



# Draw species
df_scr_sp<- as.data.frame(score$species)
colnames(df_scr_sp)<-c("x","y")
df_scr_sp

