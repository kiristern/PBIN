#https://wiki.qcbs.ca/r_workshop10

library(vegan)
library(mvpart)
library(MVPARTwrap)
library(rdaTest)
library(labdsv)
library(plyr)
library(MASS)
library(tidyverse)
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
# when tidyverse decides not to load and don't want to restart R:
# enviro_var <- rename(enviro_var, c("Cumulative_precipitation_t1_t7_mm"="Cumul_precip", "Mean_temperature_t0_t7" = "Avg_temp"))

#Remove rows with too many NAs
summary(enviro_var)
env_var <- enviro_var %>% select(1:6, 12:13)
env_var <- env_var[complete.cases(env_var), ]
env_var %>% dplyr::glimpse()
summary(env_var)

abundance <- t(data.matrix(abund_table))

#look at the species' distribution frequencies
(viral_ab <- table(unlist(abundance)))
barplot(viral_ab, las=1, xlab = "Abundance class", ylab="Frequency")

#see how many absences
sum(abundance==0)

#look at the proportion of zeros in community data
sum(abundance==0)/(nrow(abundance)*ncol(abundance))

#
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

dim(abun_norm)
str(abun_norm)
head(abun_norm)
summary(abun_norm)


##### Environmental data #######
names(env_var)
dim(env_var)
str(env_var)
head(env_var)
summary(env_var)
pairs(env_var, main="Bivariate Plots of the Environmental Data" ) 

# #Standardize the environmental data (11 variables) using the function decostand() in vegan.
# env_std <- decostand(env_var, method = "standardize")
# apply(env_std, 2, mean) #centered data (mean~0)
# apply(env_std, 2, sd) #scaled data (std dev=1)

######## RDA #######

?rda
rda_spe <- rda(abun_norm~., data = env_var)

#extract results
summary(rda_spe, display = NULL)
  #the explanatory variables (env vars) explain 38.02% of the variance in Y (species)

#select the significant explanatory variables
?ordiR2step
ordiR2step(rda(abun_norm~1, data=env_var), scope= formula(rda_spe), direction= "forward", R2scope=TRUE, pstep=1000)
  #the proportion of variation explained by the constraining variables being 0.3542

#retain significant variables only
signif_env <- subset(env_var, select = c("Months", "Date", "Period", "bloom2", "Cumul_precip"))

#rda with significant variables 
rda_spe_signif <- rda(abun_norm~., data=signif_env)
summary(rda_spe_signif, display=NULL)
  #The proportion of the variance of Y (species) explained by the X (env) variables = 35.42% (constrained) 
  #the unexplained variance of Y = 64.58% (unconstrained)

#adjusted R^2
(adjR2 <- RsquareAdj(rda_spe_signif)$adj.r.squared) 
  #strength of the relationship between X and Y corrected for the number of X variables is 27.94%
  #the explanatory variables explain 27.94% of the variance in Y (species)

#test significance of the RDA canonical axis (ie. constraining variables, ie. env var)
?anova.cca
anova.cca(rda_spe_signif, step=1000) 
  #RDA model is highly significant (p=0.001)
anova.cca(rda_spe_signif, step=1000, by="axis") 
  #axis' 1-3 v significant(p=0.001), 4 (p=0.006), 5(p=0.045), 6-11 not

#plot RDA
#Quick plots scaling 1 and 2
windows()
plot(rda_spe_signif, scaling=1, main="Triplot RDA (scaling 1)")
windows()
plot(rda_spe_signif, scaling=2, main="Triplot RDA (scaling 2)")

#Advanced plots scaling 1
  #distance among objects approximate their Euclidean distances
  #used to interpret distances among objects
windows()
plot(rda_spe_signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA3"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(rda_spe_signif, display="sites", choices=c(1,3), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(rda_spe_signif, display="species", choices=c(1), scaling=1),
       scores(rda_spe_signif, display="species", choices=c(3), scaling=1),
       col="black",length=0)
text(scores(rda_spe_signif, display="species", choices=c(1), scaling=1),
     scores(rda_spe_signif, display="species", choices=c(3), scaling=1),
     labels=rownames(scores(rda_spe_signif, display="species", scaling=1)),
     col="black", cex=0.8)    
arrows(0,0,
       scores(rda_spe_signif, display="bp", choices=c(1), scaling=1),
       scores(rda_spe_signif, display="bp", choices=c(3), scaling=1),
       col="red")
text(scores(rda_spe_signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(rda_spe_signif, display="bp", choices=c(3), scaling=1)+0.05,
     labels=rownames(scores(rda_spe_signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#Advanced plots scaling 2
  #angles ebtween var X and Y relect their correlation
  #used to interpret the relationship betwen X and Y
windows()
plot(rda_spe_signif, scaling=2, main="Triplot RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(rda_spe_signif, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(rda_spe_signif, display="species", choices=c(1), scaling=2)*2,
       scores(rda_spe_signif, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(rda_spe_signif, display="species", choices=c(1), scaling=2)*2.1,
     scores(rda_spe_signif, display="species", choices=c(2), scaling=2)*2.1,
     labels=rownames(scores(rda_spe_signif, display="species", scaling=2)),
     col="black", cex=0.8)    
arrows(0,0,
       scores(rda_spe_signif, display="bp", choices=c(1), scaling=2),
       scores(rda_spe_signif, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(rda_spe_signif, display="bp", choices=c(1), scaling=2)+0.05,
     scores(rda_spe_signif, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(rda_spe_signif, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  











