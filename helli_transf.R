# Hellinger transformation: RDA, MRT
ASV_count
dim(ASV_count)
cyano_counts
dim(cyano_counts)

#Hellinger transformation
#https://wiki.qcbs.ca/r_workshop9 
tASV_count <- t(ASV_count) #species should be in columns
dim(tASV_count)
vir_hel<-decostand(tASV_count, method="hellinger")

tcyano_counts <- t(cyano_counts)
bact_hel<-decostand(tcyano_counts, method="hellinger")

# OR

bactps_helli <- transform(bact_physeq, transform = "hellinger", target = "OTU")

bactps_helli %>% otu_table() %in% bact_hel #both methods are the same! yay!

cyano_ps_helli <- subset_taxa(bactps_helli, Phylum == "p__Cyanobacteria")
cyano_helli <- cyano_ps_helli %>% otu_table() 

##set vir_hel to same format as bact_hel
#remove sample ID at beginning
row.names(vir_hel) <- sub("*._*._*._*._*._*._*._","", row.names(vir_hel))
#change "_" to "."
rownames(vir_hel) <- gsub("_", ".", row.names(vir_hel))

#CLR transformation
library(zCompositions)
abund.no0 = tASV_count[ ,colSums(tASV_count)!=0 ]
abund.count <- t(cmultRepl(t(abund.no0), method="CZM", output="p-counts"))
abund.prop.rows <- decostand(abund.count,"total",MARGIN=1)
abund.clr <- t(apply(abund.prop.rows, 1, function(x){log(x) - mean(log(x))}))
vir.abund.clr <- as.data.frame(abund.clr)

bact.abund.no0 = tcyano_counts[ ,colSums(tcyano_counts)!=0 ]
bact.abund.count <- t(cmultRepl(t(bact.abund.no0), method="CZM", output="p-counts"))
bact.abund.prop.rows <- decostand(bact.abund.count,"total",MARGIN=1)
bact.abund.clr <- t(apply(bact.abund.prop.rows, 1, function(x){log(x) - mean(log(x))}))
bact.abund.clr <- as.data.frame(bact.abund.clr)



##### RDA TRANSFORMED HELLINGER ######
library(vegan)
transposed_asv_count_meta <- t(asv_count_meta)
meta

vir_helli2 <- decostand(transposed_asv_count_meta, method="hellinger")
env <- meta
names(env)

#remove catagorical data from env (do RDA without sites and time -- see PERMANOVA (far) below)
str(env)
temp_vars <- env[,!(colnames(env) %in% c("Total_Phosphorus_ug", "Total_Nitrogen_mg", "Dissolved_P", "Dissolved_N", "Cumulative_precipitation_t1_t7_mm", "Mean_temperature_t0_t7",
                                         "Microcystin", "Dolicho.Abundance", "Micro.Abundance", "Cyano.Abundance", "cyano.sum.helli",
                                         "micro.sum.helli" , "doli.sum.helli", "description", "Date"))]


env_vars <- env[,!(colnames(env) %in% c("description", "Date", "Months", "Years", "Site", "Period", "bloom2",
                                        "Microcystin", "Dolicho.Abundance", "Micro.Abundance", "Cyano.Abundance"))]
#"cyano.sum.helli", "doli.sum.helli", "micro.sum.helli"))]

colnames(env_vars)

#standardize environmental data
library(vegan)
env_vars.std <- env_vars
env_vars.std[1:6] <-decostand(env_vars.std[1:6], method="standardize")
head(env_vars.std, n=2)
# apply(env_vars[,8:16], 2, mean) #data are now centered (mean~0)
# apply(env_vars[,8:16], 2, sd) #data are now scaled (st. dev =1)

#check how many rows there are without any NAs
complete.cases(env_vars.std)
#rm NAs
env_keep <- env_vars.std[complete.cases(env_vars.std), ]
env_keep %>% dplyr::glimpse() 
summary(env_keep)

temp_keep <- temp_vars[complete.cases(temp_vars), ]

#### Remove viral asvs that are not present (due to removal of NA from env vars)
sp.asv <- t(vir_helli)

#rm sample rows that are not present in env_keep/temp_keep
vir.rm <- sp.asv[rownames(sp.asv) %in% rownames(env_keep),]
dim(vir.rm)
#which rows are the same
#intersect(abund_name, env_row_name)

# #specific samples that are not the same
# (row_remove <- setdiff(row.names(sp.asv), row.names(env_keep)))
# #count how many are different
# length(setdiff(row.names(sp.asv), row.names(env_keep)))  
# 
# #remove rows (samples) that aren't in env_var from abundance
# sp.asv.rm <- sp.asv[!(row.names(sp.asv) %in% row_remove), ]
# dim(sp.asv.rm)

# species and enviro data without NAs needed for RDA
head(vir.rm)
head(env_keep)

length(env_keep) #number of env vars
length(vir.rm) #nnumber asv
nrow(env_keep) #number samples
nrow(vir.rm) #numer of samples

vir.rda <- rda(vir.rm~., data = env_keep)

summary(vir.rda, display = NULL)
#These results contain: (1) the proportion of variance of Y explained by the X variables (constrained proportion, 16.17% here)
#(2) the unexplained variance of Y (unconstrained proportion, 83.83% here) and 
#(3) then summarize the eigenvalues, the proportions explained and the cumulative proportion of each canonical axis 
# (each canonical axis = each constraining variable, in this case, the environmental variables from env).

# Select the significant explanatory variables by forward selection
ordiR2step(rda(vir.rm~1, data = env_keep), scope=formula(vir.rda), direction="forward", R2scope = T, pstep=1000)

anova.cca(vir.rda, by ="terms")

temp_sig <- adonis(vir.rm~., data = temp_keep)
str(temp_keep)
temp_sig
env_sig <- adonis(vir.rm~., data = env_keep) 
env_sig

#get adjusted R2
(R2adj <- RsquareAdj(vir.rda)$adj.r.squared)

#test significance of the model and of each canonical axis
anova.cca(vir.rda, step=1000) #model
anova.cca(vir.rda, step=1000, by="axis") #canonical axes

#significant env vars
env.signif <- subset(env_keep, select = c("Cumulative_precipitation_t1_t7_mm", "doli.sum.helli",
                                          "Total_Phosphorus_ug", "Dissolved_P", "Dissolved_N"))
env.signif <- subset(env_keep, select = c("doli.sum.helli", "Cumulative_precipitation_t1_t7_mm"))

#the proportion of variation explained by the three constraining variables being 0.058

rda.signif <- rda(vir.rm~., data=env.signif)
summary(rda.signif, display=NULL)
#The explanatory variables now explain 11.67% of the variance in Y (species).

(R2adj <- RsquareAdj(rda.signif)$adj.r.squared)
#Here the strength of the relationship between X and Y corrected for the number of X variables is 0.0227.

#The significance of the model and of each canonical axis can be tested using the function anova 
#(note this is different from retaining significant variables as was done with forward selection, now we're testing the significance of the RDA axes):
anova.cca(rda.signif, step=1000)
#the RDA model is highly significant (p=0.006)
anova.cca(rda.signif, step=1000, by="axis")
# as well as both canonical axes.

#To visualize the results of the RDA, triplots can be drawn using the plot(). 
#Note that as in PCA, users can create scaling 1 and scaling 2 triplots. 
#In scaling 1, distance among objects approximate their Euclidean distances 
# in scaling 2, angles between variables X and Y reflect their correlation. 
#Thus, scaling 1 triplots can be used to interpret distances among objects 
#and scaling 2 triplots to interpret the relationships between X and Y. 

#Quick plots scaling 1 and 2
plot(rda.signif, scaling=1, main="Triplot RDA (scaling 1)")
plot(rda.signif, scaling=2, main="Triplot RDA (scaling 2)")

# #Advanced plots scaling 1
# plot(rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
# points(scores(rda.signif, display="sites", choices=c(1,2), scaling=1),
#        pch=21, col="black", bg="steelblue", cex=1.2)
# arrows(0,0,
#        scores(rda.signif, display="species", choices=c(1), scaling=1),
#        scores(rda.signif, display="species", choices=c(2), scaling=1),
#        col="black",length=0)
# text(scores(rda.signif, display="species", choices=c(1), scaling=1),
#      scores(rda.signif, display="species", choices=c(2), scaling=1),
#      labels=rownames(scores(rda.signif, display="species", scaling=1)),
#      col="black", cex=0.8)    
# arrows(0,0,
#        scores(rda.signif, display="bp", choices=c(1), scaling=1),
#        scores(rda.signif, display="bp", choices=c(2), scaling=1),
#        col="red")
# text(scores(rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
#      scores(rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
#      labels=rownames(scores(rda.signif, display="bp", choices=c(2), scaling=1)),
#      col="red", cex=1) 
# 
# #Advanced plots scaling 2
# plot(rda.signif, scaling=2, main="Triplot RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
# points(scores(rda.signif, display="sites", choices=c(1,2), scaling=2),
#        pch=20, col="steelblue", cex=1.2)
# arrows(0,0,
#        scores(rda.signif, display="species", choices=c(1), scaling=2)*2,
#        scores(rda.signif, display="species", choices=c(2), scaling=2)*2,
#        col="black",length=0)
# text(scores(rda.signif, display="species", choices=c(1), scaling=2)*2.1,
#      scores(rda.signif, display="species", choices=c(2), scaling=2)*2.1,
#      labels=rownames(scores(rda.signif, display="species", scaling=2)),
#      col="black", cex=0.8)    
# arrows(0,0,
#        scores(rda.signif, display="bp", choices=c(1), scaling=2),
#        scores(rda.signif, display="bp", choices=c(2), scaling=2),
#        col="red")
# text(scores(rda.signif, display="bp", choices=c(1), scaling=2)+0.05,
#      scores(rda.signif, display="bp", choices=c(2), scaling=2)+0.05,
#      labels=rownames(scores(rda.signif, display="bp", choices=c(2), scaling=2)),
#      col="red", cex=1)







##### db-RDA #####
# https://sites.ualberta.ca/~ahamann/teaching/renr690/Lab9b.pdf
#decide which distance measure to use by looking at the rank correlations between dissimilarity indices and gradient separation (the higher the value the better)
rankindex(env_keep, sp.asv.rm, indices = c("euc", "man", "gow", "bra", "kul"), stepacross = F, method = "spearman")

#remove site and temporal data from env_keep (do RDA without sites and time -- see PERMANOVA (far) below)
colnames(env_keep)
env_vars <- env_keep[,!(colnames(env_keep) %in% c("Months", "Years", "Site", "Period", "bloom2"))]
colnames(env_vars)

dbRDA <- capscale(vir.rm ~ ., data = env_keep, distance="bray")
plot(dbRDA)


#get the R2 for model fit for teh constrained ordinations
dbR2 <- RsquareAdj(dbRDA)$r.squared
dbR2

#adjusted R^2: adjusted R2 measures the unbiased amount of explained variation
dbR2adj <- RsquareAdj(dbRDA)$adj.r.squared
dbR2adj

# plot the RDA using ggplot (ggord package)
library(ggord)
ggord(dbRDA, env_keep$bloom2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  


# overall test of the significance of the analysis
anova(dbRDA)
#test axes for significance
anova(dbRDA, by="axis", perm.max=500)
#test for sig. environmental variables
anova(dbRDA, by="terms", permu=200)
#                                   Df SumOfSqs      F Pr(>F)    
# Months                             6   3.7920 5.0998  0.001 ***
# Years                              7   5.8803 6.7786  0.001 ***
# Site                               1   0.1032 0.8324  0.630    
# Period                             2   0.4735 1.9104  0.005 ** 
# bloom2                             1   0.1635 1.3195  0.170    
# Total_Phosphorus_ug                1   0.2037 1.6436  0.068 .  
# Total_Nitrogen_mg                  1   0.1820 1.4683  0.111    
# Dissolved_P                        1   0.1254 1.0122  0.414    
# Dissolved_N                        1   0.2277 1.8374  0.035 *  
# Cumulative_precipitation_t1_t7_mm  1   0.1861 1.5018  0.086 .  
# Mean_temperature_t0_t7             1   0.2035 1.6418  0.056 .  
# Dolicho.Abundance                  1   0.0633 0.5108  0.968    
# Micro.Abundance                    1   0.0736 0.5936  0.916    
# Cyano.Abundance                    1   0.1260 1.0170  0.437    
# Residual                          23   2.8503

#transforming negtaive eigenvalues (get rid of the -ve eigenvalues when doing the analysis)
#1) add a constant:
dbRDA_constant <- capscale(sp.asv.rm ~., env_vars, distance = "bray", add=TRUE)
plot(dbRDA_constant)
anova(dbRDA_constant)

#2) take the sqrt of dissimilarities
dbRDA_sqrt <- capscale(sp.asv.rm ~., env_vars, dist="bray", sqrt.dist=T)
plot(dbRDA_sqrt)
anova(dbRDA_sqrt)

#3) Do a square root transformation, Wisconsin double standardization (this emphasizes the environmental variables):
dbRDA_metaMDS=capscale(sp.asv.rm ~ ., env_vars, dist="bray", metaMDS=TRUE, sqrt.dist=TRUE) 
plot(dbRDA_metaMDS)
anova(dbRDA_metaMDS)



#select the significant explanatory variables
?ordiR2step
ordiR2step(rda(sp.asv.rm~1, data=env_vars), scope=formula(dbRDA), direction= "forward", R2scope=TRUE, pstep=1000)

#the proportion of variation explained by the constraining variables being 0.3542

#retain significant variables only
signif_env <- subset(env_vars, select = c("Total_Nitrogen_mg", "Mean_temperature_t0_t7", "Dissolved_N"))

#rda with significant variables 
dbrda_env_signif <- capscale(sp.asv.rm~., data=signif_env, dist="bray", sqrt.dist=T)
summary(dbrda_env_signif, display=NULL)
screeplot(dbrda_env_signif)
#The proportion of the variance of Y (species) explained by the X (env) variables = 37.93% (constrained) 
#the unexplained variance of Y = 62.07% (unconstrained)

#adjusted R^2
(adjR2 <- RsquareAdj(dbrda_env_signif)$adj.r.squared) 
#strength of the relationship between X and Y corrected for the number of X variables is 27.94%
#the explanatory variables explain 38.32% of the variance in Y (species)

#test significance of the RDA canonical axis (ie. constraining variables, ie. env var)
anova.cca(dbrda_env_signif, step=1000) 
#RDA model is highly significant (p=0.001)
anova.cca(dbrda_env_signif, step=1000, by="axis") 
#axis' 1-2 v significant(p=0.001), 3 (p=0.026)

#plot dbRDA
#Quick plots scaling 1 and 2
plot(dbrda_env_signif, scaling=2, main="Triplot RDA (scaling 2)")

#with ellipses
ggord(dbrda_env_signif, signif_env$bloom2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  







#https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html
#simple RDA
simpleRDA <- rda(vir.rm ~., env_keep)
summary(simpleRDA)
screeplot(simpleRDA)
# The total variance of the data set, partitioned into constrained and unconstrained variances, is a standard result. This result shows how much variation in your response variables was redundant with the variation in your explanatory variables. If the constrained variance is much higher than your unconstrained variance, the analysis suggests that much of the variation in the response data may be accounted for by your explanatory variables. If, however, there is a large proportion of unconstrained variation (i.e. variation in your response matrix that is non-redundant with the variation in the explanatory matrix), then the results should be interpreted with caution as only a small amount of the variation in your response matrix is displayed
# constrained axes (RDA axes)
# unconstrained axes (PCA axes)
#Object and response variable scores are often reported as "site" and "species" scores, respectively. These scores are the coordinates used to ordinate points and vectors. The coordinates of variables should be understood as the "tip" of their vector with the origin as its "tail". The direction of the vector is the direction of increase for that variable.
#Explanatory variable scores, also referred to as constraining variable scores, may be interpreted as response variable scores when the explanatory variable in question is quantitative. Scores for each state of nominal or factorial variables are the coordinates of these states' centroids and show the average position of the sites that have that state.
# first 2 axes explain most of the variation (RDA1: 0.3452 & RDA2: 0.2721) so plotting by these two axes represent the data well
# the unconstrained eigenvalue (PC1) is 0.09707, which is comparatively small, which means it does not display any important residenual structure of teh response data
# each RDA axis has an eigenvalue, which is the proportion of the variance explained by each axis
# The species and site scores = where the sites and species fall along the axes
# Partition of variance: the overall variance is partitioned into constrained and unconstrained fractions.
# constrained: the amount of variance the species by site matrix is explained by the explanatory variables (expressed as a proportion--is equivalent to R2 in a multiole regression). This R2 is biased so you have to look at the adjusted R2
# Eigenvalues and their contribute to variance: RDA1 -> RDA26 for the canonical axes, and unconstrained axes. The cumulative contribute of the variance is the proportion of the total variance of the response data explained by the RDA. The results also give the eigenvalues.
# Eigenvalues: the canonical eigenvalues are decreasing in value (in the order they are presented). But sometimes the residual structure (residual eigenvalue PC1) can be larger than the last RDA eigenvalue, which means that the residual structure has more variance than some of the structures that can be explained by the explan variables.
# Canonical eigenvalues: measure the amount of variance explained by the RDA model
# Residual eigenvalues: measure the amount of variance represented by the residual axes
# Accumulated constrained eigenvalues: cumulative amounts of variance expressed as proportions of the total explained variance
# Species scores: coordinates of the tips of vectors representing the response variables in bi or triplots
# Site scores (weighted sums of species scores): coordinates of the sites as expressed in the space of the response variables
# Site constraints (linear combinations of constraining variables): coordinates of the sites in space of the explanatory variables
# Biplot scores for constraining variables: coordinates of the tips of the vectors represnting explanatory variables
# Centroids for factor constraints: coordinates of centroids of levels of factor variables

#getting the canonical coefficient (ie. the equivalent of regression coefficients for each explanatory variable on each canonical axis)
coef(simpleRDA)

#get the R2 for model fit for teh constrained ordinations
R2 <- RsquareAdj(simpleRDA)$r.squared
R2

#adjusted R^2: adjusted R2 measures the unbiased amount of explained variation
R2adj <- RsquareAdj(simpleRDA)$adj.r.squared
R2adj #0.5084 ==> this model explains 50.8% of the variation in the data (if used the biased R2, any variable included in the explanatory responses would increase the R2, so the R2 needs to be adjusted for the number of eplanatory variables -- esp since we have 8)

#plot RDA
# two ways: using "wa" (weighted sums of species) and "lc" (fitted site scores)
# wa: more robust to noise in the environmental variables but are a step between constrained towards unconstrained. Default because more robust to noise in data
# lc: orthogonal linear combinations of the explanatory variable

#The interpretation of these plots depends on what scaling has been chosen. 
#In general, consider type I scaling if the distances between objects are of particular value or if most explanatory variables are binary or nominal.
# Consider type II scaling if the correlative relationships between variables are of more interest.

# https://mb3is.megx.net/gustame/constrained-analyses/rda
### Type I Scaling - Distance plots (object focused) ###
#Distances between object points approximate Euclidean distances. Thus, objects ordinated closer together can be expected to have similar variable values. This will not always hold true, as RDA only recovers part of the variation in the data set.
#Right-angled projections of object points onto vectors representing response variables approximate variable values for a given object.
#The angles between vectors representing response variables are meaningless.
#The angles between vectors representing response variables and those representing explanatory variables reflect their (linear) correlation.
#Note, that binary explanatory variables may be represented as points. These points are the centroids of objects which have a state "1" for a given binary variable. Projecting centroid points onto a vector representing a response variable reflects the relationship between these variables.
#Distances between centroids and between centroids and object points approximates Euclidean distances.

# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1
plot(simpleRDA, scaling=1, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")

# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(simpleRDA, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')

### Type II Scaling - Correlation plots (response variable focused) ###
#Distances between object points should not be considered to approximate Euclidean distances.
#Right-angled projections of object points onto vectors representing response variables approximate variable values for a given object.
#The angles between all vectors reflect their (linear) correlation. The correlation is equal to the cosine of the angle between vectors (e.g. a vector pair describing an angle of 90° are uncorrelated as cos(90) = 0), those describing an angle of 20° have strong, positive correlation as cos(20) =  0.94)
# Note, that binary or nominal explanatory variables may be represented as points. These points are the centroids of objects which have a state "1" for a given binary variable or realise a particular level of a nominal explanatory variable. Projecting centroid points onto a vector representing a response variable reflects the relationship between these variables.

# Scaling 2
plot(simpleRDA, main="Triplot RDA matrix ~ env - scaling 2 - wa scores")
spe2.sc <- scores(simpleRDA, choices=1:2, display="sp") # scores() choices= indicates which axes are to be selected, make sure to specify the scaling if its different than 2 
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')


# plot the RDA using ggplot (ggord package)
library(ggord)
ggord(simpleRDA, env_vars$bloom2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  


#plot results of RDA using "lc"
# site scores as linear combinations of the environmental variables 
# Scaling 1
plot(simpleRDA, scaling=1, display=c("lc", "sp", "cn"), main="Triplot RDA matrix ~ env -scaling 1- lc scores")
arrows(0,0, spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')

# scaling 2
plot(simpleRDA, display=c("sp", "lc", "cn"), main="Triplot RDA matrix ~ env -scaling2-lc scores", col=env_vars$bloom2)
arrows(0,0,spe2.sc[,1],spe2.sc[,2], length=0, lty=1,col='red')


# To choose the elements that are plotted, use the argument display=c(), sp=species, wa= site scores in the species space (weighted averages), lc= fitted site scores (linear combinations of explanatory variables) and cn= constraints (the explanatory variables)

## Fwd selection of variables ##
# variance inflation factors in the RDA
vif.cca(simpleRDA) #Anything above 10 should be examined or avoided

# RDA with all explanatory variables  
spe.rda.all <- rda(sp.asv.rm ~ Years + Site + Period + bloom2 + Total_Phosphorus_ug + Total_Nitrogen_mg + Dissolved_P + Dissolved_N + 
                     Cumulative_precipitation_t1_t7_mm + Mean_temperature_t0_t7 + Micro.Abundance, data=env_vars)

# Forward selection using ordistep
ordistep(rda(sp.asv.rm ~., data=env_vars), direction="forward", pstep=1000, R2scop=TRUE) #R2scope only accepts models with lower adjusted R2

# db-RDA
dbRDA2 <- capscale(sp.asv.rm ~ ., data=env_vars, distance = "bray")
dbRDA2
plot(dbRDA2)
summary(dbRDA2)
screeplot(dbRDA2)
#why use CAP: it allows for any dissimilarity measure & takes into account any correlation structure among the response variables (i.e. this is why its called ‘canonical’)








gen.imp <- t(vir_abund_helli)
env <- meta

#RDA requires complete data frames (i.e., no missing data)
sum(is.na(gen.imp))
#look at the structure of the df
str(env) 
#confirm that asv and enviro data are in the same order
identical(rownames(gen.imp), rownames(env))

#RDA is a regression-based method, and so can be subject to problems when using highly correlated predictors (Dormann et al., 2013). Generally, the |r| > 0.7 “rule of thumb” is a good guideline for removing correlated predictors. We will also check for multicollinearity using Variance Inflation Factors (VIF)
#Variable reduction should be guided by an ecological interpretation of the relevance of possible predictors. Here, we use the function pairs.panels to visualize correlations among our predictors. Correlation coefficients are in the upper right diagonal, with their size scaled to their |r|. The lower left shows scatter plots, while the diagonal shows histograms of the data
library(psych)
pairs.panels(env[,7:15], scale=T)
#remove total phsophorus and dissolved p
pred <- subset(env, select=-c(Dissolved_P, Total_Phosphorus_ug))

#rda
vir.rda <- rda(formula=gen.imp ~ Months + Years + Site + Period + bloom2 +
                 Total_Nitrogen_mg  + Dissolved_N + Cumulative_precipitation_t1_t7_mm + Mean_temperature_t0_t7+Cyano.Abundance, data = pred, scale=T) 





# *no RDA* for site and time. use PERMANOVA (combining adonis2 and betadisper)




##### MRT #####
# library(devtools)
# devtools::install_github("cran/mvpart")
# install_github("cran/MVPARTwrap", force=TRUE)
library(MVPARTwrap)

vir_helli #from phyloseq
vir_helli2 #from decostand
names(env)

#remove catagorical data from env (do RDA without sites and time -- see PERMANOVA (far) below)
str(env)

mrt_vars <- env[,(colnames(env) %in% c("Months", "Years", "Site", "Period"))] #test with nico
dim(mrt_vars)
env_vars <- mrt_vars

colnames(env_vars)

#rm NAs
env_keep <- env_vars[complete.cases(env_vars), ]
env_keep %>% dplyr::glimpse() 
summary(env_keep)

#### Remove viral asvs that are not present (due to removal of NA from env vars)
sp.asv <- t(vir_helli)
sp.asv <- vir_helli2

#rm sample rows that are not present in env_keep
vir.rm <- sp.asv[rownames(sp.asv) %in% rownames(env_keep),]

# species and enviro data without NAs needed for RDA
head(vir.rm)
dim(vir.rm)
head(env_keep)
dim(env_keep)

#use all data to identify which group accurately predicts 
mrt <- mvpart(as.matrix(vir.rm) ~ ., env_keep, 
              legend=F,
              margin=0.01,
              cp=0,
              xv = "pick", 
              xval=nrow(vir.rm),
              xvmult=100,
              which = 4)


mrt_month <- mvpart(vir.rm~ Months, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.975

mrt_year <- mvpart(as.matrix(vir.rm)~ Years, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.962

mrt_site <- mvpart(as.matrix(vir.rm)~ Site, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.98


mrt_period <- mvpart(as.matrix(vir.rm)~ Period, env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.976


mrt_bloom <- mvpart(as.matrix(vir.rm)~ bloom2, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.983


mrt_totP <- mvpart(as.matrix(vir.rm)~ Total_Phosphorus_ug, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.884


mrt_totN <- mvpart(as.matrix(vir.rm)~ Total_Nitrogen_mg, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.958


mrt_DisP <- mvpart(as.matrix(vir.rm)~ Dissolved_P, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.726
(DisP_R2 <- RsquareAdj(mrt_DisP)$adj.r.squared)
rpart.pca(mrt_DisP)

mrt_DisN <- mvpart(as.matrix(vir.rm)~ Dissolved_N, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.835


mrt_precip <- mvpart(as.matrix(vir.rm)~ Cumulative_precipitation_t1_t7_mm, env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.889


mrt_temp <- mvpart(as.matrix(vir.rm)~ Mean_temperature_t0_t7, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.962


mrt_cyano <- mvpart(as.matrix(vir.rm)~ cyano.sum.helli, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.836

mrt_doli <- mvpart(as.matrix(vir.rm)~ doli.sum.helli, env_keep,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.937

mrt_micro <- mvpart(as.matrix(vir.rm)~ micro.sum.helli, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.947


#significant vars based off of MRT scores
mrt_signif <- mvpart(as.matrix(vir.rm) ~ cyano.sum.helli+Cumulative_precipitation_t1_t7_mm+Dissolved_N+Dissolved_P+Total_Phosphorus_ug, 
                     env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)

# significant vars based off of ordi2step
signif_env
mrt_ordi <- mvpart(as.matrix(vir.rm) ~ Years+Months+Dissolved_N+Period+Dissolved_P, signif_env,
                   legend=T, margin=0.01, cp=0, xv="pick",
                   xval=nrow(vir.rm), xvmult=100, which=4, big.pts=T, bars=F)



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
