library(vegan)
library(mvpart)
library(MVPARTwrap)
library(rdaTest)
library(labdsv)
library(plyr)
library(MASS)

Calculate distance matrix
bc<-vegdist(abundance_removed, method="bray", binary=FALSE) 

# look at an unconstrained ordination first, it is always a good idea to look at both unconstrained and constrained ordinations
# set the seed: to reproduce the same result in the fture
set.seed(100)
bci.mds<-metaMDS(abundance_removed, distance = "bray", k = 2)

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # 0.1644892

# colour by bloom
ggplot(MDS_xy, aes(MDS1, MDS2, colour=complete_env_keep$bloom2)) +
  scale_shape_discrete(name="Bloom") + 
  geom_point() + 
  theme_bw() + 
  ggtitle('stress:0.12')

# RDA
simpleRDA <- rda(abundance_removed ~., data=complete_env_keep)
summary(simpleRDA)
screeplot(simpleRDA) #bstick not available for constrained ordinations

#canonical coefficient (i.e. the equivalent of regression coefficients for each explanatory variable on each canonical axis)
coef(simpleRDA)

# unadjusted R^2 retreived from the rda result
(R2 <- RsquareAdj(simpleRDA)$r.squared)

# adjusted R^2
(R2adj <- RsquareAdj(simpleRDA)$adj.r.squared)
#The adjusted R2 measures the unbiased amount of explained variation. So this model explains 45.6% of the variation in the data. If you used the biased R2, any variable included in the explanatory responses would increase the R2, so the R2 needs to be adjusted for the number of explanatory variables (especially since we have eight included here).

# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1
plot(simpleRDA, scaling=1, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")

# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(simpleRDA, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')

# Scaling 2
plot(simpleRDA, main="Triplot RDA matrix ~ env - scaling 2 - wa scores")
spe2.sc <- scores(simpleRDA, choices=1:2, display="sp") # scores() choices= indicates which axes are to be selected, make sure to specify the scaling if its different than 2 
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')

# plot the RDA using ggplot (ggord package)
library(ggord)
ggord(simpleRDA, complete_env_keep$bloom2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  

# plot the results of the RDA using ‘lc’
# site scores as linear combinations of the environmental variables 
# Scaling 1
plot(simpleRDA, scaling=1, display=c("lc", "sp", "cn"), main="Triplot RDA matrix ~ env -scaling 1- lc scores")
arrows(0,0, spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')


#text(simpleRDA, display = "spec", cex=0.7, col="blue")
# scaling 2
plot(simpleRDA, display=c("sp", "lc", "cn"), main="Triplot RDA matrix ~ env -scaling2-lc scores", col=complete_env_keep$bloom2)
arrows(0,0,spe2.sc[,1],spe2.sc[,2], length=0, lty=1,col='red')

# To choose the elements that are plotted, use the argument display=c(), sp=species, wa= site scores in the species space (weighted averages), lc= fitted site scores (linear combinations of explanatory variables) and cn= constraints (the explanatory variables)
#The results look very similar. We will stick with ‘wa’ because they are the suggested scores.

#When you want to reduce the number of explanatory variables.
# variance inflation factors in the RDA
vif.cca(simpleRDA)
#omit anything over 10









