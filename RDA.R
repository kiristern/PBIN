#run preprocessing from begining of Analysis.R script
abundance_removed
complete_env_keep


# dim(abundance_removed)
# str(abundance_removed)
# head(abundance_removed)
# summary(abundance_removed)


##### Environmental data #######
names(complete_env_keep)
dim(complete_env_keep)
str(complete_env_keep)
head(complete_env_keep)
summary(complete_env_keep)
#pairs(env_keep, main="Bivariate Plots of the Environmental Data" ) 

#View mean and std_dev of standardized environmental data
apply(complete_env_keep[,c(6:11)], 2, mean) #centered data (mean~0)
apply(complete_env_keep[,c(6:11)], 2, sd) #scaled data (std dev=1)

######## RDA #######

#Plotting with ordistep2: only Months and Years are significant

 # #Use adonis to find significant environmental variables
(abund_table.adonis <- adonis(abundance_removed ~ Months + Years + Site + Period + bloom2 +
                                  Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp,
                                  data=complete_env_keep, permutations = 9999))
 #Extract the best variables
(bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=0.05])
#Remove last two NA entries
(bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)])
#Only use those environmental variables in cca that were found significant
#CCA estimates species optima, regression coefficients, and site scores using a Gaussian weighted-averaging model combined with regression.
(eval(parse(text=paste("sol <- rda(abundance_removed ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=complete_env_keep)",sep=""))))
#You can use the following to use all the environmental variables
#sol<-cca(abund_table ~ ., data=meta_table2)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"),scaling=2)
#Extract site data first
df_sites<-data.frame(scrs$sites,bloom=as.factor(complete_env_keep[,5]),Site=as.factor(complete_env_keep[,3]),Months=as.factor(complete_env_keep[,1]))
colnames(df_sites)<-c("x","y","bloom","Site","Months")
#Draw sites
(p <- ggplot() +
  geom_point(data=df_sites,aes(x,y,colour=bloom),size=2) +
  labs(x="RDA1(%)",y="RDA2(%)"))
#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines
# which are the coordinates of the heads of unit length vectors. In plot these are
# scaled by their correlation (square root of the column r2) so that "weak" predictors
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths
# using command scores. The plotted (and scaled) arrows are further adjusted to the
# current graph using a constant multiplier: this will keep the relative r2-scaled
# lengths of the arrows but tries to fill the current plot. You can see the multiplier
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul.
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

(p <- p + geom_segment(data=df_arrows,
                       aes(x = 0, y = 0, xend = x, yend = y),
                       arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.6))

(p2 <- p + geom_text(data=as.data.frame(df_arrows*1.3),
                    aes(x, y, label = rownames(df_arrows),size=3,face="bold"),
                    color="black",alpha=0.6, size=2, fontface="bold"))

# Draw species -- NOT INFORMATIVE, TOO MANY CAN'T READ
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")
# Either choose text or points
# (p<-p+geom_text(data=df_species,aes(x,y,label=rownames(df_species))))
# (p<-p+geom_point(data=df_species,aes(x,y,shape="Species"))+scale_shape_manual("",values=2))

(p2<-p2+theme_bw()+geom_text(size=3, fontface="bold"))

(p3<-p2+scale_colour_manual(values = c("red","blue", "green")) + 
    theme(panel.background = element_rect(fill = "white")))

(p4<-p3+theme(axis.text.x  = element_text(vjust=0.5, size=12), 
              axis.text.y  = element_text(vjust=0.5, size=12), 
              axis.title.x = element_text(size = 15, face = "bold", color="black"),
              axis.title.y = element_text(size=15, face = "bold",color="black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1)))

RsquareAdj(sol)
anova(sol)

################################
#https://wiki.qcbs.ca/r_workshop10

env.z<-subset(complete_env_keep)

spe.rda <- rda(abundance_removed~., data=env.z)

# rda_spe <- rda(abundance_removed ~ Months + Years + Site + Period + bloom2 + Tot_P + 
#                  Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp, 
#                data = complete_env_keep)

#extract results
summary(spe.rda, display = NULL)
# summary(rda_spe, display = NULL)
  #the explanatory variables (env vars) explain 38.02% of the variance in Y (species)

#select the significant explanatory variables
?ordiR2step
ordiR2step(rda(abundance_removed~1, data=env.z), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)
  #the proportion of variation explained by the constraining variables being 0.3542

#retain significant variables only
signif_env <- subset(complete_env_keep, select = c("Months", "Years"))

#rda with significant variables 
rda_spe_signif <- rda(abundance_removed~., data=complete_env_keep)
summary(rda_spe_signif, display=NULL)
screeplot(rda_spe_signif)
  #The proportion of the variance of Y (species) explained by the X (env) variables = 50.47% (constrained) 
  #the unexplained variance of Y = 49.53% (unconstrained)

#adjusted R^2
(adjR2 <- RsquareAdj(rda_spe_signif)$adj.r.squared) 
  #strength of the relationship between X and Y corrected for the number of X variables is 27.94%
  #the explanatory variables explain 38.32% of the variance in Y (species)

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
plot(rda_spe_signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(rda_spe_signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
# arrows(0,0,
#        scores(rda_spe_signif, display="species", choices=c(1), scaling=1),
#        scores(rda_spe_signif, display="species", choices=c(2), scaling=1),
#        col="black",length=0)
# text(scores(rda_spe_signif, display="species", choices=c(1), scaling=1),
#      scores(rda_spe_signif, display="species", choices=c(2), scaling=1),
#      labels=rownames(scores(rda_spe_signif, display="species", scaling=1)),
#      col="black", cex=0.8)
arrows(0,0,
       scores(rda_spe_signif, display="bp", choices=c(1), scaling=1),
       scores(rda_spe_signif, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(rda_spe_signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(rda_spe_signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(rda_spe_signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#Advanced plots scaling 2
  #angles ebtween var X and Y relect their correlation
  #used to interpret the relationship betwen X and Y
windows()
plot(rda_spe_signif, scaling=2, main="Triplot RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(rda_spe_signif, display="sites", choices=c(1,2), scaling=2),
       pch=21, col="black", bg="steelblue", cex=1.2)
# arrows(0,0,
#        scores(rda_spe_signif, display="species", choices=c(1), scaling=2)*2,
#        scores(rda_spe_signif, display="species", choices=c(2), scaling=2)*2,
#        col="black",length=0)
# text(scores(rda_spe_signif, display="species", choices=c(1), scaling=2)*2.1,
#      scores(rda_spe_signif, display="species", choices=c(2), scaling=2)*2.1,
#      labels=rownames(scores(rda_spe_signif, display="species", scaling=2)),
#      col="black", cex=0.8)    
arrows(0,0,
       scores(rda_spe_signif, display="bp", choices=c(1), scaling=2),
       scores(rda_spe_signif, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(rda_spe_signif, display="bp", choices=c(1), scaling=2)+0.05,
     scores(rda_spe_signif, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(rda_spe_signif, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)  











