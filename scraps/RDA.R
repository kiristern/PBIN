#run preprocessing from begining of Analysis.R script
vir_abun_removed
complete_env_keep



#how to read results:
#http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html 


#the following inspired from:
#https://popgen.nescent.org/2018-03-27_RDA_GEA.html

#confirm that asv and env data are in the same order
identical(rownames(vir_abun_removed), rownames(complete_env_keep))

# #visualize correlations among predictors
# pairs.panels(complete_env_keep[,6:12], scale=T)
#   #totP and totN are correlated at 0.80; can remove one

# enviro_var<-meta
# summary(enviro_var)
# env_var <- enviro_var %>% select(1:6,11:12)
# env_var<-env_var[complete.cases(env_var),]
# dim(env_var)
# dim(abundance_removed)
# str(abundance_removed)
# head(abundance_removed)
# summary(abundance_removed)


rda1 <- rda(formula=vir_abun_removed ~ Months + Years + Site + Period + bloom2 +
  Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp, data=complete_env_keep, scale=T)
  #"proportion" col for "constrained" (0.3807) is equivalent to the R2 of a multiple regression. This will be biased and should be adjusted:
RsquareAdj(rda1)
  #the constrained ordination explains about 14% of the variation

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:
summary(eigenvals(rda1, model="constrained"))
#to visualize this information using a screeplot of the canonical eigenvalues
screeplot(rda1)
  #can see the first two constrained axes explain most of the variance

#check rda model for significance. The null hypothesis is that no linear relationship exists between ASV data and env predictors
signif.full <- anova.cca(rda1, parallel = getOption("mc.cores"), permutation=999)
signif.full
  #the full model is significant. Doesn't tell us much. can check each constrained axis for significance:
signif.axis <- anova.cca(rda1, by="axis", parallel = getOption("mc.cores"))
signif.axis
  # Forward tests for axes
  # Permutation: free
  # Number of permutations: 999
  # 
  # Model: rda(formula = vir_abun_removed ~ Months + Years + Site + Period + bloom2 + Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp + cyano_count, data = complete_env_keep, scale = T)
  # Df Variance      F Pr(>F)    
  # RDA1      1   161.49 4.2450  0.001 ***
  # RDA2      1   145.54 3.8257  0.001 ***
  # RDA3      1   105.21 2.7655  0.001 ***
  # RDA4      1   100.33 2.6373  0.001 ***
  # RDA5      1    78.84 2.0724  0.065 .  
  # RDA6      1    63.53 1.6700  0.667    
  # RDA7      1    58.97 1.5502  0.848    
  # RDA8      1    50.90 1.3380  0.995    
  # RDA9      1    48.94 1.2864  0.996    
  # RDA10     1    43.63 1.1468  1.000    
  # RDA11     1    40.83 1.0734  1.000    
  # RDA12     1    37.12 0.9756  1.000    
  # RDA13     1    33.98 0.8933  1.000    
  # RDA14     1    31.67 0.8325  1.000    
  # RDA15     1    29.51 0.7757  1.000    
  # RDA16     1    28.11 0.7390  1.000    
  # RDA17     1    21.96 0.5771  1.000    
  # RDA18     1    18.43 0.4845  0.997    
  # Residual 47  1788.01                  
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Check Variance Inflation Factors for predictor variables used in teh model:
vif.cca(rda1)
  #can remove some of the higher VIF values (ie. Months)

#Plot RDA using scaling=3 (symmetrical scaling), scales ASV and individual scores by the sqrt of the eigenvalues
plot(rda1, scaling=2) #default is axes 1 and 2
# plot(rda1, choices=c(1,3), scaling=3) #axes 1 and 3



#using adonis results (below):
rda2 <- rda(formula=vir_abun_removed ~ Months + Years  + Period + bloom2 +
              Tot_P +  Dissolved_N + Avg_temp, data=complete_env_keep, scale=T)
#"proportion" col for "constrained" (0.3153) is equivalent to the R2 of a multiple regression. This will be biased and should be adjusted:
RsquareAdj(rda2)
  #the constrained ordination explains about 14% of the variation

#The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:
summary(eigenvals(rda2, model="constrained"))
#to visualize this information using a screeplot of the canonical eigenvalues
screeplot(rda2)

plot(rda2, scaling=2) #default is axes 1 and 2
  #placement of env vars indicate their loading on the two displayed RDA axes. 
    #Dissolved_N is loading heavily on RDA1 indicating that this var explains a larger portion of the variance associated with axis 1
    #the location of the species relative tot he env vars indicates how strongly a species is associated with a particular env variable
plot(rda2, choices=c(1,3), scaling=3) #axes 1 and 3





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
(abund_table.adonis <- adonis(vir_abun_removed ~ Months + Years + Site + Period + bloom2 +
                                  Tot_P + Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp,
                                  data=complete_env_keep, permutations = 9999))
 #Extract the best variables
(bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=0.05])
#Remove last two NA entries
(bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)])
#Only use those environmental variables in cca that were found significant
#CCA estimates species optima, regression coefficients, and site scores using a Gaussian weighted-averaging model combined with regression.
bestenv <- (eval(parse(text=paste("sol <- rda(vir_abun_removed ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=complete_env_keep)",sep=""))))
#You can use the following to use all the environmental variables
sol<-cca(vir_abun_removed ~ ., data=complete_env_keep)
scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"),scaling=2)

best_scrs<-scores(bestenv,display=c("sp","wa","lc","bp","cn"),scaling=2)

#Extract site data first
df_sites<-data.frame(best_scrs$sites)#,bloom=as.factor(complete_env_keep[,5]),Site=as.factor(complete_env_keep[,3]),Months=as.factor(complete_env_keep[,1]))
colnames(df_sites)<-c("x","y")#,"bloom","Site","Months")
#Draw sites
(p <- ggplot() +
  geom_point(data=df_sites,aes(x,y),size=2) +
  labs(x="RDA1(%)",y="RDA2(%)"))
#Draw biplots
multiplier <- vegan:::ordiArrowMul(best_scrs$biplot)
# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines
# which are the coordinates of the heads of unit length vectors. In plot these are
# scaled by their correlation (square root of the column r2) so that "weak" predictors
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths
# using command scores. The plotted (and scaled) arrows are further adjusted to the
# current graph using a constant multiplier: this will keep the relative r2-scaled
# lengths of the arrows but tries to fill the current plot. You can see the multiplier
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul.
df_arrows<- best_scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

(p <- p + geom_segment(data=df_arrows,
                       aes(x = 0, y = 0, xend = x, yend = y),
                       arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.6))

(p2 <- p + geom_text(data=as.data.frame(df_arrows*1.3),
                    aes(x, y, label = rownames(df_arrows),size=3,face="bold"),
                    color="black",alpha=0.6, size=2, fontface="bold")+
  theme_classic())

# Draw species -- NOT INFORMATIVE, TOO MANY CAN'T READ
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")
# Either choose text or points
# (p<-p+geom_text(data=df_species,aes(x,y,label=rownames(df_species))))
# (p<-p+geom_point(data=df_species,aes(x,y,shape="Species"))+scale_shape_manual("",values=2))

#graph aesthetics
# (p2<-p2+theme_bw()+geom_text(size=3, fontface="bold"))
# 
# (p3<-p2+scale_colour_manual(values = c("red","blue", "green")) + 
#     theme(panel.background = element_rect(fill = "white")))
# 
# (p4<-p3+theme(axis.text.x  = element_text(vjust=0.5, size=12), 
#               axis.text.y  = element_text(vjust=0.5, size=12), 
#               axis.title.x = element_text(size = 15, face = "bold", color="black"),
#               axis.title.y = element_text(size=15, face = "bold",color="black"),
#               panel.border = element_rect(colour = "black", fill=NA, size=1)))

RsquareAdj(sol)
anova(sol)








################################
#https://wiki.qcbs.ca/r_workshop10

spe.rda <- rda(vir_abun_removed~., data=complete_env_keep)

# rda_spe <- rda(abundance_removed ~ Months + Years + Site + Period + bloom2 + Tot_P + 
#                  Tot_N + Dissolved_P + Dissolved_N + Cumul_precip + Avg_temp, 
#                data = complete_env_keep)

#extract results
summary(spe.rda, display = NULL)
# summary(rda_spe, display = NULL)
  #the explanatory variables (env vars) explain 38.02% of the variance in Y (species)

#select the significant explanatory variables
?ordiR2step
ordiR2step(rda(vir_abun_removed~1, data=complete_env_keep), scope= formula(spe.rda), direction= "forward", R2scope=TRUE, pstep=1000)

  #the proportion of variation explained by the constraining variables being 0.3542

#retain significant variables only
signif_env <- subset(complete_env_keep, select = c("Months", "Years", "bloom2", "Dissolved_N", "Tot_P", "Avg_temp"))

#rda with significant variables 
rda_spe_signif <- rda(vir_abun_removed~., data=signif_env)
summary(rda_spe_signif, display=NULL)
screeplot(rda_spe_signif)
  #The proportion of the variance of Y (species) explained by the X (env) variables = 37.93% (constrained) 
  #the unexplained variance of Y = 62.07% (unconstrained)

#adjusted R^2
(adjR2 <- RsquareAdj(rda_spe_signif)$adj.r.squared) 
  #strength of the relationship between X and Y corrected for the number of X variables is 27.94%
  #the explanatory variables explain 38.32% of the variance in Y (species)

#test significance of the RDA canonical axis (ie. constraining variables, ie. env var)
?anova.cca
anova.cca(rda_spe_signif, step=1000) 
  #RDA model is highly significant (p=0.001)
anova.cca(rda_spe_signif, step=1000, by="axis") 
  #axis' 1-2 v significant(p=0.001), 3 (p=0.026)

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





