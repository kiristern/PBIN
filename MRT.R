#https://bookdown.org/forestgeoguest/mpart/mvpart.html

#run preprocessing from Initialize.R script
head(vir_abun_removed)
head(complete_env_keep)


formula <- vir_abun_removed ~ Months + Years + Site + Period + bloom2 + Tot_P + Tot_N + Dissolved_P + 
                             Dissolved_N + Cumul_precip + Avg_temp + Cyano.Abundance + Micro.Abundance + Dolicho.Abundance
#formula_signif <- vir_abun_removed ~ Months + Years + Period +  Cumul_precip + Avg_temp + bloom2 + Site + cyano_count

mrt <- mvpart(
  form=formula,
  data=complete_env_keep,  
  legend=FALSE, 
  margin=0.01, 
  cp=0, 
  xv="pick",
  xval=nrow(vir_abun_removed), 
  xvmult=100, 
  which=4)
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


tree_month <- mvpart(as.matrix(vir_abun_removed)~ Months, complete_env_keep,
               legend=T, margin=0.01, cp=0, xv="pick",
               xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_year <- mvpart(as.matrix(vir_abun_removed)~ Years, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

rpart.pca(tree_year)

tree_site <- mvpart(as.matrix(vir_abun_removed)~ Site, complete_env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_period <- mvpart(as.matrix(vir_abun_removed)~ Period, complete_env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_bloom <- mvpart(as.matrix(vir_abun_removed)~ bloom2, complete_env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_totP <- mvpart(as.matrix(vir_abun_removed)~ Tot_P, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_totN <- mvpart(as.matrix(vir_abun_removed)~ Tot_N, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_DisP <- mvpart(as.matrix(vir_abun_removed)~ Dissolved_P, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)
(DisP_R2 <- RsquareAdj(tree_DisP)$adj.r.squared)
rpart.pca(tree_DisP)

tree_DisN <- mvpart(as.matrix(vir_abun_removed)~ Dissolved_N, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_precip <- mvpart(as.matrix(vir_abun_removed)~ Cumul_precip, complete_env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_temp <- mvpart(as.matrix(vir_abun_removed)~ Avg_temp, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_cyano <- mvpart(as.matrix(vir_abun_removed)~ Cyano.Abundance, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_doli <- mvpart(as.matrix(vir_abun_removed)~ Dolicho.Abundance, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

tree_micro <- mvpart(as.matrix(vir_abun_removed)~ Micro.Abundance, complete_env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vir_abun_removed), xvmult=100, which=4, big.pts=T, bars=F)

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

