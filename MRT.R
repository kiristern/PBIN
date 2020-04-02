#https://bookdown.org/forestgeoguest/mpart/mvpart.html

abundance_removed
complete_env_keep


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

