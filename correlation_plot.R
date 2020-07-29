library(corrplot)
library(ggpubr)

## BETWEEN 2 PTS ##

#correlation anylysis: The Spearman correlation method computes the correlation between the rank of x and the rank of y variables.
#using output of gLV_preprocess.R script
viral_bacterial <- read.csv("gLV_table.csv", header = T, row.names = 1)

#transform: hellinger
vb_helli <- decostand(viral_bacterial, method="hellinger")

#visualize data using scatter plots
ggscatter(vb_helli, x = "ASV_605", y = "micro_ASV_143",
          add = "reg.line", conf.int = T,
          cor.coef=T, cor.method = "spearman")
  #In the situation where the scatter plots show curved patterns, we are dealing with nonlinear association between the two variables.
  # Are the data from each of the 2 vars (x,y) follow a normal distribution? Use Shapiro-Wilk normality test. Null hypothesis: the data are normally distributed. Alternative hypothesis: the data are not normally distributed
  shapiro.test(vb_helli$ASV_605) # p < 2.2e-16
  shapiro.test(vb_helli$micro_ASV_143) # p < 2.2e-16
  # the two p-values are << 0.05 implying that the distribution of the data are significantly different from normal distribution. 
  # visual inspection of the data normality
  ggqqplot(vb_helli$ASV_605)
  ggqqplot(vb_helli$micro_ASV_143)
  # if the data are not normally distributed, it’s recommended to use the non-parametric correlation, including Spearman and Kendall rank-based correlation tests

# Spearman’s rho statistic is also used to estimate a rank-based measure of association. This test may be used if the data do not come from a bivariate normal distribution.
spear_cor <-cor.test(vb_helli$ASV_605, vb_helli$micro_ASV_143,  method = "spearman")
  #rho is the Spearman’s correlation coefficient.
  #the correlation coefficient between x and y is 0.02600654 and the p-value is 0.7743
  #-1 indicates a strong negative correlation : this means that every time x increases, y decreases 
  #0 means that there is no association between the two variables (x and y) 
  #1 indicates a strong positive correlation : this means that y increases with x 






## WHOLE MATRIX ##

#compute correlation
vb_cor = cor(vb_helli, method = c("spearman"), use = "complete.obs") #If data contain missing values, use = "complete.obs" to handle missing values by case-wise deletion.
head(vb_cor)
round(vb_cor, 2)

#get significance levels (p-values)
library(Hmisc)

vb_rcorr = rcorr(as.matrix(vb_helli), type=c("spearman"))
#The output of the function rcorr() is a list containing the following elements : 
# r : the correlation matrix
# n : the matrix of the number of observations used in analyzing each pair of variables
# P : the p-values corresponding to the significance levels of correlations.

# extract correlation coefficient
vb_coeff = vb_rcorr$r

#set NaN to zero
vb_coeff[is.nan(vb_coeff)] <- 0
head(vb_coeff)

# extract p-values
vb_pval = vb_rcorr$P
#adjust for multiple comparisons
pval_adj = p.adjust(vb_pval, method=c("fdr"))

head(vb_pval)


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut],
    p = pmat[ut]
  )
}

corr_table = flattenCorrMatrix(vb_coeff, vb_pval)
head(corr_table)



# #compute the matrix of p-value
# # mat : is a matrix of data
# # ... : further arguments to pass to the native R cor.test function
# cor.mtest <- function(mat, ...) {
#   mat <- as.matrix(mat)
#   n <- ncol(mat)
#   p.mat<- matrix(NA, n, n)
#   diag(p.mat) <- 0
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       tmp <- cor.test(mat[, i], mat[, j], ...)
#       p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#     }
#   }
#   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#   p.mat
# }
# # matrix of the p-value of the correlation
# p.mat <- cor.mtest(vb_helli)
# 
# head(p.mat)
head(vb_pval)

#visualize
colourpalette <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(vb_coeff,  method="color", col=colourpalette(200),  
         type="upper", 
         order="hclust", #reorder: hierarchical clustering according to the correlation coeff
        #addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         p.mat = vb_pval,     #add significance level to the correlogram
         sig.level = 0.05, #correlations with p-value > 0.05 are considered as insignificant. 
         insig = "blank", #leave blank on no significant coeff
         diag=FALSE      # hide correlation coefficient on the principal diagonal

)
