library(corrplot)
library(ggpubr)

## BETWEEN 2 PTS ##

#correlation anylysis: The Spearman correlation method computes the correlation between the rank of x and the rank of y variables.
bacteria <- cyano_samples #loaded from cyano.R hellinger transformed
dim(bacteria)
#remove rows (ASVs) that have no data (row all zero)
bact_no0 <- bacteria[apply(bacteria[,-1], 1, function(x) !all(x==0)),]
dim(bact_no0)

viral <- t(vir_abun_removed) #loaded from Initialize.R hellinger transformed
dim(viral)
#remove rows (ASVs) that have no data (row all zero)
vir_no0 <- viral[apply(viral[,-1], 1, function(x) !all(x==0)),]
dim(vir_no0)

bact_no0 <- t(bact_no0)
vir_no0 <- t(vir_no0)

#keep date only (ie. remove everything before first period)
#change _ to .
row.names(vir_no0) <- gsub("_", ".", row.names(vir_no0))
#remove everything before 1st period (just to keep date)
row.names(vir_no0) <- gsub("^.*?\\.","", row.names(vir_no0))
#add cyano for all cyano ASVs
colnames(vir_no0) <- lapply(colnames(vir_no0), function(x) paste("cyano", x, sep = "_"))

rownames(bact_no0) %in% rownames(vir_no0)
rownames(vir_no0) %in% rownames(bact_no0)

#merge bacteria and viral
bv <- merge(bact_no0, vir_no0, by="row.names")
bv <- bv[,-1]

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
vb_cor = cor(bv, method = c("spearman"), use = "complete.obs") #If data contain missing values, use = "complete.obs" to handle missing values by case-wise deletion.
head(vb_cor)
round(vb_cor, 2)

#get significance levels (p-values)
library(Hmisc)

vb_rcorr = rcorr(as.matrix(bv), type=c("spearman"))
#The output of the function rcorr() is a list containing the following elements : 
# r : the correlation matrix
# n : the matrix of the number of observations used in analyzing each pair of variables
# P : the p-values corresponding to the significance levels of correlations.

# extract correlation coefficient
vb_coeff = vb_rcorr$r

dim(vb_coeff)
colnames(vb_coeff)
#keep only bacteria-virus correlation (no virus-virus or bact-bact)
vb_coeff_rem <- vb_coeff[343:901,1:342]

#set NaN to zero
vb_coeff_rem[is.nan(vb_coeff_rem)] <- 0
head(vb_coeff_rem)

# extract p-values
vb_pval = vb_rcorr$P

#keep only bacteria-virus correlation (no virus-virus or bact-bact)
vb_pval_rem <- vb_pval[343:901,1:342]

#adjust for multiple comparisons
pval_adj = p.adjust(vb_pval_rem, method=c("fdr"))

head(vb_pval_rem)


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# # cormat : matrix of the correlation coefficients
# # pmat : matrix of the correlation p-values
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  = (cormat)[ut],
#     p = pmat[ut]
#   )
# }

# corr_table = flattenCorrMatrix(vb_coeff_rem, vb_pval_rem)
# head(corr_table)

#flatten correlation Matrix to df
corr_table = data.frame(row=rownames(vb_coeff_rem)[row(vb_coeff_rem)], col=colnames(vb_coeff_rem)[col(vb_coeff_rem)], corr=c(vb_coeff_rem))
#same for p-value
pval_tab = data.frame(row=rownames(vb_pval_rem)[row(vb_pval_rem)], col=colnames(vb_pval_rem)[col(vb_pval_rem)], pval=c(vb_pval_rem))
#add pval to correlation df
corr_table$pval <- pval_tab$pval
head(corr_table)

#remove correlations between -0.6-0.6 & p-val >0.5
library(data.table)
corr_tab_adj = as.data.frame(setDT(corr_table)[!(corr %between% c(-0.6, 0.6) | pval > c(0.05))])
head(corr_tab_adj)

#turn back into matrix
mat_corr_adj <- with(corr_tab_adj, {
  out <- matrix(nrow=nlevels(row), ncol=nlevels(col),
                dimnames=list(levels(row), levels(col)))
  out[cbind(row, col)] <- corr
  out
})

#set NaN to zero
mat_corr_adj[is.na(mat_corr_adj)] <- 0
#remove rows that are completely empty
mat_corr_no0 = mat_corr_adj[apply(mat_corr_adj[,-1], 1, function(x) !all(x==0)),]
dim(mat_corr_no0)
dim(mat_corr_adj)

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
head(vb_pval_rem)

#visualize
colourpalette <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(vb_coeff_rem,  method="color", col=colourpalette(200),  
         type="full", 
         # order="hclust", #reorder: hierarchical clustering according to the correlation coeff
        #addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         p.mat = vb_pval_rem,     #add significance level to the correlogram
         sig.level = 0.05, #correlations with p-value > 0.05 are considered as insignificant.
         insig = "blank", #leave blank on no significant coeff
         diag=FALSE      # hide correlation coefficient on the principal diagonal

)
