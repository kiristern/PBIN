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
colnames(bact_no0) <- lapply(colnames(bact_no0), function(x) paste("cyano", x, sep = "_"))

rownames(bact_no0) %in% rownames(vir_no0)
rownames(vir_no0) %in% rownames(bact_no0)

#merge bacteria and viral
bv <- merge(bact_no0, vir_no0, by="row.names")
bv <- bv[,-1]

# #visualize data using scatter plots
# ggscatter(vb_helli, x = "ASV_605", y = "micro_ASV_143",
#           add = "reg.line", conf.int = T,
#           cor.coef=T, cor.method = "spearman")
#   #In the situation where the scatter plots show curved patterns, we are dealing with nonlinear association between the two variables.
#   # Are the data from each of the 2 vars (x,y) follow a normal distribution? Use Shapiro-Wilk normality test. Null hypothesis: the data are normally distributed. Alternative hypothesis: the data are not normally distributed
#   shapiro.test(vb_helli$ASV_605) # p < 2.2e-16
#   shapiro.test(vb_helli$micro_ASV_143) # p < 2.2e-16
#   # the two p-values are << 0.05 implying that the distribution of the data are significantly different from normal distribution. 
#   # visual inspection of the data normality
#   ggqqplot(vb_helli$ASV_605)
#   ggqqplot(vb_helli$micro_ASV_143)
#   # if the data are not normally distributed, it’s recommended to use the non-parametric correlation, including Spearman and Kendall rank-based correlation tests
# 
# # Spearman’s rho statistic is also used to estimate a rank-based measure of association. This test may be used if the data do not come from a bivariate normal distribution.
# spear_cor <-cor.test(vb_helli$ASV_605, vb_helli$micro_ASV_143,  method = "spearman")
#   #rho is the Spearman’s correlation coefficient.
#   #the correlation coefficient between x and y is 0.02600654 and the p-value is 0.7743
#   #-1 indicates a strong negative correlation : this means that every time x increases, y decreases 
#   #0 means that there is no association between the two variables (x and y) 
#   #1 indicates a strong positive correlation : this means that y increases with x 






## WHOLE MATRIX ##

# #compute correlation
# vb_cor = cor(bv, method = c("spearman"), use = "complete.obs") #If data contain missing values, use = "complete.obs" to handle missing values by case-wise deletion.
# head(vb_cor)
# round(vb_cor, 2)
# 
# colnames(vb_cor)
# vb_cor_rem <- vb_cor[343:901,1:342]



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

#remove correlations between -0.6-0.6 & p-val > 0.5
library(data.table)
corr_tab_adj = setDT(corr_table)[!(corr %between% c(-0.6, 0.6) | pval > 0.05)]
head(corr_tab_adj)

dim(corr_table)
dim(corr_tab_adj)

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


#repeat for pval
mat_pval_adj <- with(corr_tab_adj, {
  out <- matrix(nrow=nlevels(row), ncol=nlevels(col),
                dimnames=list(levels(row), levels(col)))
  out[cbind(row, col)] <- pval
  out
})
#set NaN to zero
mat_pval_adj[is.na(mat_pval_adj)] <- 0
#remove rows that are completely empty
mat_pval_no0 = mat_pval_adj[apply(mat_pval_adj[,-1], 1, function(x) !all(x==0)),]
dim(mat_pval_no0)
dim(mat_pval_adj)


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

#visualize
colourpalette <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(mat_corr_no0,  method="color", col=colourpalette(200),  
         type="full", 
         # order="hclust", #reorder: hierarchical clustering according to the correlation coeff
        #addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         p.mat = mat_pval_no0,     #add significance level to the correlogram
         sig.level = 0.05, #correlations with p-value > 0.05 are considered as insignificant.
         insig = "blank", #leave blank on no significant coeff
         diag=FALSE      # hide correlation coefficient on the principal diagonal

)



#get all the ASVs with high correlation only
corr_table
top_corrs = corr_table[with(corr_table, !(corr <= 0.8))] #&& (pval >= 0.05))]
top_corrs = as.data.frame(top_corrs)
top_corrs

#set zeros to NA
top_corrs[top_corrs==0] <- NA
#remove NA rows
top_corrs <- na.omit(top_corrs)
top_corrs

#extract cyanoASV
cyanoASV <- top_corrs$col
#unique ASV only
uniqueCyano <- as.data.frame(unique(cyanoASV))
uniqueCyano

#extract viral ASV
viralASV <- top_corrs$row
#unique ASV only
uniqueVir <- as.data.frame(unique(viralASV))
uniqueVir

#keep only top unique ASV from samples df
topUniqB <- subset(bact_no0, uniqueCyano$`unique(cyanoASV)` %in% colnames(bact_no0))

#keep only top unique ASV from samples df
topUniqV <- subset(vir_no0, uniqueVir$`unique(viralASV)` %in% colnames(vir_no0))

#merge top df
topUniqVB <- merge(topUniqB, topUniqV, by="row.names")
#set rowname as sample
topVB <- topUniqVB %>% remove_rownames() %>% column_to_rownames(var="Row.names")
topVB <- t(topVB)

#remove rows that are entirely 0
topVB_no0 = topVB[apply(topVB[,-1], 1, function(x) !all(x==0)), ]
(topVB_no0 <- t(topVB_no0))

#remove everything after third period in sample date
rownames(topVB_no0) <- sub(".[^.]+$", "", rownames(topVB_no0))

topVB_ts <- ts(topVB_no0)





# libraries
library(dygraphs)
library(xts) # To make the convertion data-frame / xts format
library(lubridate) # You will love it to work with dates

# Load the data
topVB_no0.2 <- topVB_no0

# Check the format, it is not a date yet !
str(rownames(topVB_no0.2))

# The wanna-be-date column.
dmy_form <- as.Date(rownames(topVB_no0.2), "%d.%m.%Y")

rownames(topVB_no0.2) <- dmy_form
bv2 <- bv

row.names(bv2) <- dmy_form

str(dmy_form)


# Check if it worked properly!
str(data)

# It does! Let's go to the its format like seen above, and make the dygraph
don <- xts(x = data$count, order.by = data$datetime)

# Chart
p <- dygraph(don)
p










