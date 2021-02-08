### Compare alpha diversities of cyano and phage ###
dim(ba_bact_df)
dim(ba_vir_df)
head(ba_bact_df)
head(ba_vir_df)

ba_shannon
ba_shannon$sample <- ba_bact_df$sample
vir_shannon
vir_shannon$sample <- ba_vir_df$sample

#merge dfs to ensure same dims
df <- merge(ba_bact_df,ba_vir_df, by="sample")

df <- merge(ba_shannon,vir_shannon, by="sample")
dim(df)
head(df)

# #check that ba_bact_df is .x
# which(ba_bact_df$sample == "01.06.2008") 
# ba_bact_df[37,]

head(alpha_bact <- df[,c("sample", "richness.x")])
head(alpha_vir <- df[,c("sample", "richness.y")])

head(alpha_bact <- df[,c("sample", "Shannon.x")])
head(alpha_vir <- df[,c("sample", "Shannon.y")])

#cross correlation and lagged regressions
alpha_corr <- ccf(alpha_bact$richness.x, alpha_vir$richness.y)

alpha_corr <- ccf(alpha_bact$Shannon.x, alpha_vir$Shannon.y)
alpha_corr

alphadiv.corr <- df[,c("richness.x", "richness.y")]

alphadiv.corr <- df[,c("Shannon.x", "Shannon.y")]
rownames(alphadiv.corr) <- df$sample

names(alphadiv.corr)[names(alphadiv.corr) == "richness.x"] <- "cyan.rich"
names(alphadiv.corr)[names(alphadiv.corr) == "richness.y"] <- "vir.rich"

names(alphadiv.corr)[names(alphadiv.corr) == "Shannon.x"] <- "cyan.rich"
names(alphadiv.corr)[names(alphadiv.corr) == "Shannon.y"] <- "vir.rich"
head(alphadiv.corr)


#http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
library("ggpubr")

### Is the covariation linear? 
ggscatter(alphadiv.corr, y = "cyan.rich", x = "vir.rich", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab = "Cyanobacterial richness", xlab = "Viral richness", title = "Correlation between viral and bacterial Shannon diversity")
#Yes, from the plot the relationship is linear. 
#situations where the scatter plots show curved patterns, dealing with nonlinear association b/w the 2 variables.

### Are the data from each of the 2 variables (x, y) follow a normal distribution?

##Shapiro-Wilk test can be performed as follow:
#Null hypothesis: the data are normally distributed
#Alternative hypothesis: the data are not normally distributed

# Shapiro-Wilk normality test for bact
shapiro.test(alphadiv.corr$cyan.rich) # => p = 0.034
# Shapiro-Wilk normality test for vir
shapiro.test(alphadiv.corr$vir.rich) # => p = 0.016
#p-values <0.05 implying that the distribution of the data are significantly different from normal distribution. 
#In other words, we cannot assume the normality.

#Visual inspection of the data normality 
ggqqplot(alphadiv.corr$cyan.rich, ylab = "cyan.rich")
ggqqplot(alphadiv.corr$vir.rich, ylab = "vir.rich")

richness.corr <- cor.test(alphadiv.corr$cyan.rich, alphadiv.corr$vir.rich, method = "spearman")
richness.corr
# p-val = 0.02084 < 0.05, therefore bact rich and vir rich are signif correlated with a corr coeff (rho) of -0.2076






#Correlation Species-Env
spec <- abund_clean_env_hel[,2:68]
env <- abund_clean_env_hel[,71:86]
# Define data sets to cross-correlate
x <- as.matrix(spec)
y <- as.matrix(env)
# Cross correlate data sets
correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
# Or, alternatively, the same output is also available in a handy table format
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
kable(head(correlation.table))


#### Hellinger Transformed & filtered ####
#transformed viral and cyano table with hellinger:
tcyano_helli_filt <- t(cyano_helli_filt)
dim(tcyano_helli_filt)
tvir_helli_filt <- t(vir_helli_filt)
dim(tvir_helli_filt)

#set to same dims
vir_hel_filt<- tvir_helli_filt[rownames(tvir_helli_filt) %in% rownames(tcyano_helli_filt),]
dim(vir_hel_filt)
#vir_hel_filt <- read.csv("vir_hel_filt.csv", header = T, row.names = 1)
colnames(vir_hel_filt) <- paste0("vir_", colnames(vir_hel_filt))

cyano_hel_filt <- tcyano_helli_filt[rownames(tcyano_helli_filt) %in% rownames(tvir_helli_filt),]
dim(cyano_hel_filt)
#cyano_hel_filt <- read.csv("cyano_hel_filt.csv", header=T, row.names = 1)

#https://microbiome.github.io/tutorials/Heatmap.html
helli_filt_corr <- associate(vir_hel_filt, cyano_hel_filt, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
head(helli_filt_corr)
heat(helli_filt_corr)

corr.filt <- helli_filt_corr %>% filter(Correlation > 0.6)
corr.filt[order(corr.filt$X2),]

unique(corr.filt$X1)
unique(corr.filt$X2)


















###### Mantel test #####
# https://www.flutterbys.com.au/stats/tut/tut15.2.html
## Need to run fromscrach_cyano.R


### RM RARE ####
phage <- filt_vir
dim(phage)
colnames(phage)

#filter bact to rm rare
print(bact_physeq)
filt_bact <- filter_taxa(bact_physeq, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

bact <- filt_bact %>% otu_table()
dim(bact)
colnames(bact)

#remove sample ID at beginning
colnames(phage) <- sub("*._*._*._*._*._*._*._","", colnames(phage))
colnames(phage) <- gsub("_", ".", colnames(phage))


#select cols that match dates
# bact_keep <- bact[,(colnames(bact) %in% colnames(phage))]
# dim(bact_keep)
phage_keep <- phage[,(colnames(phage) %in% colnames(bact))]
dim(phage_keep)

tbact_keep <- t(bact)
tphage_keep <- t(phage_keep)
sum(is.na(tphage_keep))

dist_vir<-sqrt(vegdist(tphage_keep, method = "bray"))
dist_bac<-sqrt(vegdist(tbact_keep, method = "bray"))

plot(dist_vir, dist_bac)
abline(lm(dist_vir ~ dist_bac))

bact.mantel <- mantel(dist_vir, dist_bac, perm=1000)
bact.mantel


#plot
hist(bact.mantel$perm)
abline(v=bact.mantel$statistic)


# #generate correlogram (multivariate correlation plot)
# plot(data.dist, env.dist, type="n")
# points(data.dist, env.dist, pch=20)
# # axis(1 )
# # axis(2, las=1)
# mtext("Viral distances", 1, line = 3)
# mtext("Bacterial distances", 2, line=3)
# abline(lm(data.dist ~ env.dist))



#### MANTEL FOR CYANO ####
print(bact_physeq)
cyano_ps <- subset_taxa(bact_physeq, Phylum == "p__Cyanobacteria")

filt_cyano_ps <- filter_taxa(cyano_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
filt_cyano <- filt_cyano_ps %>% otu_table() 

cyno <- filt_cyano
dim(cyno)


# #select cols that match dates
# cyano_keep <- cyno[,(colnames(cyno) %in% colnames(phage))]
# dim(cyano_keep)

tcyano_keep <- t(cyno)

dim(tphage_keep)
dim(tcyano_keep)

# dist_vir<-sqrt(vegdist(phage_keep, method = "bray"))
dist_cyano<-sqrt(vegdist(tcyano_keep, method = "bray"))

plot(dist_vir, dist_cyano)
abline(lm(dist_vir ~ dist_cyano))

cyano.mantel <- mantel(dist_vir, dist_cyano, perm=1000)
cyano.mantel

#plot
hist(cyano.mantel$perm)
abline(v=cyano.mantel$statistic)






##### Procrustes #####

protest(dist_vir,dist_cyano)
protest(dist_vir,dist_bac)

