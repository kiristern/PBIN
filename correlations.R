### Compare alpha diversities of cyano and phage ###
dim(ba_cyano_df)
dim(ba_vir_df)
head(ba_cyano_df)
head(ba_vir_df)

ba_shannon
ba_shannon$sample <- ba_cyano_df$sample
vir_shannon
vir_shannon$sample <- ba_vir_df$sample

#merge dfs to ensure same dims
df <- merge(ba_cyano_df,ba_vir_df, by="sample")

df <- merge(ba_shannon,vir_shannon, by="sample")
dim(df)
head(df)

# #check that ba_cyano_df is .x
# which(ba_cyano_df$sample == "01.06.2008") 
# ba_cyano_df[37,]

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
hist(alphadiv.corr$cyan.rich)
hist(alphadiv.corr$vir.rich)

#Visual inspection of the data normality 
ggqqplot(alphadiv.corr$cyan.rich, ylab = "cyan.rich")
ggqqplot(alphadiv.corr$vir.rich, ylab = "vir.rich")

richness.corr <- cor.test(alphadiv.corr$cyan.rich, alphadiv.corr$vir.rich, method = "spearman", exact = F)
richness.corr
# p-val = 0.02084 < 0.05, therefore cyano rich and vir rich are signif correlated with a corr coeff (rho) of -0.2076







#### Hellinger Transformed & filtered ####
#dm_ps <- subset_taxa(bact_physeq, Genus == "g__Dolichospermum" | Genus=="g__Microcystis")
bact_helli <- transform(bact_physeq, transform = "hellinger", target = "OTU")
bact_helli_filt = filter_taxa(bact_helli, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)
b_helli_filt <- bact_helli_filt %>% otu_table()

cya_helli_filt <- subset_taxa(bact_helli_filt, !Phylum == "p__Cyanobacteria")
c_helli_filt <- cya_helli_filt %>% otu_table()

#change colID to match
virps_helli_filt
sample_names(virps_helli_filt) <- sample_data(virps_helli_filt)$description
vir_helli_filt <- virps_helli_filt %>% otu_table()

tc_helli_filt <- t(c_helli_filt) #changed to b_helli_filt to see all bacteria corr
dim(tc_helli_filt)
tvir_helli_filt <- t(vir_helli_filt)
dim(tvir_helli_filt)

#set to same dims
vir_hel_filt<- tvir_helli_filt[rownames(tvir_helli_filt) %in% rownames(tc_helli_filt),]
dim(vir_hel_filt)
colnames(vir_hel_filt) <- paste0("vir_", colnames(vir_hel_filt))
write.csv(vir_hel_filt, "vir_hel_filt.csv")
vir_hel_filt <- read.csv("vir_hel_filt.csv", header = T, row.names = 1)
dim(vir_hel_filt)

cyano_hel_filt <- tc_helli_filt[rownames(tc_helli_filt) %in% rownames(tvir_helli_filt),]
dim(cyano_hel_filt)
write.csv(cyano_hel_filt, "cyano_hel_filt.csv")
cyano_hel_filt <- read.csv("cyano_hel_filt.csv", header=T, row.names = 1)
dim(cyano_hel_filt)

#https://microbiome.github.io/tutorials/Heatmap.html
library(microbiome)
#helli_filt_corr <- microbiome::associate(vir_hel_filt, cyano_hel_filt, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
#write.csv(helli_filt_corr, "correlation_bactnoCyan-vir.csv")
corr_bactnoCyan <- read.csv("correlation_bactnoCyan-vir.csv")
corr_cyan <- read.csv("correlation_cyan-vir.csv")
head(corr_bactnoCyan)
head(corr_cyan)

#helli_filt_corr <- helli_filt_corr[order(helli_filt_corr[,"X1"], helli_filt_corr[,"X2"]),]

#calculate JSD dissimilarity matrix on full dataset
data.dist <- vegdist()

# #order the rows and cols with levels argument if needed:
# helli_filt_corr <- helli_filt_corr %>% order(helli_filt_corr$X1)
# helli_filt_corr$X2 <- factor(helli_filt_corr$X2, levels=unique(as.character(helli_filt_corr$X2)))

#set b&w theme
library(ggplot2)
# theme_set(theme_bw())
# #pick only the correlations with p<0.05
sig_bactnoCyano <- filter(corr_bactnoCyan, p.adj < 0.05)
sig_corrCyano <- filter(corr_cyan, p.adj < 0.05)

# #arrange fig
subtab %>% ggplot(aes(x = X1, y = X2, fill = Correlation))+
  geom_tile()+
  scale_fill_gradientn("Correlation",
                       breaks = seq(from=-1, to=1, by = 0.2),
                       colours = c("darkblue", "blue", "white", "red", "darkred"),
                       limits=c(-1,1))+
  theme(axis.text.x = element_text(angle=90))+
  xlab("Viral") + ylab("Bacterial")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

corr.filt <- helli_filt_corr %>% filter(Correlation > 0.6)
corr.filt[order(-corr.filt$Correlation),] #sort by descending Correlation strength

unique(corr.filt$X1)
unique(corr.filt$X2)

heat(corr.filt, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05)

filt_bactnoCyan <- sig_bactnoCyano %>% filter(Correlation > 0.6)
filt_corrCyan <- sig_corrCyano %>% filter(Correlation > 0.6)

#see which phage is the same between both dfs
idx <- which(filt_corrCyan$X1 %in% filt_bactnoCyan$X1)
unique(filt_corrCyan$X1[idx])



which(spe.top17$from %in% corr.filt$X1)
which(spe.top17$to %in% corr.filt$X2)





###### Mantel test #####
# https://www.flutterbys.com.au/stats/tut/tut15.2.html
library(vegan)

#FILTER OUT LOW ABUNDANCE.
(bact_filt_ps = filter_taxa(bact_physeq, function(x) sum(x > 10) > (0.10*length(x)), TRUE))
(cyano_filt_ps <- subset_taxa(bact_filt_ps, Phylum == "p__Cyanobacteria"))
(vir_filt_ps = filter_taxa(viral_physeq, function(x) sum(x > 1) > (0.10*length(x)), TRUE))

#make sure same sample length
bac.ps <- prune_samples(rownames(sample_data(bact_filt_ps)) %in% rownames(meta2), bact_filt_ps)

cyan.ps <- prune_samples(rownames(sample_data(cyano_filt_ps)) %in% rownames(meta2), cyano_filt_ps)
vir.ps <- prune_samples(sample_data(vir_filt_ps)$description %in% rownames(meta2), vir_filt_ps)

sample_names(vir.ps) <- sample_data(vir.ps)$description

dist_vir<-sqrt(phyloseq::distance(vir.ps, "jsd"))
dist_cyan<-sqrt(phyloseq::distance(cyan.ps, "jsd"))
dist.bact <- sqrt(phyloseq::distance(bac.ps, "jsd"))

plot(dist_vir, dist_cyan)
abline(lm(dist_vir ~ dist_cyan))

(cyano.mantel <-mantel(dist_vir, dist_cyan, method = "spearman", perm=1000))
mantel(dist_vir, dist.bact, method="spearman", perm=1000)

#plot
hist(cyano.mantel$perm)
abline(v=cyano.mantel$statistic)


# #generate correlogram (multivariate correlation plot)
# plot(data.dist, env.dist, type="n")
# points(data.dist, env.dist, pch=20)
# # axis(1 )
# # axis(2, las=1)
# mtext("Viral distances", 1, line = 3)
# mtext("Bacterial distances", 2, line=3)
# abline(lm(data.dist ~ env.dist))



##### Procrustes #####
protest(dist_vir,dist_cyan)
protest(dist_vir, dist.bact)









