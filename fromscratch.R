library(dplyr)
library(tidyverse)

setwd("~/Documents/GitHub/PBIN/data")

#### UPLOAD DATA ####
#upload viral ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
str(ASV_count)
dim(ASV_count)
range(ASV_count)
apply(ASV_count, 2, median) #get median for each col
median(unlist(ASV_count), na.rm = T) #et median for whole df
summary(ASV_count)

colnames(ASV_count)[colnames(ASV_count) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct
head(ASV_count, n=2)
meta <- read.csv("meta_cmd.csv", row.names = 1, header = T)
head(meta, n=2)
#meta$Years <- as.factor(meta$Years)
meta$Years <- as.factor(meta$Years)
str(meta)
  
#get basic meta data info for methods section of report
length(ASV_count)
nrow(meta)
(amt <- meta %>% group_by(Years) %>% summarise(amount=length(Years))) #view how many samples per year
min(amt[,2])
max(amt[,2])
median(as.numeric(unlist(amt[,2]))) #get median samples per year
meta %>% group_by(Site) %>% summarise(amount=length(Site))


#ensure same samples between ASV_count and meta
asv_count_meta <- ASV_count[,(colnames(ASV_count) %in% rownames(meta))]
length(asv_count_meta)




#### CREATE PHYLOSEQ OBJECT ####
library(phyloseq)
library(microbiome)
#VIRAL

#add ASV count table, metadata, virTree to phyloseq table
count_phy <- otu_table(asv_count_meta, taxa_are_rows=T)
sample_info <- sample_data(meta)
virTree <- read_tree("viral_tree")

fake_taxa <- read.table("fake_viral_tax.txt", header = T, row.names = 1, fill=T)
mock_taxa <- tax_table(fake_taxa)
mock_taxa[,7] <- str_remove(mock_taxa[,7], "s__")
row.names(mock_taxa) <- mock_taxa[,7]
colnames(mock_taxa)[7] <- "species"
head(mock_taxa)
#add to phyloseq object
viral_physeq <- phyloseq(count_phy, sample_info, virTree, mock_taxa)

#check data
print(viral_physeq)

#https://microbiome.github.io/tutorials/
summarize_phyloseq(viral_physeq)
  #sparsity is how populated is the data with zeros.


#view ASV count (top 10)
sort(taxa_sums(viral_physeq), decreasing = T)[1:10]






# #### ANALYZE CONDITIONALLY RARE VIRAL ASVs & CYANOBACTERIAL ASVs ####
library(otuSummary)
condrare_viral <- rareBiosphere(vir_count)

head(condrare_viralL$summaryTable) #maximum and the minimum relative abundance, the ratio of the two abundance, the grouping of rarity, whether the taxa is singleton or doubletons
head(condrare_viralL$CRT) #subset of the above "summaryTable" only for the conditionally rare taxa.
nrow(condrare_viralL$CRT)
head(condrare_viralL$PERare) #permanently rare taxa
nrow(condrare_viralL$PERare)
head(condrare_viralL$otherRare) #summarizes rare taxa outside of CRT and PERare
nrow(condrare_viralL$otherRare)
(297+4138+892)/nrow(vir_count_littoral)







######## NON RARE ANALYSIS ###########
library(microbiome)
virps_helli <- transform(viral_physeq, transform = "hellinger", target = "OTU")
otutab <- virps_helli %>% otu_table()

library(vegan)
#the sum of an OTU across all samples is greater than 0.005% of all OTUs 
minTotRelAb = 5e-5
L = taxa_sums(otutab)
keepTaxa = taxa_names(otutab)[which((L/sum(L)) > minTotRelAb)] 
#keepTaxaL = (xL/sum(xL)) > minTotRelAb #different way to compute line above
nonrare = prune_taxa(keepTaxa, otutab)


#remove taxa not seen more than 1 times in at least 10% of the samples. This protects against ASV with small mean & trivially large coef of var
filt_virseq <- filter_taxa(viral_physeq, function(x) sum(x > 1) > (0.10*length(x)), TRUE)





#### ALPHA DIV ####
library(breakaway)
package.version("breakaway")

#richness by year
ba <- breakaway(viral_physeq)
ba
plot(ba, viral_physeq, color="Years", title="Estimate of species richness by sample date")


#boxplot years
ba_year = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
                      "Years" = viral_physeq %>% sample_data %>% get_variable("Years"))
(ba_plot <-  ggplot(ba_year, aes(x = Years, y = ba_observed_richness))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()+
  ggtitle("Observed richness by year")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_y_continuous(name = "Observed richness")





#finding linear regression:
#http://r-statistics.co/Linear-Regression.html 
fit <- lm(ba_observed_richness ~ Years, data = ba_year)
coefs <- coef(fit)
summary(fit)

#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
# (pval <- signif(summary(fit)$p.value)) #wrong! coef[2,4] != p-val

(ba_plot <-  ggplot(ba_year, aes(x = Years, y = ba_observed_richness))+
    geom_point()+
    geom_abline(intercept = coefs[1], slope = coefs[2])+
    labs(title = paste("Adj R2 = ", r2)))
                       #"Intercept =", signif(coefs[1]),
                       #"Slope =", signif(coefs[2]),
                       # "p-value =", pval)))
#geom_crossbar()
ba_plot + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                       geom="crossbar", width=0.5) + theme_minimal()




#richness by bloom/no-bloom
#boxplot bloom
ba_bloom = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
                      "Bloom" = filt_virseq %>% sample_data %>% get_variable("bloom2"))
ba_bloom_yn <- na.omit(ba_bloom[1:2])

(ba_plot <-  ggplot(ba_bloom_yn, aes(x = Bloom, y = ba_observed_richness))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()

#t test
group_by(ba_bloom_yn, Bloom) %>%
  summarise(
    # count = n(),
    mean = mean(ba_observed_richness),
    stdev = sd(ba_observed_richness)
  )

#tests to check independent t-test assumptions
#Assumption 1: are the two samples independent? Yes. not taken from the same time 
#Assumption 2: does the data from each of the 2 groups follow a normal distirbution?
#Use Shapiro-Wilk normaility test
#Null hypothesis: the data are normally distributed
#Alternative hypothesis: the data are not normally distributed

# Shapiro-Wilk normality test for Bloom
with(ba_bloom_yn, shapiro.test(ba_observed_richness[Bloom == "yes"])) # p = 0.0.04121
# Shapiro-Wilk normality test for No Bloom
with(ba_bloom_yn, shapiro.test(ba_observed_richness[Bloom == "no"])) # p = 0.02033
#the two p-values are not greater than the significance level 0.05 implying that the distribution of the 
#data are significantly different from the normal distribution. Ie, we cannot assume the normality.
#if the data are not normally distributed, it’s recommended to use the non parametric two-samples Wilcoxon rank test.

#Wilcoxon test
# Question: Is there any significant changes in the richness of ASV during and not during bloom?
wilco <- wilcox.test(ba_observed_richness ~ Bloom, data = ba_bloom_yn)
wilco

#print p-value only
wilco$p.value
# 0.01906275 < 0.05 therefore significant difference between groups





#### Divnet ####
library(magrittr)
library(DivNet)
package.version("DivNet")

#find most abundant taxa
# function to find the most abundant taxa
# Goes through a phyloseq object, picks out the most abundant taxa and gives the abundance for each
# and identifies which taxa is most abundant for which sample
find.top.taxa <- function(x,taxa){
  require(phyloseq)
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa)) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j<-apply(otu,1,which.max)
  k <- j[!duplicated(j)]
  l <- data.frame(tax[k,])
  m <- data.frame(otu[,k])
  s <- as.name(taxa)
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1,which.max)]
  m[,taxa] <- n
  return(m)
}
toptaxa <- find.top.taxa(filt_virseq,"species")
head(toptaxa)
tt <- toptaxa %>% select(c("species"))
#see which taxa comes up most
tt <- as.data.frame(table(tt))
tt <- tt %>% arrange(desc(Freq))
head(tt)

## R crashes when using viral_physeq (data too large??)
#check all variables in filtered phyloseq object
sample_variables(filt_virseq)

div_ASV1_years_filt <- filt_virseq %>%
  divnet(X = "Years", ncores = 4,
         base = "ASV_1")
div_ASV1_years_filt

div_ASV1_years_filt$shannon %>% head

#to test if alpha-diversity (by default, Shannon) is equal across the values of the covariate X:
testDiversity(div_ASV1_years_filt)

filtvirps <- filt_virseq %>% otu_table()
filtvirps <- t(filtvirps)
#isolate for date only
df.filtvirps <- as.data.frame(row.names(filtvirps))
df.filtvirps$date <- gsub("^([^_]*_[^_]*_[^_]*_[^_]*)_.*$", "\\1", df.filtvirps[,1]) #removes everything after last _
df.filtvirps$date <- sub("*._*._*._*._*._*._*._","", df.filtvirps$date) #remove sample ID at beginning

#compare the plug-in Shannon with divnet estimates
library(ggplot2)
div_ASV1_years_filt$shannon %>%
  plot(filt_virseq, color = "Years") +
  scale_x_discrete(labels = df.filtvirps$date, name="Sample date")+ #change x-axis sample name to date
  ylab("Shannon diversity estimate\n(ASV1 level)")+
  ggtitle("Shannon diversity estimate between years\n(base ASV1)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5)) #center title
#only a single DivNet estimate for each year (along with error bars). For characteristics for which many samples were observed, there are smaller error bars than for samples for which there was only one sample (seems reasonable -- we had less data).




#distribution of Bray-Curtis distances between the samples
simplifyBeta(div_ASV1_years_filt, filt_virseq, "bray-curtis", "Years") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est,
             col = interaction(Covar1, Covar2))) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

merge_samples(div_ASV1_years_filt, "Years") %>%
  sample_shannon %>%
  plot()

#Shannon index using breakaway
estimates <- div_ASV1_years$shannon %>% summary %>% select("estimate")
ses <- sqrt(div_ASV1_years$`shannon-variance`)
X <- breakaway::make_design_matrix(filt_virseq, "Years")
(ba_shannon <- betta(estimates, ses, X)$table)






#remove all na samples from bloom2
yesno <- c("yes", "no")
filt_vir_omitna_bloom <- subset_samples(filt_virseq, bloom2 %in% yesno)

div_ASV1_bloom <- filt_vir_omitna_bloom %>%
  divnet(X = "bloom2", ncores = 4,
         base = "ASV_1")
div_ASV1_bloom




div_ASV1_bloom$shannon %>%
  plot(filt_vir_omitna_bloom, color = "bloom2") +
  xlab("Sample") +
  ylab("Shannon diversity estimate\n(ASV level)")
#only a single DivNet estimate for each year (along with error bars). For characteristics for which many samples were observed, there are smaller error bars than for samples for which there was only one sample (seems reasonable -- we had less data).

#distribution of Bray-Curtis distances between the samples
simplifyBeta(div_ASV1_bloom, filt_vir_omitna_bloom, "bray-curtis", "bloom2") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est,
             col = interaction(Covar1, Covar2))) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

merge_samples(filt_virseq, "Years") %>%
  sample_shannon %>%
  plot()

#Shannon index using breakaway
estimates_b <- div_ASV1_bloom$shannon %>% summary %>% select("estimate")
ses_b <- sqrt(div_ASV1_bloom$`shannon-variance`)
X_b <- breakaway::make_design_matrix(filt_vir_omitna_bloom, "bloom2")
(ba_shannon_b <- betta(estimates_b, ses_b, X_b)$table)

#Simpson index using breakaway
estimatesb2 <- div_ASV1_bloom$simpson %>% summary %>% select("estimate")
sesb2 <- sqrt(div_ASV1_bloom$`simpson-variance`)
Xb2 <- breakaway::make_design_matrix(filt_vir_omitna_bloom, "bloom2")
(ba_simpson <- betta(estimatesb2, sesb2, Xb2)$table)









#https://bioconductor.riken.jp/packages/3.8/bioc/vignettes/microbiome/inst/doc/vignette.html
#Visually-Weighted Regression curve with smoothed error bars
# Estimate Shannon diversity and add it to the phyloseq object
sample_data(filt_virseq)$ShannonDiv <- 
metadata$ShannonDiv <- filt_virseq %>% otu_table %>% microbiome::alpha() %>% select("diversity_shannon")

#compare year and microbiome shannon diversity
microbiome::plot_regression(ShannonDiv ~ Years, metadata) #doesn't work!


#visualize the core microbiota (set of taxa that are detected in a remarkable fraction of the population above a give abundance threshold)
library(ggplot2, quiet = TRUE)
p <- plot_core(transform(filt_virseq, "compositional"), 
               plot.type = "heatmap", 
               colours = gray(seq(0,1,length=5)),
               prevalences = seq(.05, 1, .05), 
               detections = 10^seq(log10(1e-3), log10(.2), length = 10), 
               horizontal = TRUE) +
  xlab("Detection Threshold (Relative Abundance (%))") 
print(p)    


# #compositional heatmap
# tmp <- plot_composition(filt_virseq, plot.type="heatmap", transform = "hellinger")
# 
# #compositional barplot
# plot_composition(transform(filt_virseq, "compositional"), 
#                  plot.type = "barplot", sample.sort = "neatmap", label=F)






#### PERMANOVA ####
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
library(vegan)
package.version("vegan")
citation("vegan")
bc <- phyloseq::distance(filt_virseq, method = "bray") #calculate bray curtis distance matrix
sampledf <- data.frame(sample_data(viral_physeq)) #made a df from the sample_data
adonis(bc ~ Years, data = sampledf)

#homogeneity of dispersion test
betadisp <- betadisper(bc, sampledf$Years)
permutest(betadisp)










##### RDA TRANSFORMED HELLINGER ######
asv_count_meta
meta

vir_helli <- decostand(asv_count_meta, method="hellinger")
env <- meta
names(env)

#remove catagorical data from env (do RDA without sites and time -- see PERMANOVA (far) below)
str(env)
#env_vars <- env[,!(colnames(env) %in% c("Date", "Months", "Years", "Site", "Period", "bloom2"))]


env_vars <- env[,!(colnames(env) %in% c("description", "Date", "Months", "Years", "Site", "Period", "bloom2",
                                        "Microcystin", "Dolicho.Abundance", "Micro.Abundance", "Cyano.Abundance"))]
                                        #"cyano.sum.helli", "doli.sum.helli", "micro.sum.helli"))]
colnames(env_vars)

#standardize environmental data
library(vegan)
env_vars.std <- env_vars
env_vars.std[1:6] <-decostand(env_vars.std[1:6], method="standardize")
head(env_vars.std, n=2)
# apply(env_vars[,8:16], 2, mean) #data are now centered (mean~0)
# apply(env_vars[,8:16], 2, sd) #data are now scaled (st. dev =1)

#check how many rows there are without any NAs
complete.cases(env_vars.std)
#rm NAs
env_keep <- env_vars.std[complete.cases(env_vars.std), ]
env_keep %>% dplyr::glimpse() 
summary(env_keep)

#### Remove viral asvs that are not present (due to removal of NA from env vars)
sp.asv <- t(vir_helli)

#rm sample rows that are not present in env_keep
vir.rm <- sp.asv[rownames(sp.asv) %in% rownames(env_keep),]
dim(vir.rm)
#which rows are the same
#intersect(abund_name, env_row_name)

# #specific samples that are not the same
# (row_remove <- setdiff(row.names(sp.asv), row.names(env_keep)))
# #count how many are different
# length(setdiff(row.names(sp.asv), row.names(env_keep)))  
# 
# #remove rows (samples) that aren't in env_var from abundance
# sp.asv.rm <- sp.asv[!(row.names(sp.asv) %in% row_remove), ]
# dim(sp.asv.rm)

# species and enviro data without NAs needed for RDA
head(vir.rm)
head(env_keep)

length(env_keep) #number of env vars
length(vir.rm) #nnumber asv
nrow(env_keep) #number samples
nrow(vir.rm) #numer of samples

vir.rda <- rda(vir.rm~., data = env_keep)

summary(vir.rda, display = NULL)
#These results contain: (1) the proportion of variance of Y explained by the X variables (constrained proportion, 16.17% here)
#(2) the unexplained variance of Y (unconstrained proportion, 83.83% here) and 
#(3) then summarize the eigenvalues, the proportions explained and the cumulative proportion of each canonical axis 
# (each canonical axis = each constraining variable, in this case, the environmental variables from env).

# Select the significant explanatory variables by forward selection
ordiR2step(rda(vir.rm~1, data = env_keep), scope=formula(vir.rda), direction="forward", R2scope = T, pstep=1000)

anova.cca(vir.rda, by ="terms")

#get adjusted R2
(R2adj <- RsquareAdj(vir.rda)$adj.r.squared)

#test significance of the model and of each canonical axis
anova.cca(vir.rda, step=1000) #model
anova.cca(vir.rda, step=1000, by="axis") #canonical axes

#significant env vars
env.signif <- subset(env_keep, select = c("Cumulative_precipitation_t1_t7_mm", "doli.sum.helli"))
env.signif <- subset(env_keep, select = c("doli.sum.helli"))

#the proportion of variation explained by the three constraining variables being 0.058

rda.signif <- rda(vir.rm~., data=env.signif)
summary(rda.signif, display=NULL)
#The explanatory variables (Tot P and Dis P) now explain 5.8% of the variance in Y (species).

(R2adj <- RsquareAdj(rda.signif)$adj.r.squared)
#Here the strength of the relationship between X and Y corrected for the number of X variables is 0.0227.

#The significance of the model and of each canonical axis can be tested using the function anova 
#(note this is different from retaining significant variables as was done with forward selection, now we're testing the significance of the RDA axes):
anova.cca(rda.signif, step=1000)
#the RDA model is highly significant (p=0.006)
anova.cca(rda.signif, step=1000, by="axis")
# as well as both canonical axes.

#To visualize the results of the RDA, triplots can be drawn using the plot(). 
#Note that as in PCA, users can create scaling 1 and scaling 2 triplots. 
#In scaling 1, distance among objects approximate their Euclidean distances 
# in scaling 2, angles between variables X and Y reflect their correlation. 
#Thus, scaling 1 triplots can be used to interpret distances among objects 
#and scaling 2 triplots to interpret the relationships between X and Y. 

#Quick plots scaling 1 and 2
plot(rda.signif, scaling=1, main="Triplot RDA (scaling 1)")
plot(rda.signif, scaling=2, main="Triplot RDA (scaling 2)")

#Advanced plots scaling 1
plot(rda.signif, scaling=1, main="Triplot RDA - scaling 1", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(rda.signif, display="sites", choices=c(1,2), scaling=1),
       pch=21, col="black", bg="steelblue", cex=1.2)
arrows(0,0,
       scores(rda.signif, display="species", choices=c(1), scaling=1),
       scores(rda.signif, display="species", choices=c(2), scaling=1),
       col="black",length=0)
text(scores(rda.signif, display="species", choices=c(1), scaling=1),
     scores(rda.signif, display="species", choices=c(2), scaling=1),
     labels=rownames(scores(rda.signif, display="species", scaling=1)),
     col="black", cex=0.8)    
arrows(0,0,
       scores(rda.signif, display="bp", choices=c(1), scaling=1),
       scores(rda.signif, display="bp", choices=c(2), scaling=1),
       col="red")
text(scores(rda.signif, display="bp", choices=c(1), scaling=1)+0.05,
     scores(rda.signif, display="bp", choices=c(2), scaling=1)+0.05,
     labels=rownames(scores(rda.signif, display="bp", choices=c(2), scaling=1)),
     col="red", cex=1) 

#Advanced plots scaling 2
plot(rda.signif, scaling=2, main="Triplot RDA - scaling 2", type="none", xlab=c("RDA1"), ylab=c("RDA2"), xlim=c(-1,1), ylim=c(-1,1))
points(scores(rda.signif, display="sites", choices=c(1,2), scaling=2),
       pch=20, col="steelblue", cex=1.2)
arrows(0,0,
       scores(rda.signif, display="species", choices=c(1), scaling=2)*2,
       scores(rda.signif, display="species", choices=c(2), scaling=2)*2,
       col="black",length=0)
text(scores(rda.signif, display="species", choices=c(1), scaling=2)*2.1,
     scores(rda.signif, display="species", choices=c(2), scaling=2)*2.1,
     labels=rownames(scores(rda.signif, display="species", scaling=2)),
     col="black", cex=0.8)    
arrows(0,0,
       scores(rda.signif, display="bp", choices=c(1), scaling=2),
       scores(rda.signif, display="bp", choices=c(2), scaling=2),
       col="red")
text(scores(rda.signif, display="bp", choices=c(1), scaling=2)+0.05,
     scores(rda.signif, display="bp", choices=c(2), scaling=2)+0.05,
     labels=rownames(scores(rda.signif, display="bp", choices=c(2), scaling=2)),
     col="red", cex=1)







##### db-RDA #####
# https://sites.ualberta.ca/~ahamann/teaching/renr690/Lab9b.pdf
#decide which distance measure to use by looking at the rank correlations between dissimilarity indices and gradient separation (the higher the value the better)
rankindex(env_keep, sp.asv.rm, indices = c("euc", "man", "gow", "bra", "kul"), stepacross = F, method = "spearman")

#remove site and temporal data from env_keep (do RDA without sites and time -- see PERMANOVA (far) below)
colnames(env_keep)
env_vars <- env_keep[,!(colnames(env_keep) %in% c("Months", "Years", "Site", "Period", "bloom2"))]
colnames(env_vars)
  
dbRDA <- capscale(vir.rm ~ ., data = env_keep, distance="bray")
plot(dbRDA)


#get the R2 for model fit for teh constrained ordinations
dbR2 <- RsquareAdj(dbRDA)$r.squared
dbR2

#adjusted R^2: adjusted R2 measures the unbiased amount of explained variation
dbR2adj <- RsquareAdj(dbRDA)$adj.r.squared
dbR2adj

# plot the RDA using ggplot (ggord package)
library(ggord)
ggord(dbRDA, env_keep$bloom2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  


# overall test of the significance of the analysis
anova(dbRDA)
#test axes for significance
anova(dbRDA, by="axis", perm.max=500)
#test for sig. environmental variables
anova(dbRDA, by="terms", permu=200)
#                                   Df SumOfSqs      F Pr(>F)    
# Months                             6   3.7920 5.0998  0.001 ***
# Years                              7   5.8803 6.7786  0.001 ***
# Site                               1   0.1032 0.8324  0.630    
# Period                             2   0.4735 1.9104  0.005 ** 
# bloom2                             1   0.1635 1.3195  0.170    
# Total_Phosphorus_ug                1   0.2037 1.6436  0.068 .  
# Total_Nitrogen_mg                  1   0.1820 1.4683  0.111    
# Dissolved_P                        1   0.1254 1.0122  0.414    
# Dissolved_N                        1   0.2277 1.8374  0.035 *  
# Cumulative_precipitation_t1_t7_mm  1   0.1861 1.5018  0.086 .  
# Mean_temperature_t0_t7             1   0.2035 1.6418  0.056 .  
# Dolicho.Abundance                  1   0.0633 0.5108  0.968    
# Micro.Abundance                    1   0.0736 0.5936  0.916    
# Cyano.Abundance                    1   0.1260 1.0170  0.437    
# Residual                          23   2.8503

#transforming negtaive eigenvalues (get rid of the -ve eigenvalues when doing the analysis)
#1) add a constant:
dbRDA_constant <- capscale(sp.asv.rm ~., env_vars, distance = "bray", add=TRUE)
plot(dbRDA_constant)
anova(dbRDA_constant)

#2) take the sqrt of dissimilarities
dbRDA_sqrt <- capscale(sp.asv.rm ~., env_vars, dist="bray", sqrt.dist=T)
plot(dbRDA_sqrt)
anova(dbRDA_sqrt)

#3) Do a square root transformation, Wisconsin double standardization (this emphasizes the environmental variables):
dbRDA_metaMDS=capscale(sp.asv.rm ~ ., env_vars, dist="bray", metaMDS=TRUE, sqrt.dist=TRUE) 
plot(dbRDA_metaMDS)
anova(dbRDA_metaMDS)



#select the significant explanatory variables
?ordiR2step
ordiR2step(rda(sp.asv.rm~1, data=env_vars), scope=formula(dbRDA), direction= "forward", R2scope=TRUE, pstep=1000)

#the proportion of variation explained by the constraining variables being 0.3542

#retain significant variables only
signif_env <- subset(env_vars, select = c("Total_Nitrogen_mg", "Mean_temperature_t0_t7", "Dissolved_N"))

#rda with significant variables 
dbrda_env_signif <- capscale(sp.asv.rm~., data=signif_env, dist="bray", sqrt.dist=T)
summary(dbrda_env_signif, display=NULL)
screeplot(dbrda_env_signif)
#The proportion of the variance of Y (species) explained by the X (env) variables = 37.93% (constrained) 
#the unexplained variance of Y = 62.07% (unconstrained)

#adjusted R^2
(adjR2 <- RsquareAdj(dbrda_env_signif)$adj.r.squared) 
#strength of the relationship between X and Y corrected for the number of X variables is 27.94%
#the explanatory variables explain 38.32% of the variance in Y (species)

#test significance of the RDA canonical axis (ie. constraining variables, ie. env var)
anova.cca(dbrda_env_signif, step=1000) 
#RDA model is highly significant (p=0.001)
anova.cca(dbrda_env_signif, step=1000, by="axis") 
#axis' 1-2 v significant(p=0.001), 3 (p=0.026)

#plot dbRDA
#Quick plots scaling 1 and 2
plot(dbrda_env_signif, scaling=2, main="Triplot RDA (scaling 2)")

#with ellipses
ggord(dbrda_env_signif, signif_env$bloom2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  







#https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html
#simple RDA
simpleRDA <- rda(vir.rm ~., env_keep)
summary(simpleRDA)
screeplot(simpleRDA)
  # The total variance of the data set, partitioned into constrained and unconstrained variances, is a standard result. This result shows how much variation in your response variables was redundant with the variation in your explanatory variables. If the constrained variance is much higher than your unconstrained variance, the analysis suggests that much of the variation in the response data may be accounted for by your explanatory variables. If, however, there is a large proportion of unconstrained variation (i.e. variation in your response matrix that is non-redundant with the variation in the explanatory matrix), then the results should be interpreted with caution as only a small amount of the variation in your response matrix is displayed
  # constrained axes (RDA axes)
  # unconstrained axes (PCA axes)
  #Object and response variable scores are often reported as "site" and "species" scores, respectively. These scores are the coordinates used to ordinate points and vectors. The coordinates of variables should be understood as the "tip" of their vector with the origin as its "tail". The direction of the vector is the direction of increase for that variable.
  #Explanatory variable scores, also referred to as constraining variable scores, may be interpreted as response variable scores when the explanatory variable in question is quantitative. Scores for each state of nominal or factorial variables are the coordinates of these states' centroids and show the average position of the sites that have that state.
# first 2 axes explain most of the variation (RDA1: 0.3452 & RDA2: 0.2721) so plotting by these two axes represent the data well
  # the unconstrained eigenvalue (PC1) is 0.09707, which is comparatively small, which means it does not display any important residenual structure of teh response data
  # each RDA axis has an eigenvalue, which is the proportion of the variance explained by each axis
  # The species and site scores = where the sites and species fall along the axes
      # Partition of variance: the overall variance is partitioned into constrained and unconstrained fractions.
        # constrained: the amount of variance the species by site matrix is explained by the explanatory variables (expressed as a proportion--is equivalent to R2 in a multiole regression). This R2 is biased so you have to look at the adjusted R2
      # Eigenvalues and their contribute to variance: RDA1 -> RDA26 for the canonical axes, and unconstrained axes. The cumulative contribute of the variance is the proportion of the total variance of the response data explained by the RDA. The results also give the eigenvalues.
      # Eigenvalues: the canonical eigenvalues are decreasing in value (in the order they are presented). But sometimes the residual structure (residual eigenvalue PC1) can be larger than the last RDA eigenvalue, which means that the residual structure has more variance than some of the structures that can be explained by the explan variables.
      # Canonical eigenvalues: measure the amount of variance explained by the RDA model
      # Residual eigenvalues: measure the amount of variance represented by the residual axes
      # Accumulated constrained eigenvalues: cumulative amounts of variance expressed as proportions of the total explained variance
      # Species scores: coordinates of the tips of vectors representing the response variables in bi or triplots
      # Site scores (weighted sums of species scores): coordinates of the sites as expressed in the space of the response variables
      # Site constraints (linear combinations of constraining variables): coordinates of the sites in space of the explanatory variables
      # Biplot scores for constraining variables: coordinates of the tips of the vectors represnting explanatory variables
      # Centroids for factor constraints: coordinates of centroids of levels of factor variables

#getting the canonical coefficient (ie. the equivalent of regression coefficients for each explanatory variable on each canonical axis)
coef(simpleRDA)

#get the R2 for model fit for teh constrained ordinations
R2 <- RsquareAdj(simpleRDA)$r.squared
R2

#adjusted R^2: adjusted R2 measures the unbiased amount of explained variation
R2adj <- RsquareAdj(simpleRDA)$adj.r.squared
R2adj #0.5084 ==> this model explains 50.8% of the variation in the data (if used the biased R2, any variable included in the explanatory responses would increase the R2, so the R2 needs to be adjusted for the number of eplanatory variables -- esp since we have 8)

#plot RDA
  # two ways: using "wa" (weighted sums of species) and "lc" (fitted site scores)
    # wa: more robust to noise in the environmental variables but are a step between constrained towards unconstrained. Default because more robust to noise in data
    # lc: orthogonal linear combinations of the explanatory variable

#The interpretation of these plots depends on what scaling has been chosen. 
#In general, consider type I scaling if the distances between objects are of particular value or if most explanatory variables are binary or nominal.
# Consider type II scaling if the correlative relationships between variables are of more interest.

# https://mb3is.megx.net/gustame/constrained-analyses/rda
### Type I Scaling - Distance plots (object focused) ###
#Distances between object points approximate Euclidean distances. Thus, objects ordinated closer together can be expected to have similar variable values. This will not always hold true, as RDA only recovers part of the variation in the data set.
#Right-angled projections of object points onto vectors representing response variables approximate variable values for a given object.
#The angles between vectors representing response variables are meaningless.
#The angles between vectors representing response variables and those representing explanatory variables reflect their (linear) correlation.
#Note, that binary explanatory variables may be represented as points. These points are the centroids of objects which have a state "1" for a given binary variable. Projecting centroid points onto a vector representing a response variable reflects the relationship between these variables.
#Distances between centroids and between centroids and object points approximates Euclidean distances.

# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1
plot(simpleRDA, scaling=1, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")

# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(simpleRDA, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')

### Type II Scaling - Correlation plots (response variable focused) ###
#Distances between object points should not be considered to approximate Euclidean distances.
#Right-angled projections of object points onto vectors representing response variables approximate variable values for a given object.
#The angles between all vectors reflect their (linear) correlation. The correlation is equal to the cosine of the angle between vectors (e.g. a vector pair describing an angle of 90° are uncorrelated as cos(90) = 0), those describing an angle of 20° have strong, positive correlation as cos(20) =  0.94)
# Note, that binary or nominal explanatory variables may be represented as points. These points are the centroids of objects which have a state "1" for a given binary variable or realise a particular level of a nominal explanatory variable. Projecting centroid points onto a vector representing a response variable reflects the relationship between these variables.

# Scaling 2
plot(simpleRDA, main="Triplot RDA matrix ~ env - scaling 2 - wa scores")
spe2.sc <- scores(simpleRDA, choices=1:2, display="sp") # scores() choices= indicates which axes are to be selected, make sure to specify the scaling if its different than 2 
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')


# plot the RDA using ggplot (ggord package)
library(ggord)
ggord(simpleRDA, env_vars$bloom2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  


#plot results of RDA using "lc"
# site scores as linear combinations of the environmental variables 
# Scaling 1
plot(simpleRDA, scaling=1, display=c("lc", "sp", "cn"), main="Triplot RDA matrix ~ env -scaling 1- lc scores")
arrows(0,0, spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')

# scaling 2
plot(simpleRDA, display=c("sp", "lc", "cn"), main="Triplot RDA matrix ~ env -scaling2-lc scores", col=env_vars$bloom2)
arrows(0,0,spe2.sc[,1],spe2.sc[,2], length=0, lty=1,col='red')


# To choose the elements that are plotted, use the argument display=c(), sp=species, wa= site scores in the species space (weighted averages), lc= fitted site scores (linear combinations of explanatory variables) and cn= constraints (the explanatory variables)

## Fwd selection of variables ##
# variance inflation factors in the RDA
vif.cca(simpleRDA) #Anything above 10 should be examined or avoided

# RDA with all explanatory variables  
spe.rda.all <- rda(sp.asv.rm ~ Years + Site + Period + bloom2 + Total_Phosphorus_ug + Total_Nitrogen_mg + Dissolved_P + Dissolved_N + 
                     Cumulative_precipitation_t1_t7_mm + Mean_temperature_t0_t7 + Micro.Abundance, data=env_vars)

# Forward selection using ordistep
ordistep(rda(sp.asv.rm ~., data=env_vars), direction="forward", pstep=1000, R2scop=TRUE) #R2scope only accepts models with lower adjusted R2

# db-RDA
dbRDA2 <- capscale(sp.asv.rm ~ ., data=env_vars, distance = "bray")
dbRDA2
plot(dbRDA2)
summary(dbRDA2)
screeplot(dbRDA2)
#why use CAP: it allows for any dissimilarity measure & takes into account any correlation structure among the response variables (i.e. this is why its called ‘canonical’)








gen.imp <- t(vir_abund_helli)
env <- meta

#RDA requires complete data frames (i.e., no missing data)
sum(is.na(gen.imp))
#look at the structure of the df
str(env) 
#confirm that asv and enviro data are in the same order
identical(rownames(gen.imp), rownames(env))

#RDA is a regression-based method, and so can be subject to problems when using highly correlated predictors (Dormann et al., 2013). Generally, the |r| > 0.7 “rule of thumb” is a good guideline for removing correlated predictors. We will also check for multicollinearity using Variance Inflation Factors (VIF)
#Variable reduction should be guided by an ecological interpretation of the relevance of possible predictors. Here, we use the function pairs.panels to visualize correlations among our predictors. Correlation coefficients are in the upper right diagonal, with their size scaled to their |r|. The lower left shows scatter plots, while the diagonal shows histograms of the data
library(psych)
pairs.panels(env[,7:15], scale=T)
#remove total phsophorus and dissolved p
pred <- subset(env, select=-c(Dissolved_P, Total_Phosphorus_ug))

#rda
vir.rda <- rda(formula=gen.imp ~ Months + Years + Site + Period + bloom2 +
                 Total_Nitrogen_mg  + Dissolved_N + Cumulative_precipitation_t1_t7_mm + Mean_temperature_t0_t7+Cyano.Abundance, data = pred, scale=T) 





# *no RDA* for site and time. use PERMANOVA (combining adonis2 and betadisper)











##### Cross Correlation function (CCF) ######
## calculate shannon diversity for cyanobacteria, microcystis, and dolichospermum, and phage, using the same samples.

#



shan.div.cyano <- filt_virseq %>%
  divnet(X = "Years", ncores = 4,
         base = "ASV_1")
div_ASV1_years


#compare the plug-in Shannon with divnet estimates
library(ggplot2)
div_ASV1_years$shannon %>%
  plot(filt_virseq, color = "Years") +
  xlab("Sample") +
  ylab("Shannon diversity estimate\n(ASV level)")
#only a single DivNet estimate for each year (along with error bars). For characteristics for which many samples were observed, there are smaller error bars than for samples for which there was only one sample (seems reasonable -- we had less data).






#####  GLMM #####
#https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html
pkgs_CRAN <- c("lme4","MCMCglmm","blme",
               "pbkrtest","coda","aods3","bbmle","ggplot2",
               "reshape2","plyr","numDeriv","Hmisc",
               "plotMCMC","gridExtra","R2admb",
               "broom.mixed","dotwhisker")
install.packages(pkgs_CRAN)
rr <- "http://www.math.mcmaster.ca/bolker/R"
install.packages("glmmADMB",type="source",repos=rr)
library("devtools")

## primary GLMM-fitting packages:
library("lme4")
library("glmmADMB")      ## (not on CRAN)
library("glmmTMB")
library("MCMCglmm")
library("blme")
library("MASS")          ## for glmmPQL (base R)
library("nlme")          ## for intervals(), tundra example (base R)
## auxiliary
library("ggplot2")       ## for pretty plots generally
## ggplot customization:
theme_set(theme_bw())
scale_colour_discrete <- function(...,palette="Set1") {
  scale_colour_brewer(...,palette=palette)
}
scale_colour_orig <- ggplot2::scale_colour_discrete
scale_fill_discrete <- function(...,palette="Set1") {
  scale_fill_brewer(...,palette=palette)
}
## to squash facets together ...
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))
library("gridExtra")     ## for grid.arrange()
library("broom.mixed")
## n.b. as of 25 Sep 2018, need bbolker github version of dotwhisker ...
library("dotwhisker")
library("coda")      ## MCMC diagnostics
library("aods3")     ## overdispersion diagnostics
library("plotMCMC") ## pretty plots from MCMC fits  
library("bbmle")     ## AICtab
library("pbkrtest")  ## parametric bootstrap
library("Hmisc")
## for general-purpose reshaping and data manipulation:
library("reshape2")
library("plyr")
## for illustrating effects of observation-level variance in binary data:
library("numDeriv")
library("glmmADMB")






















##### MRT #####
library(mvpart)
library(plyr)

sp.asv.rm
colnames(env_keep)

#use all data to identify which group accurately predicts 
mrt <- mvpart(data.matrix(sp.asv.rm) ~ ., env_keep, xv = "pick", xvmult=1000)

mrt_month <- mvpart(sp.asv.rm~ Months, env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.923

mrt_year <- mvpart(as.matrix(sp.asv.rm)~ Years, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.825

mrt_site <- mvpart(as.matrix(sp.asv.rm)~ Site, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.991


mrt_period <- mvpart(as.matrix(sp.asv.rm)~ Period, env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.943


mrt_bloom <- mvpart(as.matrix(sp.asv.rm)~ bloom2, env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.965


mrt_totP <- mvpart(as.matrix(sp.asv.rm)~ Total_Phosphorus_ug, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.955


mrt_totN <- mvpart(as.matrix(sp.asv.rm)~ Total_Nitrogen_mg, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.865


mrt_DisP <- mvpart(as.matrix(sp.asv.rm)~ Dissolved_P, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.957
(DisP_R2 <- RsquareAdj(mrt_DisP)$adj.r.squared)
rpart.pca(mrt_DisP)

mrt_DisN <- mvpart(as.matrix(sp.asv.rm)~ Dissolved_N, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.925


mrt_precip <- mvpart(as.matrix(sp.asv.rm)~ Cumulative_precipitation_t1_t7_mm, env_keep,
                      legend=T, margin=0.01, cp=0, xv="pick",
                      xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.947


mrt_temp <- mvpart(as.matrix(sp.asv.rm)~ Mean_temperature_t0_t7, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.751


mrt_cyano <- mvpart(as.matrix(sp.asv.rm)~ Cyano.Abundance, env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.946

mrt_doli <- mvpart(as.matrix(sp.asv.rm)~ Dolicho.Abundance, env_keep,
                    legend=T, margin=0.01, cp=0, xv="pick",
                    xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
# error = 0.923


mrt_micro <- mvpart(as.matrix(sp.asv.rm)~ Micro.Abundance, env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)
#error = 0.916


#significant vars based off of MRT scores
mrt_signif <- mvpart(as.matrix(sp.asv.rm) ~ Years+Total_Nitrogen_mg+Mean_temperature_t0_t7, env_keep,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)

# significant vars based off of ordi2step
signif_env
mrt_ordi <- mvpart(as.matrix(sp.asv.rm) ~ Years+Months+Dissolved_N+Period+Dissolved_P, signif_env,
                     legend=T, margin=0.01, cp=0, xv="pick",
                     xval=nrow(sp.asv.rm), xvmult=100, which=4, big.pts=T, bars=F)



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














#### ML ####
#RF
library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 

rf_samp <- ASV_count

#upload cyano ASV data
cyano_counts <- read.table("cyano/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
colnames(cyano_counts)
bact_meta <- cyano_counts[1:135] #select only first 135 cols
# duplicated(colnames(rf_meta)) #check if there are duplicates
colnames(bact_meta) <- substring(colnames(bact_meta), 2) #remove X at beginning of the date


### ENSURE SAME SAMPLES FOR BOTH BACTERIAL AND CYANO
#match sample dates
sampID <- colnames(rf_samp)
#remove sample ID at beginning
sampID_date <- sub("*._*._*._*._*._*._*._","", sampID)
#change "_" to "."
(sampDate <- gsub("_", ".", sampID_date))

# ie. select cols that match dates
library(tidyverse)
resp_var <- bact_meta[(colnames(bact_meta) %in% sampDate),]

colnames(rf_samp)[colnames(rf_samp) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct
colnames(rf_samp) <- sub("*._*._*._*._*._*._*._","", colnames(rf_samp)) #remove everything before 1st _ (just to keep date)
colnames(rf_samp) <- gsub("_", ".", colnames(rf_samp)) #change _ to .

length(resp_var)
length(rf_samp)

pred_var <- rf_samp[,(colnames(rf_samp) %in% colnames(resp_var))]
colnames(pred_var)
colnames(resp_var)

colnames(pred_var) %in% colnames(resp_var)
length(pred_var)
length(resp_var)
colnames(resp_var) %in% colnames(pred_var)
resp_var <- resp_var[,(colnames(resp_var) %in% colnames(pred_var))]

length(pred_var)
length(resp_var)

pred_var
resp_var

resp_var_transposed <- t(resp_var)
resp_var_transposed <- as.data.frame(resp_var_transposed)
rownames(resp_var_transposed) %in% colnames(pred_var)

ps_predict <- otu_table(pred_var, taxa_are_rows=T)
ps_response <- sample_data(resp_var_transposed)

#make phyloseq object
RF_ps <- merge_phyloseq(ps_predict, ps_response)

#hellinger transform species dataset: gives low weight to rare species
RF_ps_helli <- transform(RF_ps, transform = "hellinger", target = "OTU")
RF_ps_helli %>% otu_table() %>% head

# filtered such that only OTUs with a sum greater than 10^-5 in at least 10% of samples are kept 
(RF_ps_helli_filt = filter_taxa(RF_ps_helli, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE))

#make a df of training data with OTUs as cols and samples as rows
predictors <- as.data.frame(t(otu_table(RF_ps_helli_filt)))
dim(predictors)

# Make one col for outcome/response variable
# response <- as.factor(sample_data(RF_ps))

#bacterial phyloseq objects (run fromscratch_cyano.R)
micro_ps
#micro_ps %>% tax_table()
doli_ps
cyano_ps



#make one column for outcome/response variable
micro_response <- micro_ps %>% otu_table() %>% t()

# micro_df <- data.frame()
# for (i in 1:ncol(micro_response)){
#   micro_df[,i] <- merge(micro_response[,i], predictors, by="row.names")
# }

# m7 <- merge(micro_response[,1], predictors, by="row.names")
# rownames(m7) <- m7[,1]
# m7$Row.names <- NULL
# names(m7)[1] <- "micro_ASV_7"
# 
# m8 <- merge(micro_response[,2], predictors, by="row.names")
# rownames(m8) <- m8[,1]
# m8$Row.names <- NULL
# names(m8)[1] <- "micro_ASV_8"
# 
# m143 <- merge(micro_response[,3], predictors, by="row.names")
# rownames(m143) <- m143[,1]
# m143$Row.names <- NULL
# names(m143)[1] <- "micro_ASV_143"

# rownames(micro_response) == rownames(m7)
# head(rownames(micro_response))
# head(rownames(m7))
# 
micro_response <- micro_response[order(row.names(micro_response)), ]
# m7 <- m7[order(row.names(m7)),]
# 
# m8 <- m8[order(row.names(m8)),]
# m143 <- m143[order(row.names(m143)),]


#put in the same order for rownames
rownames(predictors) %in% rownames(micro_response)
rownames(predictors) == rownames(micro_response)

predictors <- predictors[order(row.names(predictors)),]

rownames(predictors) == rownames(micro_response)

head(micro_response)
#make one column for outcome/response variable
micro_ASV_7 <- as.double(micro_response[,1])
micro_ASV_8 <- as.double(micro_response[,2])
micro_ASV_143 <- as.double(micro_response[,3])
str(micro_ASV_143)

#combine into one df
m7 <- data.frame(micro_ASV_7, predictors)
m8 <- data.frame(micro_ASV_8, predictors)
m143 <- data.frame(micro_ASV_143, predictors)

m143
#write.csv(m143, "micro143.csv")
head(colnames(m143))
str(m143)

# names(m7)[1] <- "micro_ASV_7"
# names(m8)[1] <- "micro_ASV_8"
# names(m143)[1] <- "micro_ASV_143"

# # doli
# doli_response <- doli_ps %>% otu_table() %>% t()
# head(doli_response)
# 
# d10 <- merge(doli_response[,1], predictors, by="row.names")
# rownames(d10) <- d10[,1]
# d10$Row.names <- NULL
# names(d10)[1] <- "doli_ASV_10"
# 
# d11 <- merge(doli_response[,2], predictors, by="row.names")
# rownames(d11) <- d11[,1]
# d11$Row.names <- NULL
# names(d11)[1] <- "doli_ASV_11"
# 
# d16 <- merge(doli_response[,3], predictors, by="row.names")
# rownames(d16) <- d16[,1]
# d16$Row.names <- NULL
# names(d16)[1] <- "doli_ASV_16"
# 
# d18 <- merge(doli_response[,4], predictors, by="row.names")
# rownames(d18) <- d18[,1]
# d18$Row.names <- NULL
# names(d18)[1] <- "doli_ASV_18"
# 
# d50 <- merge(doli_response[,5], predictors, by="row.names")
# rownames(d50) <- d50[,1]
# d50$Row.names <- NULL
# names(d50)[1] <- "doli_ASV_50"
# 
# d85 <- merge(doli_response[,6], predictors, by="row.names")
# rownames(d85) <- d85[,1]
# d85$Row.names <- NULL
# names(d85)[1] <- "doli_ASV_85"
# 
# d110 <- merge(doli_response[,7], predictors, by="row.names")
# rownames(d110) <- d110[,1]
# d110$Row.names <- NULL
# names(d110)[1] <- "doli_ASV_110"
# 
# d129 <- merge(doli_response[,8], predictors, by="row.names")
# rownames(d129) <- d129[,1]
# d129$Row.names <- NULL
# names(d129)[1] <- "doli_ASV_129"
# 
# d356 <- merge(doli_response[,9], predictors, by="row.names")
# rownames(d356) <- d356[,1]
# d356$Row.names <- NULL
# names(d356)[1] <- "doli_ASV_356"
# 
# d683 <- merge(doli_response[,10], predictors, by="row.names")
# rownames(d683) <- d683[,1]
# d683$Row.names <- NULL
# names(d683)[1] <- "doli_ASV_683"
# 
# d735 <- merge(doli_response[,11], predictors, by="row.names")
# rownames(d735) <- d735[,1]
# d735$Row.names <- NULL
# names(d735)[1] <- "doli_ASV_735"


# plot 20 most important viral ASV
plottop20 <- function(RFm){
  imp <- varImp(RFm)
  imp <- data.frame(imp$importance)
  imp <- data.frame(predictors = rownames(imp), imp)
  #order the predictor levels by importance
  imp.sort <- arrange(imp, desc(Overall))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  #select top 10 predictors
  imp.top20 <- imp.sort[1:20,]
  
  #ggplot
  imp.top20 %>%
    arrange(Overall) %>%
    mutate(predictors=factor(predictors, levels=predictors)) %>%
    ggplot(
      aes(x=predictors, y = Overall)) +
    geom_bar(stat = "identity", fill = "indianred") +
    coord_flip()+
    #ggtitle("Most important viral ASV for predicting\nMicrocystis ASV7") +
    theme(plot.title = element_text(hjust = 0.5)) -> plottop20
  return(plottop20)
}

top20 <- function(RFm){
  imp <- varImp(RFm)
  imp <- data.frame(imp$importance)
  imp <- data.frame(predictors = rownames(imp), imp)
  #order the predictor levels by importance
  imp.sort <- arrange(imp, desc(Overall))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  #select top 10 predictors
  imp.top20 <- imp.sort[1:20,]
  return(imp.top20)
}



library(randomForest)
library(caret)
library(e1071)

#tune the algorithm
dataset <- m143
dim(dataset)
samp <- createDataPartition(dataset$micro_ASV_143, p = 0.9, list = F)
head(dataset)[,1:5]
training <- dataset[samp,]
testing <- dataset[-samp,]

length(dataset)
x = training[,2:577]
y = training[,1]


# Random Search
fitControl <- trainControl(method="repeatedcv", 
                           number = 10,
                           repeats = 10, 
                           search="random",
                           selectionFunction = "oneSE")

#random generate 15 mtry values with tunelength=15
rf_random <- train(micro_ASV_143~., #predictor data object
                   data=training, #outcome data object
                   method="rf", #specifies type of model: rf = random forest
                   metric="RMSE",
                   tuneLength=15,
                   trControl = fitControl,
                   importance=T)  
rf_random$results #if metric isn't specified, will see all results
print(rf_random)
# RMSE was used to select the optimal model using  the one SE rule.
# The final value used for the model was mtry = 59
plot(rf_random)
plottop20(rf_random)
rf.rand20<- top20(rf_random)
rf.rand20

#tuneGrid
#Create control function for training with 10 (number) folds and keep 3 (repeats) folds for training. search method is grid.
fitControl <- trainControl(method="repeatedcv", 
                           number = 10,
                           repeats = 10, 
                           search="grid")

#create tunegrid with 10 values from 55:65 for mtry to tunning model. Our train function will change number of entry variable at each split according to tunegrid. 
tunegrid <- expand.grid(.mtry = (55:75)) 

rf_grid <- train(micro_ASV_143~., #predictor data object
                   data=training, #outcome data object
                   method="rf", #specifies type of model: rf = random forest
                   metric="RMSE",
                   trControl = fitControl,
                   tuneGrid = tunegrid,
                   importance=T)#,
# returnResamp="all",
# savePredictions="all")
rf_grid$results
print(rf_grid)
# mtry 23 is the optimal
plot(rf_grid)
plottop20(rf_grid)
rf.grid20 <- top20(rf_grid)


#tuneRF
bestMtry <- tuneRF(x, y, stepFactor = 1.5, improve = 1e-5, ntree = 500)
print(bestMtry)


#custom tunning with multiple parameters: mtry and ntree
customRF <- list(type = "Regression",
                 library = "randomForest",
                 loop = NULL)

customRF$parameters <- data.frame(parameter = c("mtry", "ntree"),
                                  class = rep("numeric", 2),
                                  label = c("mtry", "ntree"))

customRF$grid <- function(x, y, len = NULL, search = "grid") {}

customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
  randomForest(x, y,
               mtry = param$mtry,
               ntree=param$ntree)
}

#Predict label
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)

#Predict prob
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")

customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

library(doParallel)
cores <- makeCluster(detectCores()-1)
registerDoParallel(cores = cores)
start_time <- Sys.time() #start timer

# train model
control <- trainControl(method="repeatedcv", 
                        number=10, 
                        repeats=3,
                        allowParallel = TRUE)

tunegrid <- expand.grid(.mtry=c(68:75),.ntree=c(100,500, 1000))

custom <- train(micro_ASV_143~., data=training, 
                method=customRF, 
                metric="RMSE", 
                tuneGrid=tunegrid, 
                trControl=control)

end_time <- Sys.time() #end timer
end_time - start_time # Display time

summary(custom)
plot(custom)
stopCluster(cores)


#RF with selected parameters
caret.rf <- train(micro_ASV_143~.,
                  data=training,
                  method="rf",
                  ntree=1000,
                  metric="RMSE",
                  tuneGrid=data.frame(mtry=71),
                  importance=T,
                  trControl=trainControl(method="repeatedcv"))
print(caret.rf)
plottop20(caret.rf)
(rf20 <- top20(caret.rf))



#PLS
pls.mod <- train(micro_ASV_143~.,
                 data=training,
                 method="pls",
                 scale=T,
                 trControl=trainControl("repeatedcv", 
                                        number=10,
                                        repeats = 10,
                                        selectionFunction = "oneSE"),
                 importance=T,
                 tuneLength=15
  
)
plot(pls.mod)
print(pls.mod)
pls.mod$bestTune
top20(pls.mod)
plottop20(pls.mod)
summary(pls.mod$finalModel)

# #random search in a specified grid:
# control <- trainControl("repeatedcv", 
#                         repeats = 10,
#                         number = 10,
#                         selectionFunction = "oneSE",
#                         search="grid")
# ncompf <- min(0.75*nrow(training) - 1, 0.75*nrow(training))
# grid <- expand.grid(ncomp=seq(0,ncompf,1))
# samp <- sample(1:nrow(grid), 45)
# ncompT <- data.frame(.ncomp = sort(c(samp)))
# fit.pls <- train(micro_ASV_143 ~. , 
#                  data=training, 
#                  method="pls", 
#                  metric="RMSE", 
#                  tuneGrid=ncompT, 
#                  tuneLength = 20,
#                  returnResamp = "all",
#                  savePredictions = "all",
#                  trControl=control,
#                  importance=T)
# plot(fit.pls)
# fit.pls$bestTune
# #plot(varImp(fit.pls), 20, main = "PLS")
# plottop20(fit.pls)
# (pls20 <- top20(fit.pls))



#svmLinear
svmlinear <- train(micro_ASV_143~.,
                   data=training,
                   method="svmLinear",
                   trControl=trainControl(method = "repeatedcv",
                                          number = 10,
                                          repeats=10),
                   tuneGrid = expand.grid(C = seq(0,2, length=15)),
                   importnace=T)
print(svmlinear)
svmlinear$bestTune
plot(svmlinear)
top20(svmlinear)
plottop20(svmlinear)

#Create control function for training with 10 (number) folds and keep 3 (repeats) folds for training. search method is grid.
fitControl <- trainControl(method="repeatedcv",
                           number = 10,
                           repeats = 10,
                           selectionFunction = "oneSE",
                           search="grid")
grid <- expand.grid(C = seq(0,ncompf,1))
svmlin <- train(micro_ASV_143~., #predictor data object
                  data=training, #outcome data object
                  method="svmLinear", #specifies type of model: rf = random forest
                  tuneLength=10,
                  tuneGrid=grid,
                  importance=T,
                  returnResamp="all",
                  savePredictions="all",
                  trControl = fitControl)
svmlin$results
print(svmlin)
plot(svmlin)
plottop20(svmlin)
svmlin20 <- top20(svmlin)
# Print the best tuning parameter C that maximizes model accuracy
svmlin$bestTune



#svmRadial
#Create control function for training with 10 (number) folds and keep 3 (repeats) folds for training. search method is grid.
fitControl <- trainControl(method="repeatedcv",
                           number = 10,
                           repeats = 3,
                           savePred=T, 
                           classProb=T)

svmrad <- train(micro_ASV_143~., #predictor data object
                data=training, #outcome data object
                method="svmRadial", #specifies type of model: rf = random forest
                trControl = fitControl,
                preProcess = c("center", "scale"),
                tuneLength=10,
                importance=T)#,
# returnResamp="all",
# savePredictions="all")
svmrad$results
print(svmrad)
# C = 16 is the optimal
plot(svmrad)
plottop20(svmrad)
svmrad20 <- top20(svmrad)
# Print the best tuning parameter C that maximizes model accuracy
svmrad$bestTune




#knn
fitControl <- trainControl(method="repeatedcv",
                           number = 10,
                           repeats = 3,
                           savePred=T, 
                           classProb=T)

KNN <- train(micro_ASV_143~., #predictor data object
                data=training, #outcome data object
                method="knn", #specifies type of model: rf = random forest
                trControl = fitControl,
                preProcess = c("center", "scale"),
                tuneLength=10,
                importance=T)#,
# returnResamp="all",
# savePredictions="all")
KNN$results
print(KNN)
# K = 5 is the optimal
plot(KNN)
plottop20(KNN)
knn20 <- top20(KNN)
# Print the best tuning parameter C that maximizes model accuracy
KNN$bestTune


comp.mods <- data.frame(rf.rand20, rf.grid20, svmlin20, svmrad20, knn20)

pred_var_helli <- pred_var %>% decostand(method = "hellinger")

asv.for.m143<- pred_var_helli[c("ASV_310", 
                          "ASV_180", 
                          "ASV_153",
                          "ASV_416",
                          "ASV_654",
                          "ASV_275", 
                          "ASV_139",
                          "ASV_715",
                          "ASV_633",
                          "ASV_625",
                          "ASV_93",
                          "ASV_4"),]
head(asv.for.m143 <- t(asv.for.m143))

head(respM143 <- micro_response[,3])

head(timeseriedf <- merge(respM143, asv.for.m143, by="row.names"))
rownames(timeseriedf) <- timeseriedf[,1] #set col1 as rownames
orgDate <- timeseriedf[,1]
orgDate2 <- sub("^(.*)[.].*", "\\1", orgDate[-c(1,3)]) #remove everything after last period. 1st and 3rd entries don't have same dims so omit
orgDate[-c(1,3)] <- orgDate2 
timeseriedf[,1] <- orgDate

# try to change as date doesn't work!
# timeseriedf$Row.names <- gsub(".", "-", timeseriedf$Row.names) #change _ to .
# timeseriedf[order(as.Date(timeseriedf$Row.names, format="%d/%m/%Y"))]

#break up into own cols
for (i in 1:nrow(timeseriedf)){
  timeseriedf$day[i] <- str_extract_all(timeseriedf$Row.names, "[^.]+")[[i]][[1]]
  timeseriedf$month[i] <- str_extract_all(timeseriedf$Row.names, "[^.]+")[[i]][[2]]
  timeseriedf$year[i] <- str_extract_all(timeseriedf$Row.names, "[^.]+")[[i]][[3]]
}

timeseriedf[,1] <- NULL #remove col1 can also call "timeseriesdf$Row.names"

# library(lubridate)
# timeseriedf$date <- timeseriedf %>%
#   select(day, month, year) %>%
#   mutate(date = make_date(year, month, day))

timeseriedf$date <- as.Date(with(timeseriedf, paste(year, month, day, sep="-")), "%Y-%m-%d")
timeseriedf <- timeseriedf[order(as.Date(timeseriedf$date, format="%Y-%m-%d")),] #order by date
ts.df <- timeseriedf[,1:13] #keep only ASV cols
ts <- log(ts.df)

#plot
ts %>% 
  rownames_to_column() %>% 
  gather(key = key, value = value, ASV_143:ASV_4) %>% 
  mutate(rowname = factor(rowname)) %>% 
  ggplot(aes(as.numeric(rowname), value, color = key)) + 
  geom_point() + 
  geom_line() +
  scale_x_continuous(name = "order by date") +
  scale_y_continuous(name = "log(abondance)")+
  #scale_color_manual(values=c("black", "red", "blue"))+
  theme_bw()

library(gghighlight)
ts %>% 
  rownames_to_column() %>% 
  gather(key = key, value = value, ASV_143:ASV_4) %>% 
  mutate(rowname = factor(rowname)) %>% 
  gghighlight(aes(as.numeric(rowname), value, color = key)) + 
  geom_point() + 
  geom_line() +
  scale_x_continuous(name = "order by date") +
  scale_y_continuous(name = "log(abondance)")+
  #scale_color_manual(values=c("black", "red", "blue"))+
  theme_bw()



bestMtry <- tuneRF(x,y, stepFactor = 1.5, improve = 1e-5, ntree = 500)
print(bestMtry)



str(training)

library(caretEnsemble)
library(pls)
library(kernlab)
control <- trainControl(method="repeatedcv", 
                        repeats = 10,
                        number = 10,
                        selectionFunction = "oneSE",
                        returnResamp = "all", 
                        savePredictions = "final")
model_list <- caretList(micro_ASV_143 ~ .,
                        data=training,
                        methodList = c("pls","svmRadial", "rf", "knn"),
                        tuneList = NULL,
                        continue_on_fail = FALSE, 
                        tuneLength=10, 
                        returnResamp = "final", 
                        savePredictions = "all",  
                        importance=T,
                        trControl = control)
results <- resamples(model_list)
summary(results)
dotplot(results)
modelCor(resamples(model_list))
options(digits = 3)
model_results <- data.frame(
  PLS=min(model_list$pls$results$RMSE),
  SVM = min(model_list$svmRadial$results$RMSE),
  RF = min(model_list$rf$results$RMSE),
  KNN = min(model_list$knn$results$RMSE)
)
print(model_results)
options(digits = 3)
model_results2 <- data.frame(
  PLS=min(model_list$pls$results$Rsquared),
  SVM = min(model_list$svmRadial$results$Rsquared),
  RF = min(model_list$rf$results$Rsquared),
  KNN = min(model_list$knn$results$Rsquared)
)
print(model_results2)
set.seed(222)
ensemble_1 <- caretEnsemble(model_list,
                            metric = "RMSE",
                            trControl = control)
summary(ensemble_1)







### Graph temporal series ###
head(respM7 <- micro_response[,3])
impVir <- pred_var[c("ASV_855", "ASV_4", "ASV_243"),]
head(impVir <- t(impVir))

head(timeseriedf <- merge(respM7, impVir, by="row.names"))
rownames(timeseriedf) <- timeseriedf[,1] #set col1 as rownames
orgDate <- timeseriedf[,1]
orgDate2 <- sub("^(.*)[.].*", "\\1", orgDate[-c(1,3)]) #remove everything after last period. 1st and 3rd entries don't have same dims so omit
orgDate[-c(1,3)] <- orgDate2 
timeseriedf[,1] <- orgDate

# try to change as date doesn't work!
# timeseriedf$Row.names <- gsub(".", "-", timeseriedf$Row.names) #change _ to .
# timeseriedf[order(as.Date(timeseriedf$Row.names, format="%d/%m/%Y"))]

#break up into own cols
for (i in 1:nrow(timeseriedf)){
  timeseriedf$day[i] <- str_extract_all(timeseriedf$Row.names, "[^.]+")[[i]][[1]]
  timeseriedf$month[i] <- str_extract_all(timeseriedf$Row.names, "[^.]+")[[i]][[2]]
  timeseriedf$year[i] <- str_extract_all(timeseriedf$Row.names, "[^.]+")[[i]][[3]]
}

timeseriedf[,1] <- NULL #remove col1 can also call "timeseriesdf$Row.names"

# library(lubridate)
# timeseriedf$date <- timeseriedf %>%
#   select(day, month, year) %>%
#   mutate(date = make_date(year, month, day))

timeseriedf$date <- as.Date(with(timeseriedf, paste(year, month, day, sep="-")), "%Y-%m-%d")
timeseriedf <- timeseriedf[order(as.Date(timeseriedf$date, format="%Y-%m-%d")),] #order by date
timeseriedf <- timeseriedf[,1:3] #keep only ASV cols
ts <- log(timeseriedf)

#plot
ts %>% 
  rownames_to_column() %>% 
  gather(key = key, value = value, ASV_143:ASV_4) %>% 
  mutate(rowname = factor(rowname)) %>% 
  ggplot(aes(as.numeric(rowname), value, color = key)) + 
  geom_point() + 
  geom_line() +
  scale_x_continuous(name = "order by date") +
  scale_y_continuous(name = "log(abondance)")+
  scale_color_manual(values=c("black", "red", "blue"))+
  theme_bw()












