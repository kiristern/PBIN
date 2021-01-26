# No transformations: Alpha and Beta diversity, Mantel, Procruste

### Functions ###
##function to set vir_hel to same format as bact_hel
colsamp2date <- function(tab2format){
  colnames(tab2format) <- sub("*._*._*._*._*._*._*._","", colnames(tab2format))
  colnames(tab2format) <- gsub("_", ".", colnames(tab2format))
  return(tab2format)
}


##################################

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
meta$Date <- as.Date(meta$Date)
str(meta)

#order meta by date
meta <- meta[order(meta$Date),]

#get basic meta data info for methods section of report
length(ASV_count)
nrow(meta)
(amt <- meta %>% group_by(Years) %>% summarise(amount=length(Years))) #view how many samples per year
min(amt[,2])
max(amt[,2])
median(as.numeric(unlist(amt[,2]))) #get median samples per year
meta %>% group_by(Site) %>% summarise(amount=length(Site))


#ensure same samples between ASV_count and meta
asv_count<- ASV_count[,(colnames(ASV_count) %in% rownames(meta))]
length(asv_count)
nrow(meta)



#### CREATE PHYLOSEQ OBJECT ####
library(phyloseq)
library(microbiome)
#VIRAL

#add ASV count table, metadata, virTree to phyloseq table
count_phy <- otu_table(asv_count, taxa_are_rows=T)
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
viral_physeq %>% otu_table( ) %>% dim

vir_abun <- viral_physeq %>% otu_table()

#check data
print(viral_physeq)

#https://microbiome.github.io/tutorials/
summarize_phyloseq(viral_physeq)
#sparsity is how populated is the data with zeros.





#### ALPHA DIV ####
library(breakaway)

ba.dates <- meta %>% select(Date)

#richness by year
ba <- breakaway(viral_physeq)
ba

ba_vir_df = data.frame("richness" = (ba %>% summary)$estimate,
                        #"sample" = (ba %>% summary)$sample_names,
                        "error" = (ba %>% summary)$error,
                        "Years" = viral_physeq %>% sample_data %>% get_variable("Years"),
                        "Upper" = (ba %>% summary)$upper,
                        "Lower" = (ba %>% summary)$lower,
                       "sample"= viral_physeq %>% sample_data %>% get_variable("description"))
head(ba_vir_df)
ggplot(ba_vir_df, aes(x = fct_inorder(sample), y = richness, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point()+
  geom_pointrange(aes(ymin=richness-abs(richness-Lower), ymax=richness+abs(richness-Upper)))+ #vs. geom_errorbar
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Breakaway richness by sample")+
  scale_x_discrete(labels = viral_physeq %>% sample_data %>% get_variable("Months"), name="Month")#change x-axis sample name to Month



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
                      "Site" = filt_virseq %>% sample_data %>% get_variable("Site"))
ba_bloom_yn <- na.omit(ba_bloom[1:2])

(ba_plot <-  ggplot(ba_bloom, aes(x = Site, y = ba_observed_richness))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()

#t test
#http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r
group_by(ba_bloom, Site) %>%
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
with(ba_bloom, shapiro.test(ba_observed_richness[Site == "Littoral"])) # p = 0.0.04121 for bloom; p = 0.037 for littoral
# Shapiro-Wilk normality test for No Bloom
with(ba_bloom, shapiro.test(ba_observed_richness[Site == "Pelagic"])) # p = 0.02033; p = 0.64 for pelagic
#the two p-values are not greater than the significance level 0.05 implying that the distribution of the 
#data are significantly different from the normal distribution. Ie, we cannot assume the normality.
#if the data are not normally distributed, it’s recommended to use the non parametric two-samples Wilcoxon rank test.

#Wilcoxon test
# Question: Is there any significant changes in the richness of ASV during and not during bloom?
wilco <- wilcox.test(ba_observed_richness ~ Site, data = ba_bloom)
wilco

#print p-value only
wilco$p.value
# 0.01906275 < 0.05 therefore significant difference between groups




##BETA ORDINATION###
#heatmap_distance
dist = sqrt(phyloseq::distance(viral_physeq, "bray"))

#ordination_betadiversity_PCOA
#PCOA need to be done with Euclidean metric distance
pcoa=ordinate(viral_physeq, "PCoA", distance=dist)

plot_ordination(viral_physeq, dist, color  = "Years") + 
  theme_bw() + 
  scale_colour_manual(values = c("red","blue", "green","brown","purple","yellow","black","grey","pink", "orange")) + 
  geom_point(size = 2) + 
  scale_shape_manual(values=c(8, 16, 6)) + 
  theme(axis.text.x  = element_text(vjust=0.5, size=12), 
        axis.text.y  = element_text(vjust=0.5, size=12), 
        axis.title.x = element_text(size = 15, face="bold", color="black"),
        axis.title.y = element_text(size=15,face="bold",color="black"))






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



#TODO: make new for bact and viral divnets with the same samples/dates (taken from mantel)





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


