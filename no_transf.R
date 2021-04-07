# No transformations: Alpha and Beta diversity, Mantel, Procruste

packageVersion("Phyloseq")
citation("phyloseq")

### Functions ###
##function to set vir_hel to same format as bact_hel
colsamp2date <- function(tab2format){
  colnames(tab2format) <- sub("*._*._*._*._*._*._*._","", colnames(tab2format))
  colnames(tab2format) <- gsub("_", ".", colnames(tab2format))
  return(tab2format)
}


##################################

library(dplyr)
library(stringr)
library(vegan)

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

# separate into pelagic and littoral phyloseq objects
vir_ps_lit <- subset_samples(viral_physeq, Site == "Littoral")
vir_ps_pel <- subset_samples(viral_physeq, Site == "Pelagic")



###### RELATIVE ABUNDANCE ######

### Top 20 ###
ps_lit_relab <- transform_sample_counts(vir_ps_lit, function(OTU) OTU/sum(OTU))

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("gmteunisse/Fantaxtic")
library("fantaxtic")
lit_relab_top20 <- get_top_taxa(ps_lit_relab, 20, relative = TRUE, discard_other = F,
             other_label = "Other")
taxa_abun_tab_lit <- psmelt(lit_relab_top20)

#create date col
library(lubridate)
m <- month(taxa_abun_tab_lit$Date)
day <- day(taxa_abun_tab_lit$Date)
md <- paste(day, m, sep="-")


#ensure colours are the same across both plots (need to run pelagic script below)
dd <- union(taxa_abun_tab_lit$species, taxa_abun_tab_pel$species)

#generate distinct colours for each asv
library("RColorBrewer")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
dd.col=sample(col_vector, length(dd))
names(dd.col) <- dd
dd.col[names(dd.col) == "Other"] <- "grey"


rel_ab_plot_lit <- taxa_abun_tab_lit %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("ASV", values = dd.col)+
  labs(x = "",
       y = "Relative Abundance",
       title = "Relative Abundance (littoral)") +
  facet_grid(~ Years, scales = "free") +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )+
  scale_x_discrete(labels = md, name="Sample date")
rel_ab_plot_lit



### repeat for pelagic ###
ps_pel_relab <- transform_sample_counts(vir_ps_pel, function(OTU) OTU/sum(OTU))

pel_relab_top20 <- get_top_taxa(ps_pel_relab, 20, relative = TRUE, discard_other = F,
                                other_label = "Other")
taxa_abun_tab_pel <- psmelt(pel_relab_top20)

m.p <- month(taxa_abun_tab_pel$Date)
day.p <- day(taxa_abun_tab_pel$Date)
md.p <- paste(day, m, sep="-")

rel_ab_plot_pel <- taxa_abun_tab_pel %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("ASV", values = dd.col)+
  labs(x = "",
       y = "Relative Abundance",
       title = "Relative Abundance (pelagic)") +
  facet_grid(~ Years, scales = "free") +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )+
  scale_x_discrete(labels = md.p, name="Sample date")
rel_ab_plot_pel




#### ALPHA DIV ####
library(breakaway)

ba.dates <- meta %>% dplyr::select(Date)

vir_ps_lit 
vir_ps_pel 

#richness by year
ba <- breakaway(vir_ps_lit)
ba

ymd <- vir_ps_lit %>% sample_data %>% get_variable("Date")
library(lubridate)
m <- month(ymd)
d <- day(ymd)
md <- paste( d, m, sep="-")

ba_vir_df = data.frame("richness" = (ba %>% summary)$estimate,
                        #"sample" = (ba %>% summary)$sample_names,
                        "error" = (ba %>% summary)$error,
                        "Years" = vir_ps_lit %>% sample_data %>% get_variable("Years"),
                        "Upper" = (ba %>% summary)$upper,
                        "Lower" = (ba %>% summary)$lower,
                       "sample"= vir_ps_lit %>% sample_data %>% get_variable("description"))
head(ba_vir_df)
ggplot(ba_vir_df, aes(x = forcats::fct_inorder(sample), y = richness, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point(size=3)+
  geom_errorbar(aes(ymin=richness-abs(richness-Lower), ymax=richness+abs(richness-Upper), width=0.05))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Breakaway richness of pelagic samples")+
  scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
  scale_y_continuous(name="Richness")


#boxplot years
ba_year = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
                     "Years" = vir_ps_lit %>% sample_data %>% get_variable("Years"))
(ba_plot <-  ggplot(ba_year, aes(x = Years, y = ba_observed_richness))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()+
  ggtitle("Observed richness by year (pelagic)")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_y_continuous(name = "Observed richness")


#boxplot years
means <- aggregate(ba_observed_richness ~ Years, ba_year, mean)
means$step <- 1:nrow(means)

#finding linear regression:
#http://r-statistics.co/Linear-Regression.html 
fit <- lm(ba_observed_richness ~ step , data = means)
coefs <- coef(fit)
summary(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

(ba_plot <-  ggplot(ba_year, aes(x = Years, y = ba_observed_richness))+
    geom_point()+
    geom_abline(intercept = coefs[1], slope = coefs[2])+
    annotate(geom="text", x = 3.5, y=620, label= paste("Adj R2 = ", r2,
                                                     "p-val = ", pval))+
    labs(title = "Observed richness by year (littoral)")+
    stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                       geom="crossbar", width=0.5)+ 
    theme_minimal())




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




### sampling depth (reads per sample)
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
summary(sample_sums(viral_physeq)) #large difference in number of reads, min=22; max=130183
sort(sample_sums(viral_physeq))

rich_depth <- data.frame("total_reads" =  phyloseq::sample_sums(viral_physeq),
                         "richness" = (ba %>% summary)$estimate)

#finding linear regression:
fit <- lm(richness ~ total_reads , data = rich_depth)
coefs <- coef(fit)
summary(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

ggplot(data = rich_depth,
       aes(x = total_reads, y = richness)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x = "\nTotal Reads", y = "Richness\n",
       title = "Observed richness by sampling depth")+
    geom_abline(intercept = coefs[1], slope = coefs[2])+
    annotate(geom="text", x = 2e+04, y=700, size = 3,
             label= paste("Adj R2 = ", r2,
                           "p-val = ", pval))



ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(viral_physeq),
                         "Years" = viral_physeq %>% sample_data %>% get_variable("Years")),
       aes(x = total_reads, y = Years)) +
  geom_bar(stat = "identity") +
  labs(x = "\nTotal Reads", y = "Year\n")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Number of reads per Year")+
  coord_flip()










#### Shannon diversity ####
vir_ps_pel
vir_ps_lit

vir_shannon <- estimate_richness(vir_ps_pel, measures="Shannon")

vir_shannon$sample <- rownames(vir_shannon)
vir_shannon$Years <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon)) #rm everything after 4th _
vir_shannon$Years <- sub(".*[/_]", "", vir_shannon$Years) #remove everything before 3rd _
# gsub("^.*\\_","", vir_shannon$Years) #another way to do same thing
vir_shannon$date <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', rownames(vir_shannon))
vir_shannon$date <- gsub("^[^_]*_", "",vir_shannon$date) #remove sample name -- just keep date
vir_shannon$date <- as.Date(vir_shannon$date, format="%d_%m_%Y")
  
head(vir_shannon)

library(lubridate)
m <- month(vir_shannon$date)
d <- day(vir_shannon$date)
md <- paste(d, m, sep="-")

ggplot(vir_shannon, aes(x = forcats::fct_inorder(sample), y = Shannon, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point(size=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Shannon diversity of pelagic samples")+
  scale_x_discrete(labels = md, name="Sample date")+ #change x-axis sample name to Month-Day
  scale_y_continuous(name="Shannon diversity")
  #facet_grid(~ Years, scales = "free")


#boxplot years
means <- aggregate(Shannon ~ Years, vir_shannon, mean)
means$step <- 1:nrow(means)

#finding linear regression:
#http://r-statistics.co/Linear-Regression.html 
fit <- lm(Shannon ~ step , data = means)
coefs <- coef(fit)
summary(fit)
names(summary(fit)) #see calls you can make

#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

#plotting with reg
(shannon_plot <-  ggplot(vir_shannon, aes(x = Years, y = Shannon))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()+
  ggtitle("Shannon diversity by year (pelagic)")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_y_continuous(name = "Shannon diversity")+
  geom_abline(intercept = fit$coefficients[1], slope =  fit$coefficients[2])+
  annotate(geom="text", x = 3, y=4.5, label= paste("Adj R2 = ", r2,
                                                  "p-val = ", pval))+
  labs(title = "Shannon diversity by year")

#### shannon vs. depth ####
shan <- estimate_richness(viral_physeq, measures="shannon")
head(shan)

viral_df = data.frame("shannon" = shan$Shannon,
                      "sample" = (viral_physeq %>% sample_data)$description,
                      "Years" = viral_physeq %>% sample_data %>% get_variable("Years"),
                      "depth" = sample_sums(viral_physeq))
head(viral_df)
str(viral_df)
ggplot(viral_df, aes(x = forcats::fct_inorder(sample), y = shannon, color = Years))+ #fct_inorder ensures plotting in order of sample date
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Shannon diversity by sample")+
  scale_x_discrete(labels = viral_physeq %>% sample_data %>% get_variable("Months"), name="Month")#change x-axis sample name to Month

### sampling depth (reads per sample)
summary(sample_sums(viral_physeq)) #large difference in number of reads, min=22; max=130183
sort(sample_sums(viral_physeq))

fit <- lm(shannon ~ depth , data = viral_df)
coefs <- coef(fit)
summary(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coefficients[2,4]))

ggplot(data = viral_df, aes(x = depth, y = shannon)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x = "\nTotal Reads", y = "Shannon\n")+
  geom_abline(intercept = fit$coefficients[1], slope =  fit$coefficients[2])+
  annotate(geom="text", x = 2e+04, y=4.75, label= paste("Adj R2 = ", r2,
                                                        "p-val = ", pval))


ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(viral_physeq),
                         "Years" = viral_physeq %>% sample_data %>% get_variable("Years")),
       aes(x = total_reads, y = Years)) +
  geom_bar(stat = "identity") +
  labs(x = "\nTotal Reads", y = "Year\n")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Number of reads per Year")+
  coord_flip()







##BETA ORDINATION###
#heatmap_distance
# dist = sqrt(phyloseq::distance(viral_physeq, "bray"))
# 
# #ordination_betadiversity_PCOA
# #PCOA need to be done with Euclidean metric distance
# pcoa=ordinate(viral_physeq, "PCoA", distance=dist)
# 
# plot_ordination(viral_physeq, dist, color  = "Years") + 
#   theme_bw() + 
#   scale_colour_manual(values = c("red","blue", "green","brown","purple","yellow","black","grey","pink", "orange")) + 
#   geom_point(size = 2) + 
#   scale_shape_manual(values=c(8, 16, 6)) + 
#   theme(axis.text.x  = element_text(vjust=0.5, size=12), 
#         axis.text.y  = element_text(vjust=0.5, size=12), 
#         axis.title.x = element_text(size = 15, face="bold", color="black"),
#         axis.title.y = element_text(size=15,face="bold",color="black"))


#### PERMANOVA ####
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
library(vegan)
package.version("vegan")
citation("vegan")
viral_physeq %>% sample_data() %>% head
jsd <- sqrt(phyloseq::distance(viral_physeq, method = "jsd")) #jsd is more robust to sampling depth
sampledf <- data.frame(sample_data(viral_physeq)) #make a df from the sample_data
adonis(jsd ~ Site, data = sampledf)

#homogeneity of dispersion test
betadisp <- betadisper(jsd, sampledf$Years)
permutest(betadisp)



### NMDS ###
nmds=metaMDS(comm = sqrt(phyloseq::distance(viral_physeq, "jsd")), k=3, trymax = 100)
# nmds <- ordinate(physeq = viral_physeq,
#                  method = "NMDS",
#                  distance = "jsd",
#                  trymax=500)
plot_ordination(physeq = viral_physeq,
                ordination = nmds,
                color = "Years",
                #shape = "Site",
                title = "NMDS of Lake Champlain viral Communities")+
  geom_point(aes(color=Years))+
  scale_color_brewer(palette = "Paired")
 # geom_point(colour="grey90", size=1.5)
nmds$stress #0.1432 good -- convergence. high stress value means that the algorithm had a hard time representing the distances between samples in 2 dimensions (anything <0.2 is considered acceptable)













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










