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
library(stringr)

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
#filter by site
meta_pelagic <- meta %>%
  filter(Site == "Pelagic")
nrow(meta_pelagic)

meta_littoral <- meta %>%
  filter(Site == "Littoral")
nrow(meta_littoral)

#filter viral count by littoral / pelagic (ie. select cols that match rownames from litt. or pel.)
vir_count_littoral <- as.data.frame(ASV_count[, (colnames(ASV_count) %in% rownames(meta_littoral))])
#write.table(vir_count_littoral, "ASV_count_littoral.txt", row.names = T, quote = F, sep = "\t")

vir_count_pelagic <- ASV_count[, (colnames(ASV_count) %in% rownames(meta_pelagic))]
#write.table(vir_count_pelagic, "ASV_count_pelagic.txt", row.names = T, quote = F, sep = "\t")

#transpose ASV count table
vir_lit <- t(vir_count_littoral)
vir_pel <- t(vir_count_pelagic)

#get total count of asv
asv_tot.l <- colSums(vir_lit)
asv_tot.p <- colSums(vir_pel)

#transform to df
asv_tot.l <- as.data.frame(asv_tot.l)
asv_tot.p <- as.data.frame(asv_tot.p)

#extract asv_id
asv_id.l <- row.names(asv_tot.l)
asv_id.p <- row.names(asv_tot.p)

#add asv_ids to df
asv_tot.l$ID=asv_id.l
asv_tot.p$ID=asv_id.p

#add new empty column
newcol <- "rel_ab"
asv_tot.l[,newcol] <- NA
asv_tot.p[,newcol] <- NA

# asv_tot <- mutate(asv_tot, sample=row.names(asv_tot))

#get relative abundance function
get_rel_abun <- function(x){
  x / sum(x)
}

#apply function to the first col of df asv_tot and put into rel_ab col of df asv_tot
asv_tot.l[3] <- get_rel_abun(asv_tot.l[1])
asv_tot.p[3] <- get_rel_abun(asv_tot.p[1])

sum(asv_tot.l$rel_ab)
sum(asv_tot.p$rel_ab)


#get top 20
topN = 20
most_abundant_taxa <- sort(taxa_sums(vir_ps_pel), TRUE)[1:topN]
most_abundant_taxa

top20.pel <- prune_taxa(names(most_abundant_taxa), vir_ps_pel)

vir_ps_lit_relab <- transform_sample_counts(top20.lit, function(x) x / sum(x))
vir_ps_pel_relab <- transform_sample_counts(top20.pel, function(x) x / sum(x))

top20.l <- top20.lit %>% otu_table()
top20.p <- top20.pel %>% otu_table()


# asv_tot.l <- arrange(asv_tot.l, desc(rel_ab))
asv_tot.p <- arrange(asv_tot.p, desc(rel_ab))
#
# top20.l <- head(asv_tot.l, 20)
# top20.l
#
top20.p <- head(asv_tot.p, 20)
# top20.p

# sum(top20.l$rel_ab)
# sum(top20.p$rel_ab)

## write.csv(top20.l, "top20L_oct29.csv")
## write.csv(top20.p, "top20P_oct29.csv")

# sort(top20.l$ID)
#
# top20L <- top20.l
# top20P <- top20.p

top20L <- read.csv("top20L_oct29.csv")
head(top20L, n=2)
top20P <- read.csv("top20P_oct29.csv")

#select top20 ASVs from full df
#asv_tax.l <- vir_lit[,(colnames(vir_lit) %in% top20.l$ID)]
asv_tax.l <- vir_lit[,(colnames(vir_lit) %in% rownames(top20.l))]

nrow(asv_tax.l)
ncol(asv_tax.l)

asv_tax.p <- vir_pel[,(colnames(vir_pel) %in% rownames(top20.p))]
nrow(asv_tax.p)
ncol(asv_tax.p)

#relative abundance matrix
asv_rel_abun.l <- decostand(asv_tax.l, method="total")
asv_rel_abun.l <- as.data.frame(asv_rel_abun.l)

asv_rel_abun.p <- decostand(asv_tax.p, method="total")
asv_rel_abun.p <- as.data.frame(asv_rel_abun.p)

#change col names to taxa ID
head(name_tax.l <- select(top20L, ID, rel_ab))
head(name_tax.p <- select(top20P, ID, rel_ab))

#if select doesn't work
name_tax.l <- top20L[, c(2, 3)]
name_tax.p <- top20P[,c(2,3)]
#add brackets around ID
name_tax.l$ID_brackets <- with(name_tax.l, paste0("(", ID, ")"))
name_tax.p$ID_brackets <- with(name_tax.p, paste0("(", ID, ")"))

#merge ASV name with ID in new col
head(virus_ID.l <- paste(top20L$Description, name_tax.l$ID_brackets, sep=" "))
head(virus_ID.p <- paste(top20P$Description, name_tax.p$ID_brackets, sep=" "))

#add col to name_tax df
name_tax.l <- tibble::add_column(name_tax.l, virus_ID.l)
name_tax.p <- add_column(name_tax.p, virus_ID.p)
head(name_tax.p)

tvir <- t(viral_physeq %>% otu_table())
nrow(tvir)

#change ASV_ to real taxonomic name
names(asv_rel_abun.l) <- name_tax.l$virus_ID.l[match(names(asv_rel_abun.l), name_tax.l$ID)]
colnames(asv_rel_abun.l)

names(asv_rel_abun.p) <- name_tax.p$virus_ID.p[match(names(asv_rel_abun.p), name_tax.p$ID)]
colnames(asv_rel_abun.p)

#duplicate each sample 20 times (number of unique ASVs)
asv_rel_abun_dup.l <- asv_rel_abun.l[rep(seq_len(nrow(asv_rel_abun.l)), each = 20), ]
asv_rel_abun_dup.p <- asv_rel_abun.p[rep(seq_len(nrow(asv_rel_abun.p)), each = 20), ]

nrow(asv_rel_abun.l)
#create new row of repeated ASV_IDs 50 times (# of samples)
ASV_ID.l <- data.frame(ASV_ID.l = c(colnames(asv_rel_abun_dup.l)))
n=94
replicate.l <- do.call("rbind", replicate(n, ASV_ID.l, simplify = FALSE))


nrow(asv_rel_abun.p)
ASV_ID.p <- data.frame(ASV_ID.p = c(colnames(asv_rel_abun_dup.p)))
n=72
replicate.p <- do.call("rbind", replicate(n, ASV_ID.p, simplify = FALSE))
#write.csv(replicate, "replicate_filt_removed.csv")
#replicate <- read.csv("replicate_test.csv")

#add ASV_ID col to asv_rel_abun df
asv_rel_abun_dup.l <- tibble::add_column(asv_rel_abun_dup.l, replicate.l, .before=1)
asv_rel_abun_dup.p <- add_column(asv_rel_abun_dup.p, replicate.p, .before=1)

#create new col named samples which duplicates rownames
samples.l <- row.names(asv_rel_abun_dup.l)
samples.p <- row.names(asv_rel_abun_dup.p)

#add col to df
asv_rel_abun_dup.l <- tibble::add_column(asv_rel_abun_dup.l, samples.l, .before=1)
asv_rel_abun_dup.p <- add_column(asv_rel_abun_dup.p, samples.p, .before=1)

#remove everything after "."
asv_rel_abun_dup.l$samples.l <- gsub("\\..*", "", asv_rel_abun_dup.l$samples.l)
asv_rel_abun_dup.p$samples.p <- gsub("\\..*", "", asv_rel_abun_dup.p$samples.p)
# write.csv(asv_rel_abun_dup, file = "repeat_sample_names_test.csv")

#transpose
sample_abundance.l <- as.data.frame(t(asv_rel_abun.l))
sample_abundance.p <- as.data.frame(t(asv_rel_abun.p))

#one col of all abundances
stacked.l <- stack(sample_abundance.l)
stacked.p <- stack(sample_abundance.p)

#add stacked to asv_rel_abun_dup df
asv_rel_abun_dup.l <- tibble::add_column(asv_rel_abun_dup.l, stacked.l$values, .after=1)
asv_rel_abun_dup.p <- add_column(asv_rel_abun_dup.p, stacked.p$values, .after=1)

#select certain cols only
df.l <- asv_rel_abun_dup.l %>% select("samples.l", "stacked.l$values")
df.p <- asv_rel_abun_dup.p %>% select("samples.p", "stacked.p$values")

#add replicate names to df
df.l <- tibble::add_column(df.l, replicate.l, .before=1)
df.p <- add_column(df.p, replicate.p, .before=1)

#rename cols
names(df.l)[names(df.l) == "stacked.l$values"] <- "abundance"
names(df.p)[names(df.p) == "stacked.p$values"] <- "abundance"

# write.csv(df, file="asvID_filt_rem_rel_abun.csv")

#make a taxa col
df.l <- df.l %>% mutate(taxa = ASV_ID.l)
df.p <- df.p %>% mutate(taxa = ASV_ID.p)

#remove brackets
df.l$taxa <- gsub("\\s*\\([^\\)]+\\)","",as.character(df.l$taxa))
df.p$taxa <- gsub("\\s*\\([^\\)]+\\)","",as.character(df.p$taxa))

#create date col
df.l$Years <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', df.l$samples.l) #rm everything after 4th _
df.l$Years <- sub(".*[/_]", "", df.l$Years) #remove everything before 3rd _
# gsub("^.*\\_","", vir_shannon$Years) #another way to do same thing
date <- unique(df.l$samples.l)
date <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', date)
date <- gsub("^[^_]*_", "",date) #remove sample name -- just keep date
date <- as.Date(date, format="%d_%m_%Y")
library(lubridate)
m <- month(date)
d <- day(date)
md <- paste(d, m, sep="-")

df.p$Years <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', df.p$samples.p) #rm everything after 4th _
df.p$Years <- sub(".*[/_]", "", df.p$Years) #remove everything before 3rd _
# gsub("^.*\\_","", vir_shannon$Years) #another way to do same thing
datep <- unique(df.p$samples.p)
datep <- sub("^([^_]*.[^_]*.[^_]*.[^_]*).*$",'\\1', datep)
datep <- gsub("^[^_]*_", "",datep) #remove sample name -- just keep date
datep <- as.Date(datep, format="%d_%m_%Y")
mp <- month(datep)
dp <- day(datep)
mdp <- paste(dp, mp, sep="-")

head(df.l)
head(df.p)


#ensure colours are the same across both plots
dd <- union(df.p$ASV_ID.p, df.p$ASV_ID.p)
dd.col <- rainbow(length(dd))
names(dd.col) <- dd

ggplot(df.p, aes(x = samples.p, y = abundance, fill = ASV_ID.p))+
  geom_bar(stat="identity", show.legend = T)+
  scale_fill_manual("ASV", values = dd.col)+
  theme_minimal()+
  facet_grid(~ Years, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5),#center title
       # legend.text = element_text(size = 5),
        legend.key.size = unit(5, 'mm'))+ 
  ggtitle("Relative abundance of top 20 ASV in pelagic zone")+
  ylab("Relative Abundance")+
  scale_x_discrete(labels = mdp, name="Sample date")

ggplot(df.p, aes(x = samples.p, y = abundance, fill = ASV_ID.p))+
  geom_bar(stat="identity", show.legend = T)+
  scale_fill_manual("ASV", values = dd.col)+
  theme_minimal()+
  facet_grid(~ Years, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5),#center title
        # legend.text = element_text(size = 5),
        legend.key.size = unit(5, 'mm'))+ 
  ggtitle("Relative abundance of top 20 ASV in littoral zone")+
  ylab("Relative Abundance")+
  scale_x_discrete(labels = mdp, name="Sample date")









#### ALPHA DIV ####
library(breakaway)

ba.dates <- meta %>% dplyr::select(Date)

vir_ps_lit 
vir_ps_pel 

#richness by year
ba <- breakaway(viral_physeq)
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
adonis(jsd ~ Years, data = sampledf)

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










