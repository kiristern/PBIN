#upload cyano ASV data
cyano_counts <- read.table("cyano/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
head(cyano_counts)
cyano_taxa <- read.csv("cyano/ASVs_taxonomy_Champ_Greengenes.csv", header = T, row.names = 1, fill=T)
head(cyano_taxa)

nrow(meta)
length(cyano_counts)

cyano_taxa$Species <- rownames(cyano_taxa)

colnames(cyano_counts)
#remove X at beginning of date
colnames(cyano_counts)[1:135] <- substring(colnames(cyano_counts)[1:135], 2)
#select dates cols only
cyano_counts <- cyano_counts[1:135]

#match sample dates
meta2 <- meta

rownames(meta2)[rownames(meta2) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct

#remove sample ID at beginning
row.names(meta2) <- sub("*._*._*._*._*._*._*._","", row.names(meta2))
#change "_" to "."
rownames(meta2) <- gsub("_", ".", row.names(meta2))

#make sure meta matches cyano samples
nrow(meta2)
length(cyano_counts)

#select cols that match dates
library(tidyverse)
bact_counts <- cyano_counts[,(colnames(cyano_counts) %in% rownames(meta2))]
length(bact_counts)

meta2 <- meta2[rownames(meta2) %in% colnames(bact_counts),]
nrow(meta2)

# row.names(bact_count_meta) <- bact_count_meta$refDate
# duplicated(bact_count_meta$refDate) #check to make sure no duplicates. no dups. WHY ROWNAMES NOT WORKING?!?
# #change _ to .
# row.names(bact_count_meta) <- gsub("_", ".", row.names(bact_count_meta))
# #remove everything before 1st _ (just to keep date)
# row.names(bact_count_meta) <- gsub("^.*?\\.","", row.names(bact_count_meta))


# meta2 <- meta
# 
# #add col with the dates.only
# meta2$refDate <- dates.only
# #select dates and Site only
# meta2 <- select(meta2, refDate, Site)
# 
# meta2_P <- meta2 %>% 
#   filter(Site == "Pelagic")
# 
# meta2_L <- meta2 %>% 
#   filter(Site == "Littoral")

# #filter bacterial count by littoral / pelagic (ie. select cols that match dates from litt. or pel.)
# bact_count_littoral <- cyano_counts[, (colnames(cyano_counts) %in% meta2_L$refDate)]
# colnames(bact_count_littoral) 
# # write.table(bact_count_littoral, "bact_count_littoral.txt", row.names = T, quote = F, sep = "\t")
# 
# bact_count_pelagic <- cyano_counts[, (colnames(cyano_counts) %in% meta2_P$refDate)]
# colnames(bact_count_pelagic)
# # write.table(bact_count_pelagic, "bact_count_pelagic.txt", row.names = T, quote = F, sep = "\t")

#Phyloseq
library(phyloseq)

bac_count <- otu_table(bact_counts, taxa_are_rows = T)
cyano_taxa_ps <- tax_table(cyano_taxa)
rownames(cyano_taxa_ps) <- rownames(cyano_taxa)

#add to phyloseq object
sample_info_cyano <- sample_data(meta2)
bact_physeq <- phyloseq(bac_count, cyano_taxa_ps, sample_info_cyano)
print(bact_physeq)

#quick check
bact_physeq %>% tax_table %>% head()

bact_abun <- bact_physeq %>% otu_table()


#rename cols
colnames(tax_table(bact_physeq)) <- c("Kingdom", "Phylum", "Class",
                                      "Order", "Family", "Genus", "ASV")

library(microbiome)
bactps_helli <- transform(bact_physeq, transform = "hellinger", target = "OTU")
bact_helli <- bactps_helli %>% otu_table()

library(vegan)
#remove taxa not seen more than 10 times in at least 10% of the samples. 
bact_filt = filter_taxa(bact_physeq, function(x) sum(x > 1) > (0.10*length(x)), TRUE)


(bact_helli_filt = filter_taxa(bactps_helli, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE))
dim(bact_helli_filt)




# #### ANALYZE CONDITIONALLY RARE VIRAL ASVs & CYANOBACTERIAL ASVs ####
library(otuSummary)
condrare_viral_bact <- rareBiosphere(cyano_counts)
(nrow(condrare_viral_bact$CRT)+nrow(condrare_viral_bact$PERare)+nrow(condrare_viral_bact$otherRare))/nrow(cyano_counts)

cyno_fil <- cyano_ps_filt %>% otu_table()
condrare_viral_bact <- rareBiosphere(cyno_fil)
(nrow(condrare_viral_bact$CRT)+nrow(condrare_viral_bact$PERare)+nrow(condrare_viral_bact$otherRare))/nrow(cyano_counts)




# bact_pa <- bact_physeq %>% otu_table() %>% decostand("pa") #pa: scale x to presence/absence scale (1/0) #ensure against taking ASV with large amounts in just a few samples
# bact_pa_10 <- bact_abun[(rowSums(bact_pa)/ncol(bact_pa) >= 0.10),] #select only asv present in more than 10% of samples
# bact_helli_filt2 <- bact_helli[rownames(bact_pa_10),] #select asvs present in more than 10% of samples from helli transformed
# dim(bact_helli_filt2)
# 
# bact_helli_filt %>% otu_table() %>% rownames() == rownames(bact_helli_filt2)



# micro_ps <- subset_taxa(bact_physeq, Genus == "g__Microcystis")
# micro_ps %>% otu_table()
micro_ps_helli_filt <- subset_taxa(bact_helli_filt, Genus == "g__Microcystis")
micro_ps_helli_filt %>% otu_table()

doli_ps_helli_filt <- subset_taxa(bact_helli_filt, Genus == "g__Dolichospermum")
doli_ps_helli_filt %>% otu_table() %>% row.names()

cyano_ps_helli_filt <- subset_taxa(bact_helli_filt, Phylum == "p__Cyanobacteria")
cyano_ps_helli_filt %>% otu_table() %>% row.names()

cyano_ps_filt <- subset_taxa(bact_filt, Phylum == "p__Cyanobacteria")



# #prune OTUs that are not present in any of the samples – BUT DON’T TRIM MORE THAN THAT! it is tempting to trim noise right away, but many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate.
# BP <- prune_taxa(taxa_sums(bact_physeq) > 0, bact_physeq)


# bact_physeq_L <- phyloseq(bac_count_L, cyano_taxa_ps)
# print(bact_physeq_L)
# 
# bac_count_P <- otu_table(bact_count_pelagic, taxa_are_rows = T)
# bact_physeq_P <- phyloseq(bac_count_P, cyano_taxa_ps)
# print(bact_physeq_P)

#default graphic
plot_richness(BP)

#select alpha-diversity measure we want only:
plot_richness(BP, measure="Shannon")




# Estimate Shannon diversity and add it to the phyloseq object
sample_data(bact_physeq)$micro_ps <- get_variable(bact_physeq, "Micro.Abundance")



  
micro_ps <- subset_taxa(bact_physeq, Genus == "g__Microcystis")
micro_shannon <- estimate_richness(micro_ps, measure="Shannon")

m_meta <- micro_ps %>% sample_data
m_meta_keep <- m_meta[complete.cases(m_meta), ] #rm NAs

doli_ps <- subset_taxa(bact_physeq, Genus == "g__Dolichospermum")
d_meta <- doli_ps %>% sample_data
d_meta_keep <- d_meta[complete.cases(d_meta), ] #rm NAs






### Conditionally rare ####
condrare_bact <- rareBiosphere(cyano_counts)
(346+4613+1820)/nrow(cyano_counts)








#non rare

filt_cyanops <- filter_taxa(bact_physeq, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

#rm non-real cyano (from nico's paper): OTUs that were not prokaryotes but still present in the database
filt_cyanops %>% tax_table %>% head()
sort(unique(tax_table(filt_cyanops)[, 4]))

rmTaxa <- c("o__Cryptophyta", "o__Chlorophyta", "o__Stramenopiles") #no Streptophyta in phyloseq object
cyano_rm_ps <- subset_taxa(filt_cyanops, ta4 != "rmTaxa")





#### Divnet ####
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
toptaxa <- find.top.taxa(cyano_ps_filt,"ASV")
head(toptaxa)
tt <- toptaxa %>% select(c("ASV"))
#see which taxa comes up most
tt <- as.data.frame(table(tt))
tt <- tt %>% arrange(desc(Freq))
head(tt)

## R crashes when using bact_physeq (data too large??). reduce dim by filtering more
# check all variables in filtered phyloseq object
sample_variables(cyano_ps_filt)

div_cyan_ASV20_years_filt <- cyano_ps_filt %>%
  divnet(X = "Years", ncores = 4,
         base = "ASV_20")
div_cyan_ASV20_years_filt

div_cyan_ASV20_years_filt$shannon %>% head

#to test if alpha-diversity (by default, Shannon) is equal across the values of the covariate X:
testDiversity(div_cyan_ASV20_years_filt)

filtbactps <- reorder_ps %>% otu_table() #run lines below for re-ordered ps
filtbactps <- t(filtbactps)
#isolate for date only
df.filtbactps <- as.data.frame(row.names(filtbactps))

df.filtbactps$date <- df.filtbactps[,1]
df.filtbactps$date<- sub('(.*)[.](.*)', "\\1", df.filtbactps$date) #removes everything after last .
df.filtbactps #check which rows need to be renamed!
df.filtbactps$date[38] <- "01.06.2008"
df.filtbactps$date[41] <- "02.07.2008"
df.filtbactps$year <- as.factor(sub('.*\\.', '', df.filtbactps$date))
df.filtbactps$months <- gsub("(?:[^.]+\\.){1}([^.]+).*", "\\1", df.filtbactps$date) 

#compare the plug-in Shannon with divnet estimates
library(ggplot2)
div_cyan_ASV20_years_filt$shannon %>%
  plot(reorder_ps, color = "Years") + #run lines below to order properly!
  scale_x_discrete(labels = df.filtbactps$months, name="Month")+ #change x-axis sample name to date
  ylab("Shannon diversity estimate\n(ASV20 level)")+
  ggtitle("Shannon diversity estimate between years\n(base ASV20)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5)) #center title
#only a single DivNet estimate for each year (along with error bars). For characteristics for which many samples were observed, there are smaller error bars than for samples for which there was only one sample (seems reasonable -- we had less data).

#reorder samples by date
meta_cyano <- cyano_ps_filt %>% sample_data()
asv_cyano <- cyano_ps_filt %>% otu_table()

all(colnames(asv_cyano) %in% rownames(meta_cyano))
str(meta_cyano)
meta_cyano$Date <- as.Date(meta_cyano$Date)
meta_cyano <- meta_cyano[order(meta_cyano$Date),]
str(meta_cyano)
meta_cyano$Years <- as.factor(meta_cyano$Years)
ordered <- rownames(meta_cyano)

reorder_ps <- phyloseq(otu_table(asv_cyano, taxa_are_rows = T),
                       sample_data(meta_cyano))
otu_table(reorder_ps) <- otu_table(reorder_ps)[,ordered]


#Shannon index using breakaway
estimates <- div_cyan_ASV20_years_filt$shannon %>% summary %>% select("estimate")
ses <- sqrt(div_cyan_ASV20_years_filt$`shannon-variance`)
X <- breakaway::make_design_matrix(cyano_ps_filt, "Years")
(ba_shannon <- betta(estimates, ses, X)$table)


#distribution of Bray-Curtis distances between the samples
simplifyBeta(div_cyan_ASV20_years_filt, reorder_ps, "bray-curtis", "Years") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est,
             col = interaction(Covar1, Covar2))) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

merge_samples(div_bact_ASV11_years_filt, "Years") %>%
  sample_shannon %>%
  plot()

#Shannon index using breakaway
estimates_bact <- div_bact_ASV11_years_filt$shannon %>% summary %>% select("estimate")
ses <- sqrt(div_bact_ASV11_years_filt$`shannon-variance`)
X <- breakaway::make_design_matrix(bact_filt, "Years")
(ba_shannon <- betta(estimates_bact, ses, X)$table)






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











