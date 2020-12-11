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
#remove taxa not seen more than 3 times in at least 10% of the samples. 
(bact_helli_filt = filter_taxa(bact_helli, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE))
dim(bact_helli_filt)

bact_pa <- bact_physeq %>% otu_table() %>% decostand("pa") #pa: scale x to presence/absence scale (1/0) #ensure against taking ASV with large amounts in just a few samples
bact_pa_10 <- bact_abun[(rowSums(bact_pa)/ncol(bact_pa) >= 0.10),] #select only asv present in more than 10% of samples
bact_helli_filt2 <- bact_helli[rownames(bact_pa_10),] #select asvs present in more than 10% of samples from helli transformed
dim(bact_helli_filt2)

rownames(bact_helli_filt) == rownames(bact_helli_filt2)



# micro_ps <- subset_taxa(bact_physeq, Genus == "g__Microcystis")
# micro_ps %>% otu_table()
micro_ps <- subset_taxa(bact_helli_filt, Genus == "g__Microcystis")
micro_ps %>% otu_table()

doli_ps <- subset_taxa(bact_helli_filt, Genus == "g__Dolichospermum")
doli_ps %>% otu_table() %>% row.names()

cyano_ps <- subset_taxa(bact_helli_filt, Phylum == "p__Cyanobacteria")



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

# condrare_bact_L <- rareBiosphere(bact_count_littoral)
# (174+5240+1364)/nrow(bact_count_littoral)
# 
# condrare_bact_P <- rareBiosphere(bact_count_pelagic)
# (191+5311+1275)/nrow(bact_count_pelagic)







#non rare

filt_cyanops <- filter_taxa(bact_physeq, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

#rm non-real cyano (from nico's paper): OTUs that were not prokaryotes but still present in the database
filt_cyanops %>% tax_table %>% head()
sort(unique(tax_table(filt_cyanops)[, 4]))

rmTaxa <- c("o__Cryptophyta", "o__Chlorophyta", "o__Stramenopiles") #no Streptophyta in phyloseq object
cyano_rm_ps <- subset_taxa(filt_cyanops, ta4 != "rmTaxa")




### Divnet ####
cyano_rm_ps %>% tax_table() %>% head()

#repeat for cyano
toptaxa_bac <- find.top.taxa(cyano_rm_ps,"ta7")
head(toptaxa_bac)
tt_bact <- toptaxa_bac %>% select(c("ta7"))
#see which taxa comes up most
tt_bact <- as.data.frame(table(tt_bact))
tt_bact <- tt_bact %>% arrange(desc(Freq))
head(tt_bact)

divnet_cyanosp_years <- divnet(tax_glom(cyano_rm_ps, taxrank = "ta7"),
                               base = "ASV_1",
                               X = "Years",
                               ncores = 4)
divnet_cyanosp_years

divnet_cyanosp_bloom <- divnet(tax_glom(cyano_rm_ps, taxrank = "ta7"),
                               base = "ASV_1",
                               X = "bloom2",
                               ncores = 4)
divnet_cyanosp_bloom










