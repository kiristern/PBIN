setwd("~/Documents/GitHub/PBIN/data")

metadata <- read.csv("metadata3.csv", row.names = 1, header = T)
#create col with sample names in order to merge with second df and keep original sample name
metadata$Sample <- row.names(metadata)
metadata$OriginalSample <- row.names(metadata)
#extract data
metadata$Sample <- sub("*._*._*._*._*._*._*._","", metadata$Sample)
#change "_" to "."
metadata$Sample <- gsub("_", ".", metadata$Sample)

cyano_counts <- read.table("cyano/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
cyano_taxa <- read.csv("cyano/ASVs_taxonomy_Champ_Greengenes.csv", header = T, row.names = 1, fill=T)

#select cyno ASVs only
cyano_asv <- cyano_counts[,1:135]
#merge cyano_taxa to cyano_asv
cyanodf <- merge(cyano_asv, cyano_taxa, by="row.names")
#remove X at beginning of date
colnames(cyanodf)[1:136] <- substring(colnames(cyanodf)[1:136], 2)

#### group by cyanobacteria, microcystis, dolichospermum ####
cyanobacteria <- cyanodf[grep("Cyanobacteria", cyanodf$Phylum), ]
cyanobacteria <- cyanobacteria[,1:136]
colsum_cy <- as.data.frame(colSums(cyanobacteria[,-1], na.rm=T))
#rename col
names(colsum_cy)[names(colsum_cy) == "colSums(cyanobacteria[, -1], na.rm = T)"] <- "Cyano Abundance"
#create col with Samples in order to merge with second df
colsum_cy$Sample <- row.names(colsum_cy)

###add cyanobacteria abundance to metadata###
#check similar samples
row.names(colsum_cy) %in% metadata$Sample
#merge by similar samples
meta_cyano <- merge(colsum_cy, metadata, by="Sample", all=T)

###add microcystis to meta_cyano
microcystis <- cyanodf[grep("Microcystaceae", cyanodf$Family), ]
microcystis <- microcystis[,1:136]
colsum_mi <- as.data.frame(colSums(microcystis[,-1], na.rm=T))
names(colsum_mi)[names(colsum_mi) == "colSums(microcystis[, -1], na.rm = T)"] <- "Micro Abundance"
colsum_mi$Sample <- row.names(colsum_mi)

row.names(colsum_mi) %in% meta_cyano$Sample
meta_cm <- merge(colsum_mi, meta_cyano, by="Sample", all=T)

###add dolichospermum to meta_cm
dolichospermum <- cyanodf[grep("Dolichospermum", cyanodf$Genus), ]
dolichospermum <- dolichospermum[,1:136]
colsum_do <- as.data.frame(colSums(dolichospermum[,-1], na.rm=T))
names(colsum_do)[names(colsum_do) == "colSums(dolichospermum[, -1], na.rm = T)"] <- "Dolicho Abundance"
colsum_do$Sample <- row.names(colsum_do)

row.names(colsum_do) %in% meta_cm$Sample
meta_cmd <- merge(colsum_do, meta_cm, by="Sample", all=T)

#put rownames back to original sample names
colnames(meta_cmd)
meta <- meta_cmd[,-17]
row.names(meta) <- meta_cmd[,17]

#write.csv(meta_cmd, "meta_cmd2.csv") #manually removed NA rows in OriginalSample and Sample col




