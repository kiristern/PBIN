setwd("~/Documents/GitHub/PBIN/data")

metadata <- read.csv("metadata3.csv", row.names = 1, header = T)
#create col with sample names in order to merge with second df and keep original sample name
metadata$Sample <- row.names(metadata)
metadata$OriginalSample <- row.names(metadata)
#extract data
metadata$Sample <- sub("*._*._*._*._*._*._*._","", metadata$Sample)
#change "_" to "."
metadata$Sample <- gsub("_", ".", metadata$Sample)

bact_counts <- read.table("cyano/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
bact_taxa <- read.csv("cyano/ASVs_taxonomy_Champ_Greengenes.csv", header = T, row.names = 1, fill=T)

#select cyno ASVs only
bact_asv <- bact_counts[,1:135]
#merge cyano_taxa to cyano_asv
bactdf <- merge(bact_asv, bact_taxa, by="row.names")
#remove X at beginning of date
colnames(bactdf)[2:136] <- substring(colnames(bactdf)[2:136], 2)
rownames(bactdf) <- bactdf$Row.names
bactdf <- bactdf[,-1]

bact_helli <- bactdf

bact_helli[,1:135] <- decostand(bactdf[1:135], method = "hellinger")

#### group by cyanobacteria, microcystis, dolichospermum ####
cyanobac_helli <- bact_helli[grep("Cyanobacteria", cyanodf$Phylum), ]
cyanobac_helli <- cyanobac_helli[,1:135]
colsum_cy <- as.data.frame(colSums(cyanobac_helli, na.rm=T))
#rename col
names(colsum_cy)[names(colsum_cy) == "colSums(cyanobacteria, na.rm = T)"] <- "cyano.sum.helli"
#create col with Samples in order to merge with second df
colsum_cy$Sample <- row.names(colsum_cy)

###add cyanobacteria abundance to metadata###

meta

#check similar samples
row.names(colsum_cy) %in% meta$description
#merge by similar samples
meta_cyano <- merge(meta, colsum_cy, by.y="Sample", by.x = "description", all.x = T)

#rename col
meta_cyano <- dplyr::rename(meta_cyano, cyano.sum.helli = "colSums(cyanobac_helli, na.rm = T)")
colnames(meta_cyano)
rownames(meta_cyano) <- rownames(meta)
head(meta_cyano)


###add microcystis to meta_cyano
micro.helli <- bact_helli[grep("Microcystaceae", cyanodf$Family), ]
micro.helli <- micro.helli[,1:135]
colsum_mi <- as.data.frame(colSums(micro.helli, na.rm=T))
colsum_mi$Sample <- row.names(colsum_mi)

row.names(colsum_mi) %in% meta_cyano$description
meta_cm <- merge(meta_cyano, colsum_mi, by.x="description", by.y="Sample", all.x =T)

meta_cm <- dplyr::rename(meta_cm, micro.sum.helli = "colSums(micro.helli, na.rm = T)")
colnames(meta_cm)
rownames(meta_cm) <- rownames(meta)
head(meta_cm, n=4)

###add dolichospermum to meta_cm
doli.helli <- bact_helli[grep("Dolichospermum", cyanodf$Genus), ]
doli.helli <- doli.helli[,1:135]
colsum_do <- as.data.frame(colSums(doli.helli, na.rm=T))
colsum_do$Sample <- row.names(colsum_do)

row.names(colsum_do) %in% meta_cm$description
meta_cmd <- merge(meta_cm, colsum_do, by.x = "description", by.y="Sample", all.x=T)

#put rownames back to original sample names
colnames(meta_cmd)
meta_cmd <- dplyr::rename(meta_cmd, doli.sum.helli = "colSums(doli.helli, na.rm = T)")
rownames(meta_cmd) <- rownames(meta)

meta_cmd
write.csv(meta_cmd, "meta_cmd.csv") #manually removed NA rows in OriginalSample and Sample col




