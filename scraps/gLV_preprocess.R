
### run Initialize.R first ####
#removed taxa not seen more than 3 times in at least 5% of the samples
remove_rare_vir <- t(filt_vir)

bacterial <- read.csv("data/cyano/bacteria_metagenonet.csv", row.names=1, header=T)
bact_taxa <- read.csv("data/cyano/ASVs_taxonomy_Champ_Greengenes.csv")

#### merge taxa ID to ASVs ####
#assign ASV_X as rownames in bact_taxa df
bact_taxa_rownames <- bact_taxa[,-1] 
rownames(bact_taxa_rownames) <- bact_taxa[,1]
#merge df by rownames
bact <- merge(bacterial, bact_taxa_rownames, by="row.names")

#### group by microcystis, dolichospermum ####

microcystis <- bact[grep("Microcystaceae", bact$Family), ]
microcystis <- microcystis[1:136]
micro <- microcystis[,-1]
rownames(micro) <- microcystis[,1]
micro <- t(micro)
dim(micro)

dolichospermum <- bact[grep("Dolichospermum", bact$Genus), ]
dolichospermum <- dolichospermum[1:136]
doli <- dolichospermum[,-1]
rownames(doli) <- dolichospermum[,1]
doli <- t(doli)
dim(doli)


#remove X in front of rownames
row.names(micro) <- substring(row.names(micro), 2)
row.names(doli) <- substring(row.names(doli), 2)
row.names(micro)

#keep date only (ie. remove everything before first period)
#change _ to .
row.names(remove_rare_vir) <- gsub("_", ".", row.names(remove_rare_vir))
#remove everything before 1st period (just to keep date)
row.names(remove_rare_vir) <- gsub("^.*?\\.","", row.names(remove_rare_vir))

remove_rare_vir <- as.data.frame(remove_rare_vir)
row.names(remove_rare_vir) <-substring(row.names(remove_rare_vir), 2)

#add micro before ASV to ID from phage
colnames(micro) <- lapply(colnames(micro), function(x) paste("micro", x, sep = "_"))
#add doli before ASV to ID from phage
colnames(doli) <- lapply(colnames(doli), function(x) paste("doli", x, sep="_"))


#merge bacterial df
merge_bacteria <- merge(micro, doli, by="row.names")
head(merge_bacteria)
merge_bact <- merge_bacteria[,-1]
rownames(merge_bact) <- merge_bacteria[,1]
head(merge_bact)

merge_bact$date <- row.names(merge_bact)
remove_rare_vir$date <- row.names(remove_rare_vir)

head(remove_rare_vir)

row.names(merge_bact)
row.names(remove_rare_vir)


#merge bacterial and viral 
bact_phage_tab <- merge(merge_bact, remove_rare_vir, by="date", na.rm=F)
View(bact_phage_tab)

write.csv(bact_phage_tab, "gLV_table.csv", sep="\t")




