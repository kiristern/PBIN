viral <- read.table("../ASVs_counts_copy.tsv", row.names = 1, header=T)
viral <- t(viral)
bacterial <- read.csv("bacteria_metagenonet.csv", row.names=1, header=T)
bact_taxa <- read.csv("ASVs_taxonomy_Champ_Greengenes.csv")

#### transform bacterial into relative abundance ####
#function to get relative abundance
get_rel_abund <- function(x){
  x / sum(x)
}
#apply function to bacterial df 
bacterial_relab <- get_rel_abund(bacterial)

#### merge taxa ID to ASVs ####
#assign ASV_X as rownames in bact_taxa df
bact_taxa_rownames <- bact_taxa[,-1] 
rownames(bact_taxa_rownames) <- bact_taxa[,1]
#merge df by rownames
bact <- merge(bacterial_relab, bact_taxa_rownames, by="row.names")

#### group by cyanobacteria, microcystis, dolichospermum ####
cyanobacteria <- bact[grep("Cyanobacteria", bact$Phylum), ]
cyanobacteria <- cyanobacteria[1:136]
cyano <- cyanobacteria[,-1]
rownames(cyano) <- cyanobacteria[,1]
cyano <- t(cyano)

microcystis <- bact[grep("Microcystaceae", bact$Family), ]
microcystis <- microcystis[1:136]
micro <- microcystis[,-1]
rownames(micro) <- microcystis[,1]
micro <- t(micro)

dolichospermum <- bact[grep("Dolichospermum", bact$Genus), ]
dolichospermum <- dolichospermum[1:136]
doli <- dolichospermum[,-1]
rownames(doli) <- dolichospermum[,1]
doli <- t(doli)

#Get sums
sum_cyano <- as.data.frame(rowSums(cyano))
sum_micro <- as.data.frame(rowSums(micro))
sum_doli <- as.data.frame(rowSums(doli))

#remove X in front of rownames
row.names(sum_cyano) <- substring(row.names(sum_cyano), 2)
row.names(sum_micro) <- substring(row.names(sum_micro), 2)
row.names(sum_doli) <- substring(row.names(sum_doli), 2)

#keep date only (ie. remove everything before first period)
row.names(viral) <- sub("*._*._*._*._*._*._*._","", row.names(viral))
#change all _ to .
row.names(viral) <- gsub("_", ".", row.names(viral))

#check rows from bacterial that are in viral
row.names(viral) %in% row.names(sum_cyano)
row.names(viral) %in% row.names(sum_micro)
row.names(viral) %in% row.names(sum_doli)

#specific samples that are not the same
vir_rowrem_cyano <- setdiff(row.names(viral), row.names(sum_cyano))
vir_rowrem_micro <- setdiff(row.names(viral), row.names(sum_micro))
vir_rowrem_doli <- setdiff(row.names(viral), row.names(sum_doli))

#remove rows (samples) that aren't in bacteria from viral
viral_cyano <- viral[!(row.names(viral) %in% vir_rowrem_cyano), ]
viral_micro <- viral[!(row.names(viral) %in% vir_rowrem_micro), ]
viral_doli <- viral[!(row.names(viral) %in% vir_rowrem_doli), ]

row.names(viral_cyano)
row.names(viral)



#### Dolichospermum ####

#merge dolico with viral
doli <- merge(bacterial[1], viral, by="row.names", all=T)
#transform col1 into row.names
doli2 <- doli[,-1]
rownames(doli2)<-doli[,1]
doli_rem_samp <- doli[,-1]
#transform data
dim(doli)
doli_trans <-decostand(doli_rem_samp, method="hellinger", na.rm = T)
doli_trans_nona <- na.omit(doli_trans)
dim(doli_trans_nona)

#divide into training/test sets
samp_size <- floor(0.70 * nrow(doli_trans_nona))
set.seed(1234)
train_idx <- sample(seq_len(nrow(doli_trans_nona)), size = samp_size)
# train_doli <- doli_trans_nona[train_idx, ]
# test_doli <- doli_trans_nona[-train_idx, ]

(doli_rf <- randomForest(Dolicho.Abundance ~ ., data = doli_trans_nona, subset = train_idx))

#Plotting the Error vs Number of Trees Graph
plot(doli_rf)

#look at the importance that the classifier has assigned to each variable
varImpPlot(doli_rf)


#### Microcystis ####

#merge micro with viral
mic <- merge(bacterial[2], viral, by="row.names", all=T)
#transform col1 into row.names
mic2 <- mic[,-1]
rownames(mic2)<-mic[,1]
#transform data
mic_trans <-decostand(mic2, method="hellinger", na.rm = T)

