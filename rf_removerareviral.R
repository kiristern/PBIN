#removed taxa not seen more than 3 times in at least 5% of the samples
remove_rare_vir <- t(filt_vir)

bacterial <- read.csv("cyano/bacteria_metagenonet.csv", row.names=1, header=T)
bact_taxa <- read.csv("cyano/ASVs_taxonomy_Champ_Greengenes.csv")

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
dim(cyano)


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

#Get sums
(sum_cyano <- as.data.frame(rowSums(cyano)))
sum_micro <- as.data.frame(rowSums(micro))
sum_doli <- as.data.frame(rowSums(doli))

#remove X in front of rownames
row.names(sum_cyano) <- substring(row.names(sum_cyano), 2)
row.names(sum_micro) <- substring(row.names(sum_micro), 2)
row.names(sum_doli) <- substring(row.names(sum_doli), 2)
row.names(sum_cyano)

#add col with sample 
sum_cyano$date <- row.names(sum_cyano)
head(sum_cyano)
sum_micro$date <- row.names(sum_micro)
sum_doli$date <- row.names(sum_doli)

#keep date only (ie. remove everything before first period)
row.names(remove_rare_vir) <- sub("*._*._*._*._*._*._*._","", row.names(remove_rare_vir))
#change all _ to .
row.names(remove_rare_vir) <- gsub("_", ".", row.names(remove_rare_vir))
row.names(remove_rare_vir)
#transform data
dim(remove_rare_vir)
viral_helli <-decostand(remove_rare_vir, method="hellinger", na.rm = T)
viral_helli <- na.omit(viral_helli)
dim(viral_helli)

#assign new name to the duplicated sample/date at position 100
# rownames(viral_helli)[99]
# rownames(viral_helli)[100]
# rownames(viral_helli)[100] <- "15.05.2011.2"

#check if this exists in sum_cyano df
sum_cyano["15.05.2011.2", ] #it does not
#remove
# viral_helli_rem <- viral_helli[-100,]
# row.names(viral_helli_rem)

# dim(viral_helli_rem)
dim(sum_cyano)

#check rows from bacterial that are in viral
row.names(viral_helli) %in% row.names(sum_cyano)
row.names(viral_helli) %in% row.names(sum_micro)
row.names(viral_helli) %in% row.names(sum_doli)

#specific samples that are the same
(keep_cyano <- which(row.names(viral_helli) %in% row.names(sum_cyano)))
keep_micro <- which(row.names(viral_helli) %in% row.names(sum_micro))
keep_doli <- which(row.names(viral_helli) %in% row.names(sum_doli))
#viral samples to keep
viral4cyano_keep <- viral_helli[keep_cyano, ]
viral4micro_keep <- viral_helli[keep_micro, ]
viral4doli_keep <- viral_helli[keep_doli, ]

dim(viral4cyano_keep)
dim(sum_cyano)
row.names(viral4cyano_keep)
row.names(sum_cyano)

#rows that are in sum_bacterial but not in viral_keep
row.names(sum_cyano) %in% row.names(viral4cyano_keep)
row.names(sum_micro) %in% row.names(viral4micro_keep)
row.names(sum_doli) %in% row.names(viral4doli_keep)

#rows to remove BETTER TP KEEP ROWS RATHER THAN DROPPING TO AVOID ACCIDENTALLY RE-RUNNING CODE TWICE
(rowrem_cyano <- setdiff(row.names(sum_cyano), row.names(viral4cyano_keep)))
rowrem_micro <- setdiff(row.names(sum_micro), row.names(viral4micro_keep))
rowrem_doli <- setdiff(row.names(sum_doli), row.names(viral4doli_keep))
#remove rows (samples) that aren't in bacteria from viral
sumcyano_keep <- as.data.frame(sum_cyano[!(row.names(sum_cyano) %in% rowrem_cyano), ])
summicro_keep <- as.data.frame(sum_micro[!(row.names(sum_micro) %in% rowrem_micro), ])
sumdoli_keep <- as.data.frame(sum_doli[!(row.names(sum_doli) %in% rowrem_doli), ])
# sumcyano_keep <- sum_cyano[!(row.names(sum_cyano) %in% rowrem_cyano), ]
# summicro_keep <- sum_micro[!(row.names(sum_micro) %in% rowrem_micro), ]
# sumdoli_keep <- sum_doli[!(row.names(sum_doli) %in% rowrem_doli), ]

# #rename col names
names(sumcyano_keep)[names(sumcyano_keep)=="rowSums(cyano)"] <- "sum"
names(summicro_keep)[names(summicro_keep)=="rowSums(micro)"] <- "sum"
names(sumdoli_keep)[names(sumdoli_keep)=="rowSums(doli)"] <- "sum"



dim(viral4cyano_keep)
dim(sumcyano_keep)
row.names(sumcyano_keep) %in% row.names(viral4cyano_keep)
row.names(viral4cyano_keep) %in% row.names(sumcyano_keep)

duprownambac <- row.names(sumcyano_keep) 
duplicate(duprownambac)

#### script returns ####
sumcyano_keep
summicro_keep
sumdoli_keep

viral4cyano_keep
viral4micro_keep
viral4doli_keep



#### Random Forest ####

#divide into training/test sets
samp_size <- floor(0.70 * nrow(viral4cyano_keep))
set.seed(1234)
train_idx <- sample(seq_len(nrow(viral4cyano_keep)), size = samp_size)
# train <- viral_helli[train_idx, ]
# test <- viral_helli[-train_idx, ]

(cyano_rf <- randomForest(sumcyano_keep ~ ., 
                          data = viral4cyano_keep, 
                          subset = train_idx))

#Plotting the Error vs Number of Trees Graph
plot(cyano_rf)

#look at the importance that the classifier has assigned to each variable
varImpPlot(cyano_rf)


### micro RF  ####
(micro_rf <- randomForest(summicro_keep$sum ~ ., #remove $sum if didn't as.df
                          data = viral4micro_keep, 
                          subset = train_idx))

#Plotting the Error vs Number of Trees Graph
plot(micro_rf)

#look at the importance that the classifier has assigned to each variable
varImpPlot(micro_rf)



### doli RF  ####
(doli_rf <- randomForest(sumdoli_keep$sum ~ ., #remove $sum if didn't as.df
                         data = viral4doli_keep, 
                         subset = train_idx))

#Plotting the Error vs Number of Trees Graph
plot(doli_rf)

#look at the importance that the classifier has assigned to each variable
varImpPlot(doli_rf)
