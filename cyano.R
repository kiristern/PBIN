setwd("~/Documents/GitHub/PBIN/data/cyano")

cyano_counts <- read.table("Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
cyano_taxa <- read.table("Champ_ASVs_taxonomy.txt", sep = "", fill = TRUE, header=T)
cyano <- read.table("Champlain_cyano.txt", fill = TRUE)
vir_data <- abundance_removed

#extract viral data sample dates
row.names(vir_data) <- sub("*._*._*._*._*._*._*._","", row.names(vir_data))
#change "_" to "."
row.names(vir_data) <- gsub("_", ".", row.names(vir_data))
#remove everything after last "."
row.names(vir_data)<-sub(".[^.]+$", "", row.names(vir_data))

#remove X at beginning of date
colnames(cyano_counts)[1:135] <- substring(colnames(cyano_counts)[1:135], 2)
#remove everything after last period in colnames
colnames(cyano_counts)[c(1, 5:63, 65:135)] <- sub(".[^.]+$", "", colnames(cyano_counts)[c(1, 5:63, 65:135)])

#compare samples in common between vir and cyano
row.names(vir_data) %in% colnames(cyano_counts)

#specific samples that are not the same
(cols_remove <- setdiff(colnames(cyano_counts), row.names(vir_data)))
#count how many are different
length(setdiff(colnames(cyano_counts), row.names(vir_data)))

#remove rows (samples) that aren't in env_var from abundance
cyano_removed <- cyano_counts[, !(colnames(cyano_counts) %in% cols_remove)]
dim(cyano_removed)

head(cyano_removed)
#write.csv(cyano_removed, "relevantsamples.csv")

cyano_names <- read.csv("relevantsamples.csv")
head(cyano_names)
head(cyano_taxa)

#merge dfs
count_taxa <- left_join(cyano_names, cyano_taxa, "ASV")

#group by cyanobacteria, microcystis, dolichospermum
class_keep <- "Cyanobacteria"
genus_keep <- "Microcystaceae"
sp_keep <- "Dolichospermum"

head(cyanobacteria <- count_taxa[count_taxa$Class %in% class_keep, ])
head(microcystis <- count_taxa[count_taxa$Genus %in% genus_keep, ])


#transpose table
cyano_transpose <- t(cyano_removed)
#get total count of asv
cyano_tot <- colSums(cyano_transpose)
#transform to df
cyano_tot <- as.data.frame(cyano_tot)
#extract asv_id 
cyano_asv_id <- row.names(cyano_tot)
#add asv_ids to df
cyano_tot$ID=cyano_asv_id
#add new empty column
newcol <- "rel_ab"
cyano_tot[,newcol] <- NA

#get relative abundance function
get_rel_abun <- function(x){
  x / sum(asv_tot[1])
}



