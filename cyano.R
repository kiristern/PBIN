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

#tranform cyano counts to relative abundance
cyano_removed <-decostand(cyano_removed, method="hellinger")










