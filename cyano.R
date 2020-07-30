setwd("~/Documents/GitHub/PBIN/data/cyano")

cyano_counts <- read.table("Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
cyano_taxa <- read.csv("ASVs_taxonomy_Champ_Greengenes.csv", header = T, row.names = 1, fill=T)
cyano <- read.table("Champlain_cyano.txt", fill = TRUE)
#viral taxa not seen more than 3 times in at least 10% of the samples, transformed hellinger, and removed samples that aren't in env_var(meta)  
vir_data <- vir_abun_removed

#extract viral data sample dates
row.names(vir_data)
row.names(vir_data) <- sub("*._*._*._*._*._*._*._","", row.names(vir_data))
#change "_" to "."
row.names(vir_data) <- gsub("_", ".", row.names(vir_data))
row.names(vir_data)
#remove everything after last "."
# row.names(vir_data)<-sub(".[^.]+$", "", row.names(vir_data))

colnames(cyano_counts)
#remove X at beginning of date
colnames(cyano_counts)[1:135] <- substring(colnames(cyano_counts)[1:135], 2)
cyano_counts <- cyano_counts %>% select(1:135)
#remove everything after last period in colnames
# colnames(cyano_counts)[c(1, 5:63, 65:135)] <- sub(".[^.]+$", "", colnames(cyano_counts)[c(1, 5:63, 65:135)])

#compare samples in common between vir and cyano
colnames(cyano_counts) %in% row.names(vir_data)

#specific samples that are not the same
(cols_remove <- setdiff(colnames(cyano_counts), row.names(vir_data)))
#count how many are different
length(setdiff(colnames(cyano_counts), row.names(vir_data)))

#remove rows (samples) that aren't in env_var from abundance
cyano_removed <- cyano_counts[, !(colnames(cyano_counts) %in% cols_remove)]
dim(cyano_removed)
head(cyano_removed)
#write.csv(cyano_removed, "relevantsamples.csv")

cyano_names <- cyano_removed #read.csv("relevantsamples.csv", row.names = 1, header = T)
head(cyano_names)
head(cyano_taxa)

#create col with ASV in cyano_taxa in order to merge with second df
cyano_taxa$ASV <- row.names(cyano_taxa)

#transform asv density as a proportion of the sum of all densities
cyano_rel_abundance <-decostand(cyano_removed, method="hellinger")
# #transpose table
# cyano_transpose <- t(cyano_rel_abundance)
# 
# #get total count of asv
# cyano_tot <- colSums(cyano_transpose)
# #transform to df
# cyano_tot <- as.data.frame(cyano_tot)
# #extract asv_id
# cyano_asv_id <- row.names(cyano_tot)
# #add asv_ids to df
# cyano_tot$ASV=cyano_asv_id
# #add new empty column
# newcol <- "tot_rel_ab"
# cyano_tot[,newcol] <- NA
# 
# #get relative abundance function
# get_rel_abund <- function(x){
#   x / sum(x)
# }
# 
# #apply function to the first col of df asv_tot and put into rel_ab col of df asv_tot
# cyano_tot[3] <- get_rel_abun(cyano_tot[1])
# 
# #merge dfs
# count_taxa <- left_join(cyano_tot, cyano_taxa, "ASV")
# head(count_taxa)
# 
# #tranpose cyano_rel_abundance back
# cyano_trans <- t(cyano_rel_abundance)
# cyano_trans <- as.data.frame(cyano_trans)
# #create new col with replicated ASV name
# cyano_trans$ASV=cyano_asv_id

head(cyano_rel_abundance)
#extract asv_id
cyano_asv_id <- row.names(cyano_rel_abundance)
#add asv_ids to df
cyano_rel_abundance$ASV=cyano_asv_id

#merge dfs
cyano_df <- left_join(cyano_rel_abundance, count_taxa, "ASV")

#set rownames back to ASV
cyano_df %>% remove_rownames() %>% column_to_rownames(var="ASV")

#group by cyanobacteria, microcystis, dolichospermum
head(cyanobacteria <- cyano_df[grep("Cyanobacteria", cyano_df$Phylum), ])

colnames(cyanobacteria)
#keep only certain cols
cyano_samples <- cyanobacteria %>% select(1:50)
#set rownames back to ASV
cyano_samples %>% remove_rownames() %>% column_to_rownames(var="ASV")


#########  CYANO SAMPLES TABLE TRANSFORMED BY HELLINGER AND REMOVED SAMPLES <10% ##########
cyano_samples



#remove everything after last period in colnames
#check which ones end in 1; "$" anchors 1 to the end -- ATTN: incl dates that end in 2011, need to sort
grep(".1$", colnames(cyano_samples))
# colnames(cyano_samples)[c(2,4,7,9,11,13,
                          # 15,18,20,22,24,
                          # 26,28,30,33,35,
                          # 37,40,42,45,47,
                          # 49,51,53,55,57,
                          # 59,61,64,67,69,71)] <- sub(".1$", "", colnames(cyano_samples)[c(2,4,7,9,11,13,
                          #                                                                 15,18,20,22,24,
                          #                                                                 26,28,30,33,35,
                          #                                                                 37,40,42,45,47,
                          #                                                                 49,51,53,55,57,
                          #                                                                 59,61,64,67,69,71)])
# #remove X at beginning of date
# colnames(cyano_samples) <- sub(".","",colnames(cyano_samples))

#format as date
colnames(cyano_samples) <- as.Date(colnames(cyano_samples), "%d.%m.%Y")
#write.csv(cyano_samples, "cyano_samples.csv")

microcystis <- cyano_df[grep("Microcystaceae", cyano_df$Family), ]
colnames(microcystis)
micro_samples <- microcystis %>% select(1:66)
grep(".1$", colnames(micro_samples))
# colnames(micro_samples)[c(2,4,7,9,11,13,
                          # 15,18,20,22,24,
                          # 26,28,30,33,35,
                          # 37,40,42,45,47,
                          # 49,51,53,55,57,
                          # 59,61,64,67,69,71)] <- sub(".1$", "", colnames(micro_samples)[c(2,4,7,9,11,13,
                          #                                                                 15,18,20,22,24,
                          #                                                                 26,28,30,33,35,
                          #                                                                 37,40,42,45,47,
                          #                                                                 49,51,53,55,57,
                          #                                                                 59,61,64,67,69,71)])
colnames(micro_samples) <- as.Date(colnames(micro_samples), "%d.%m.%Y")
# write.csv(micro_samples, "micro_samples.csv")



dolichospermum <- cyano_df[grep("Dolichospermum", cyano_df$Genus), ]
doli_samples <- dolichospermum %>% select(1:66)
# colnames(doli_samples)[c(2,4,7,9,11,13,
                          # 15,18,20,22,24,
                          # 26,28,30,33,35,
                          # 37,40,42,45,47,
                          # 49,51,53,55,57,
                          # 59,61,64,67,69,71)] <- sub(".1$", "", colnames(doli_samples)[c(2,4,7,9,11,13,
                          #                                                                 15,18,20,22,24,
                          #                                                                 26,28,30,33,35,
                          #                                                                 37,40,42,45,47,
                          #                                                                 49,51,53,55,57,
                          #                                                                 59,61,64,67,69,71)])
colnames(doli_samples) <- as.Date(colnames(doli_samples), "%d.%m.%Y")
# write.csv(doli_samples, "doli_samples.csv")

