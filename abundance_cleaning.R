library(dplyr)
library(tibble)
library(tidyr)
library(vegan)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table
ASV_count <- t(vir_abun_removed)

#transpose table
asv <- t(ASV_count)
#get total count of asv
asv_tot <- colSums(asv)
#transform to df
asv_tot <- as.data.frame(asv_tot)
#extract asv_id 
asv_id <- row.names(asv_tot)
#add asv_ids to df
asv_tot$ID=asv_id
#add new empty column
newcol <- "rel_ab"
asv_tot[,newcol] <- NA
# asv_tot <- mutate(asv_tot, sample=row.names(asv_tot))

#get relative abundance function
get_rel_abun <- function(x){
  x / sum(x)
}

#apply function to the first col of df asv_tot and put into rel_ab col of df asv_tot
asv_tot[3] <- get_rel_abun(asv_tot[1])

#sort relative abundance from largest to smallest
asv_tot <- arrange(asv_tot, desc(rel_ab))
top20 <- head(asv_tot, 20)
top20
#write.csv(top20, file="top20_filt_removed.csv")

top20 <- read.csv("top20_filt_removed.csv")

#select asvs 1-21 excluding asv17
asv_tax <- as.data.frame(asv) %>% select(
  "ASV_2",
  "ASV_1",
  "ASV_3",
  "ASV_5",
  "ASV_4",
  "ASV_7",
  "ASV_8",
  "ASV_12",
  "ASV_6",
  "ASV_21",
  "ASV_15",
  "ASV_33",
  "ASV_11",
  "ASV_18",
  "ASV_17",
  "ASV_19",
  "ASV_48",
  "ASV_13",
  "ASV_25"
)
nrow(asv_tax)
ncol(asv_tax)
#relative abundance matrix
asv_rel_abun <- decostand(asv_tax, method="total")
asv_rel_abun <- as.data.frame(asv_rel_abun)

#change col names to taxa ID
#name_tax <- select(top17, ID, virus)
  #if select doesn't work
  name_tax <- top20[, c(2, 4)]

#add brackets around ID
name_tax$ID_brackets <- with(name_tax, paste0("(", ID, ")"))
  
#merge ASV name with ID in new col
(virus_ID <- paste(name_tax$taxa, name_tax$ID_brackets, sep=" "))
virus_ID[1] <- "Unknown ASV (ASV_2)"
virus_ID[3] <- "Unknown ASV (ASV_3)"


#add col to name_tax df
name_tax <- add_column(name_tax, virus_ID)

#change ASV_ to real taxonomic name
names(asv_rel_abun) <- name_tax$virus_ID[match(names(asv_rel_abun), name_tax$ID)]
head(asv_rel_abun)
####Combine ASV_14, 16, 2, 3 into unknown ASV

# asv_rel_abun$Unknown_ASV <- asv_rel_abun %>% select(ASV_2, ASV_3, ASV_14, ASV_16) %>% rowSums()
# 
# #remove asv_X cols
# asv_rel_abun_unk <- asv_rel_abun %>% select(-c(ASV_2, ASV_3, ASV_14, ASV_16))
# ncol(asv_rel_abun_unk)
# #duplicate each sample 17 times
# asv_rel_abun_unk_dup <- asv_rel_abun_unk[rep(seq_len(nrow(asv_rel_abun_unk)), each = 17), ]
# 
# #create new row of repeated ASV_IDs 168 times (# of samples) - make sure order is from ASV_1 - end (not order of abundance)
# ASV_ID_unk <- data.frame(ASV_ID = c("Uncultured cyanophage clone KRB1008M5 (ASV_1)", 
#                                 "Uncultured cyanomyovirus clone SZCPS18 (ASV_4)",
#                                 "Uncultured Myoviridae g20_82_56_1%_NEQ (ASV_5)", 
#                                 "Uncultured cyanophage clone BwC24 (ASV_6)", 
#                                 "Uncultured Myoviridae g20_94_44_1%_NEQ (ASV_7)", 
#                                 "Uncultured cyanomyovirus clone 88268CPSCC8 (ASV_8)",
#                                 "Uncultured cyanophage clone LAB_g20_b26_C12 (ASV_9)",
#                                 "Uncultured cyanophage clone KRB1008M5 (ASV_10)",
#                                 "Uncultured cyanophage clone LAB_g20_b28_D9 (ASV_11)",
#                                 "Uncultured cyanomyovirus clone 45202CPSCC1 (ASV_12)",
#                                 "Uncultured cyanophage clone KRA0808M1 (ASV_13)",
#                                 "Cyanophage LLM-cg20-1 clone LLM-cg20-1 (ASV_15)",
#                                 "Cyanophage S-RIM isolate S-RIM5 (ASV_18)",
#                                 "Uncultured Myoviridae clone KUSW41 (ASV_19)",
#                                 "Synechococcus phage S-LKM3 (ASV_20)",
#                                 "Uncultured Myoviridae g20_61_47_1%_UW (ASV_21)",
#                                 "Unknown_ASV"
# ))
# n=168
# replicate_unk <- do.call("rbind", replicate(n, ASV_ID_unk, simplify = FALSE))
# 
# write.csv(replicate_unk, "replicate_unk.csv")
# replicate_unk <- read.csv("replicate_unk.csv")
# 
# #create new col named samples which duplicates rownames
# samples_unk <- row.names(asv_rel_abun_unk_dup)
# #add col to df
# asv_rel_abun_unk_dup <- add_column(asv_rel_abun_unk_dup, samples_unk, .before=1)
# 
# #remove everything after "."
# asv_rel_abun_unk_dup$samples_unk <- gsub("\\..*", "", asv_rel_abun_unk_dup$samples_unk)
# #write.csv(asv_rel_abun_dup, file = "repeat_sample_names.csv")
# 
# #transpose
# sample_abundance_unk <- as.data.frame(t(asv_rel_abun_unk))
# 
# #one col of all abundances
# stacked_unk <- stack(sample_abundance_unk)
# 
# #add stacked to asv_rel_abun_dup df
# asv_rel_abun_unk_dup <- add_column(asv_rel_abun_unk_dup, stacked_unk$values, .after=1)
# write.csv(asv_rel_abun_unk_dup, file="asv_rel_abun_unk.csv")
# 






#duplicate each sample 20 times (number of unique ASVs)
asv_rel_abun_dup <- asv_rel_abun[rep(seq_len(nrow(asv_rel_abun)), each = 20), ]

#create new row of repeated ASV_IDs 168 times (# of samples) - make sure order is from ASV_1 - end (not order of abundance)
# ASV_ID <- data.frame(ASV_ID = c("Uncultured cyanophage clone KRB1008M5 (ASV_1)", 
#                                 "Unknown ASV (ASV_2)", 
#                                 "Unknown ASV (ASV_3)",
#                                 "Uncultured cyanomyovirus clone SZCPS18 (ASV_4)",
#                                 "Uncultured Myoviridae g20_82_56_1%_NEQ (ASV_5)", 
#                                 "Uncultured cyanophage clone BwC24 (ASV_6)", 
#                                 "Uncultured Myoviridae g20_94_44_1%_NEQ (ASV_7)", 
#                                 "Uncultured cyanomyovirus clone 88268CPSCC8 (ASV_8)",
#                                 "Uncultured cyanophage clone LAB_g20_b26_C12 (ASV_9)",
#                                 "Uncultured cyanophage clone KRB1008M5 (ASV_10)",
#                                 "Uncultured cyanophage clone LAB_g20_b28_D9 (ASV_11)",
#                                 "Uncultured cyanomyovirus clone 45202CPSCC1 (ASV_12)",
#                                 "Uncultured cyanophage clone KRA0808M1 (ASV_13)",
#                                 "Unknown ASV (ASV_14)",
#                                 "Cyanophage LLM-cg20-1 clone LLM-cg20-1 (ASV_15)",
#                                 "Unknown ASV (ASV_16)",
#                                 "Cyanophage S-RIM isolate S-RIM5 (ASV_18)",
#                                 "Uncultured Myoviridae clone KUSW41 (ASV_19)",
#                                 "Synechococcus phage S-LKM3 (ASV_20)",
#                                 "Uncultured Myoviridae g20_61_47_1%_UW (ASV_21)"
# ))

ASV_ID <- data.frame(ASV_ID = c("Uncultured cyanophage clone KRB1008M5 (ASV_1)", 
                                "ASV_2 (ASV_2)", 
                                "ASV_3 (ASV_3)",
                                "Uncultured cyanomyovirus clone SZCPS18 (ASV_4)",
                                "Uncultured Myoviridae g20_82_56_1%_NEQ (ASV_5)", 
                                "Uncultured cyanophage clone BwC24 (ASV_6)", 
                                "Uncultured Myoviridae g20_94_44_1%_NEQ (ASV_7)", 
                                "Uncultured cyanomyovirus clone 88268CPSCC8 (ASV_8)",
                                "Uncultured cyanophage clone LAB_g20_b26_C12 (ASV_9)",
                                "Uncultured cyanophage clone KRB1008M5 (ASV_10)",
                                "Uncultured cyanophage clone LAB_g20_b28_D9 (ASV_11)",
                                "Uncultured cyanomyovirus clone 45202CPSCC1 (ASV_12)",
                                "Uncultured cyanophage clone KRA0808M1 (ASV_13)",
                                "ASV_14 (ASV_14)",
                                "Cyanophage LLM-cg20-1 clone LLM-cg20-1 (ASV_15)",
                                "ASV_16 (ASV_16)",
                                "Cyanophage S-RIM isolate S-RIM5 (ASV_18)",
                                "Uncultured Myoviridae clone KUSW41 (ASV_19)",
                                "Synechococcus phage S-LKM3 (ASV_20)",
                                "Uncultured Myoviridae g20_61_47_1%_UW (ASV_21)"
))
n=168
replicate <- do.call("rbind", replicate(n, ASV_ID, simplify = FALSE))

#write.csv(replicate, "replicate_test.csv")
#replicate <- read.csv("replicate_test.csv")

#add ASV_ID col to asv_rel_abun df
#asv_rel_abun_dup <- add_column(asv_rel_abun_dup, replicate, .before=1)

#create new col named samples which duplicates rownames
samples <- row.names(asv_rel_abun_dup)
#add col to df
asv_rel_abun_dup <- add_column(asv_rel_abun_dup, samples, .before=1)

#remove everything after "."
asv_rel_abun_dup$samples <- gsub("\\..*", "", asv_rel_abun_dup$samples)
# write.csv(asv_rel_abun_dup, file = "repeat_sample_names_test.csv")

#transpose
sample_abundance <- as.data.frame(t(asv_rel_abun))

#one col of all abundances
stacked <- stack(sample_abundance)
head(stacked)
#add stacked to asv_rel_abun_dup df
asv_rel_abun_dup <- add_column(asv_rel_abun_dup, stacked$values, .after=1)
head(asv_rel_abun_dup)

#select certain cols only
df <- asv_rel_abun_dup %>% select("samples", "stacked$values")
#add replicate names to df
df <- add_column(df, replicate, .before=1)
head(df)
#rename cols
df <- df %>%
  rename(
    abundance = "stacked$values"
  )
# write.csv(asv_rel_abun_dup, file="asv_rel_abun_rewrite.csv")

### TODO ADD REPLICATE COLUMN TO ASV_REL_ABUN FILE AND DELETE UNNECESSARY COLS

#final_asv_rel_abun <- select(asv_rel_abun_dup, "samples", "ASV_ID", "stacked$values")
# names(final_asv_rel_abun)[names(final_asv_rel_abun)=="stacked$values"] <- "abundance"
# 
# write.csv(final_asv_rel_abun, file="relative_abundance.csv")








#redo with filtered ASV (in at least 5% of samples)
filt_ASV_count <- filt_vir_seq %>% otu_table()

#transpose table
asv_filt <- t(filt_ASV_count)
#get total count of asv
asv_tot_filt <- colSums(asv_filt)
#transform to df
asv_tot_filt <- as.data.frame(asv_tot_filt)
#extract asv_id 
asv_id_filt <- row.names(asv_tot_filt)
#add asv_ids to df
asv_tot_filt$ID=asv_id_filt
#add new empty column
newcol <- "rel_ab"
asv_tot_filt[,newcol] <- NA
# asv_tot <- mutate(asv_tot, sample=row.names(asv_tot))

#apply function to the first col of df asv_tot and put into rel_ab col of df asv_tot
asv_tot_filt[3] <- get_rel_abun(asv_tot_filt[1])

#sort relative abundance from largest to smallest
asv_tot_filt <- arrange(asv_tot_filt, desc(rel_ab))
top20_filt <- head(asv_tot_filt, 20) #same as before, as predicted!! good. 


