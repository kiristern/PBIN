library(dplyr)
library(tibble)

setwd("~/Documents/GitHub/PBIN/data")

#upload ASV count table and metadata
ASV_count <- read.table("ASVs_counts_copy.tsv", row.names = 1, header=T)
meta <- read.csv("metadata3.csv", row.names=1, header=T)

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
  x / sum(asv_tot[1])
}

#apply function to the first col of df asv_tot and put into rel_ab col of df asv_tot
asv_tot[3] <- get_rel_abun(asv_tot[1])

#sort relative abundance from largest to smallest
asv_tot <- arrange(asv_tot, desc(rel_ab))
top17 <- head(asv_tot, 17)

#write.csv(top17, file="top17_new.csv")

top17 <- read.csv("top17.csv")

asv_tax <- (asv)[,1:18]
asv_tax <- asv_tax[,-17]
nrow(asv_tax)
ncol(asv_tax)
#relative abundance matrix
asv_rel_abun <- decostand(asv_tax, method="total")
asv_rel_abun <- as.data.frame(asv_rel_abun)

#change col names to taxa ID
name_tax <- select(top17, ID, virus)

#change ASV_ to real taxonomic name
names(asv_rel_abun) <- name_tax$virus[match(names(asv_rel_abun), name_tax$ID)]

#duplicate each sample 17 times (number of unique ASVs)
asv_rel_abun_dup <- asv_rel_abun[rep(seq_len(nrow(asv_rel_abun)), each = 17), ]

#create new row of repeated ASV_IDs 168 times (# of samples) - make sure order is from ASV_1 - end (not order of abundance)
ASV_ID <- data.frame(ASV_ID = c("Uncultured cyanophage clone KRB1008M5", 
                                "ASV_2", 
                                "ASV_3",
                                "Uncultured cyanomyovirus clone SZCPS18",
                                "Uncultured Myoviridae g20_82_56_1%_NEQ", 
                                "Uncultured cyanophage clone BwC24", 
                                "Uncultured Myoviridae g20_94_44_1%_NEQ", 
                                "Uncultured cyanomyovirus clone 88268CPSCC8",
                                "Uncultured cyanophage clone LAB_g20_b26_C12",
                                "Uncultured cyanophage clone KRB1008M5",
                                "Uncultured cyanophage clone LAB_g20_b28_D9",
                                "Uncultured cyanomyovirus clone 45202CPSCC1",
                                "Uncultured cyanophage clone KRA0808M1",
                                "ASV_14",
                                "Cyanophage LLM-cg20-1 clone LLM-cg20-1",
                                "ASV_16",
                                "Cyanophage S-RIM isolate S-RIM5"
))
n=168
replicate <- do.call("rbind", replicate(n, ASV_ID, simplify = FALSE))

#write.csv(replicate, "replicate.csv")
replicate <- read.csv("replicate.csv")

#add ASV_ID col to asv_rel_abun df
#asv_rel_abun_dup <- add_column(asv_rel_abun_dup, replicate, .before=1)

#create new col named samples which duplicates rownames
samples <- row.names(asv_rel_abun_dup)
#add col to df
asv_rel_abun_dup <- add_column(asv_rel_abun_dup, samples, .before=1)

#remove everything after "."
asv_rel_abun_dup$samples <- gsub("\\..*", "", asv_rel_abun_dup$samples)
#write.csv(asv_rel_abun_dup, file = "repeat_sample_names.csv")

#transpose
sample_abundance <- as.data.frame(t(asv_rel_abun))

#one col of all abundances
stacked <- stack(sample_abundance)

#add stacked to asv_rel_abun_dup df
asv_rel_abun_dup <- add_column(asv_rel_abun_dup, stacked$values, .after=1)
#write.csv(asv_rel_abun_dup, file="asv_rel_abun.csv")

# final_asv_rel_abun <- select(asv_rel_abun_dup, "samples", "ASV_ID", "stacked$values")
# names(final_asv_rel_abun)[names(final_asv_rel_abun)=="stacked$values"] <- "abundance"
# 
# write.csv(final_asv_rel_abun, file="relative_abundance.csv")