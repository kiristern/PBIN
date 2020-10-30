library(dplyr)
library(tibble)
library(tidyr)
library(vegan)

setwd("~/Documents/GitHub/PBIN/data")


#filter by site
meta_pelagic <- meta %>%
  filter(Site == "Pelagic")
nrow(meta_pelagic)

meta_littoral <- meta %>%
  filter(Site == "Littoral")
nrow(meta_littoral)

#filter viral count by littoral / pelagic (ie. select cols that match rownames from litt. or pel.)
vir_count_littoral <- as.data.frame(ASV_count[, (colnames(ASV_count) %in% rownames(meta_littoral))])
#write.table(vir_count_littoral, "ASV_count_littoral.txt", row.names = T, quote = F, sep = "\t")

vir_count_pelagic <- ASV_count[, (colnames(ASV_count) %in% rownames(meta_pelagic))]
#write.table(vir_count_pelagic, "ASV_count_pelagic.txt", row.names = T, quote = F, sep = "\t")




#transpose ASV count table
vir_lit <- t(vir_count_littoral)
vir_pel <- t(vir_count_pelagic)

#get total count of asv
asv_tot.l <- colSums(vir_lit)
asv_tot.p <- colSums(vir_pel)

#transform to df
asv_tot.l <- as.data.frame(asv_tot.l)
asv_tot.p <- as.data.frame(asv_tot.p)

#extract asv_id 
asv_id.l <- row.names(asv_tot.l)
asv_id.p <- row.names(asv_tot.p)

#add asv_ids to df
asv_tot.l$ID=asv_id.l
asv_tot.p$ID=asv_id.p

#add new empty column
newcol <- "rel_ab"
asv_tot.l[,newcol] <- NA
asv_tot.p[,newcol] <- NA

# asv_tot <- mutate(asv_tot, sample=row.names(asv_tot))

#get relative abundance function
get_rel_abun <- function(x){
  x / sum(x)
}

#apply function to the first col of df asv_tot and put into rel_ab col of df asv_tot
asv_tot.l[3] <- get_rel_abun(asv_tot.l[1])
asv_tot.p[3] <- get_rel_abun(asv_tot.p[1])

sum(asv_tot.l$rel_ab)
sum(asv_tot.p$rel_ab)


#sort relative abundance from largest to smallest
asv_tot.l <- arrange(asv_tot.l, desc(rel_ab))
asv_tot.p <- arrange(asv_tot.p, desc(rel_ab))

top20.l <- head(asv_tot.l, 20)
top20.l

top20.p <- head(asv_tot.p, 20)
top20.p

sum(top20.l$rel_ab)
sum(top20.p$rel_ab)

# write.csv(top20.l, "top20L_oct29.csv")
# write.csv(top20.p, "top20P_oct29.csv")

sort(top20.l$ID)

top20L <- read.csv("top20L_oct29.csv")
top20P <- read.csv("top20P_oct29.csv")

#select top20 ASVs from full df
asv_tax.l <- vir_lit[,(colnames(vir_lit) %in% top20.l$ID)]
nrow(asv_tax.l)
ncol(asv_tax.l)

asv_tax.p <- vir_pel[,(colnames(vir_pel) %in% top20.p$ID)]
nrow(asv_tax.p)
ncol(asv_tax.p)

#relative abundance matrix
asv_rel_abun.l <- decostand(asv_tax.l, method="total")
asv_rel_abun.l <- as.data.frame(asv_rel_abun.l)

asv_rel_abun.p <- decostand(asv_tax.p, method="total")
asv_rel_abun.p <- as.data.frame(asv_rel_abun.p)

#change col names to taxa ID
head(name_tax.l <- select(top20.l, ID, rel_ab))
head(name_tax.p <- select(top20.p, ID, rel_ab))

  #if select doesn't work
  #name_tax <- top20[, c(2, 3)]

#add brackets around ID
name_tax.l$ID_brackets <- with(name_tax.l, paste0("(", ID, ")"))
name_tax.p$ID_brackets <- with(name_tax.p, paste0("(", ID, ")"))

#merge ASV name with ID in new col
head(virus_ID.l <- paste(top20.l$Description, name_tax.l$ID_brackets, sep=" "))
head(virus_ID.p <- paste(top20.p$Description, name_tax.p$ID_brackets, sep=" "))

#add col to name_tax df
name_tax.l <- add_column(name_tax.l, virus_ID)
name_tax.p <- add_column(name_tax.p, virus_ID)


tvir <- t(virps)
nrow(tvir)

#change ASV_ to real taxonomic name
names(asv_rel_abun.l) <- name_tax.l$virus_ID[match(names(asv_rel_abun.l), name_tax.l$ID)]
colnames(asv_rel_abun.l)

names(asv_rel_abun.p) <- name_tax.p$virus_ID[match(names(asv_rel_abun.p), name_tax.p$ID)]
colnames(asv_rel_abun.p)

#duplicate each sample 20 times (number of unique ASVs)
asv_rel_abun_dup.l <- asv_rel_abun.l[rep(seq_len(nrow(asv_rel_abun.l)), each = 20), ]
asv_rel_abun_dup.p <- asv_rel_abun.p[rep(seq_len(nrow(asv_rel_abun.p)), each = 20), ]

nrow(asv_rel_abun.l)
#create new row of repeated ASV_IDs 50 times (# of samples)
ASV_ID.l <- data.frame(ASV_ID.l = c(colnames(asv_rel_abun_dup.l)))
n=94
replicate.l <- do.call("rbind", replicate(n, ASV_ID.l, simplify = FALSE))


nrow(asv_rel_abun.p)
ASV_ID.p <- data.frame(ASV_ID.p = c(colnames(asv_rel_abun_dup.p)))
n=72
replicate.p <- do.call("rbind", replicate(n, ASV_ID.p, simplify = FALSE))
#write.csv(replicate, "replicate_filt_removed.csv")
#replicate <- read.csv("replicate_test.csv")

#add ASV_ID col to asv_rel_abun df
asv_rel_abun_dup.l <- add_column(asv_rel_abun_dup.l, replicate.l, .before=1)
asv_rel_abun_dup.p <- add_column(asv_rel_abun_dup.p, replicate.p, .before=1)

#create new col named samples which duplicates rownames
samples.l <- row.names(asv_rel_abun_dup.l)
samples.p <- row.names(asv_rel_abun_dup.p)

#add col to df
asv_rel_abun_dup.l <- add_column(asv_rel_abun_dup.l, samples.l, .before=1)
asv_rel_abun_dup.p <- add_column(asv_rel_abun_dup.p, samples.p, .before=1)

#remove everything after "."
asv_rel_abun_dup.l$samples.l <- gsub("\\..*", "", asv_rel_abun_dup.l$samples.l)
asv_rel_abun_dup.p$samples.p <- gsub("\\..*", "", asv_rel_abun_dup.p$samples.p)
# write.csv(asv_rel_abun_dup, file = "repeat_sample_names_test.csv")

#transpose
sample_abundance.l <- as.data.frame(t(asv_rel_abun.l))
sample_abundance.p <- as.data.frame(t(asv_rel_abun.p))

#one col of all abundances
stacked.l <- stack(sample_abundance.l)
stacked.p <- stack(sample_abundance.p)

#add stacked to asv_rel_abun_dup df
asv_rel_abun_dup.l <- add_column(asv_rel_abun_dup.l, stacked.l$values, .after=1)
asv_rel_abun_dup.p <- add_column(asv_rel_abun_dup.p, stacked.p$values, .after=1)

#select certain cols only
df.l <- asv_rel_abun_dup.l %>% select("samples.l", "stacked.l$values")
df.p <- asv_rel_abun_dup.p %>% select("samples.p", "stacked.p$values")

#add replicate names to df
df.l <- add_column(df.l, replicate.l, .before=1)
df.p <- add_column(df.p, replicate.p, .before=1)

#rename cols
names(df.l)[names(df.l) == "stacked.l$values"] <- "abundance"
names(df.p)[names(df.p) == "stacked.p$values"] <- "abundance"

# write.csv(df, file="asvID_filt_rem_rel_abun.csv")

#make a taxa col
df.l <- df.l %>% mutate(taxa = ASV_ID.l)
df.p <- df.p %>% mutate(taxa = ASV_ID.p)

#remove brackets
df.l$taxa <- gsub("\\s*\\([^\\)]+\\)","",as.character(df.l$taxa))
df.p$taxa <- gsub("\\s*\\([^\\)]+\\)","",as.character(df.p$taxa))


head(df.l)
head(df.p)









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


