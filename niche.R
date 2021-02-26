####### clusters ASVs to form new OTUs #######
#https://github.com/benjjneb/dada2/issues/947
#https://github.com/adriaaulaICM/bbmo_niche_sea/blob/master/src/analysis/nucdist_otuclustering.R
library(tibble)
library(DECIPHER)
library(Biostrings)

# BiocManager::install("DECIPHER")
# library("DECIPHER")
seqs <- readDNAStringSet("ASVs.fa")
seqs.or <- OrientNucleotides(seqs)
aligned <- AlignSeqs(seqs.or)
d.align <- DECIPHER::DistanceMatrix(aligned)
d.align[1:5]

# gendist <- tidyr::gather(as.data.frame(d.align))
# as.data.frame(d.align)[1:5, 1:5]
# head(gendist)
# #write.csv(gendist, "gendist.csv")
# gendist$to <- paste0("ASV_", rownames(gendist))
# names(gendist)[names(gendist) == "key"] <- "from"
# dim(gendist)
dim(vir.corr.tab)
# dim(gendist.table)

ind <- which(upper.tri(d.align, diag = F), arr.ind = TRUE)
nn <- dimnames(d.align)
gendist.tab <- data.frame(from = nn[[1]][ind[, 1]],
           to = nn[[2]][ind[, 2]],
           gen.dist = d.align[ind])
head(gendist.tab, n=20)
dim(gendist.tab)
df3 <- inner_join(gendist.tab, vir.corr.tab, by = c("to", "from"))
df4 <- inner_join(gendist.tab, vir.corr.tab, by = c("to" = "from", "from" = "to"))
gendist.covar <- rbind(df3,df4) %>% unique()
dim(gendist.covar)
#write.csv(gendist.covar, "gendist.covar.csv")





cluster99 <- DECIPHER::IdClusters(
  d.align,
  method = "complete",
  cutoff = 0.01)
colnames(cluster99) <- "cluster99"
cluster98<- DECIPHER::IdClusters(
  d.align,
  method = "complete",
  cutoff = 0.02)
colnames(cluster98) <- "cluster98"
cluster97<- DECIPHER::IdClusters(
  d.align,
  method = "complete",
  cutoff = 0.03)
colnames(cluster97) <- "cluster97"
cluster96<- DECIPHER::IdClusters(
  d.align,
  method = "complete",
  cutoff = 0.04)
colnames(cluster96) <- "cluster96"
cluster95 <- DECIPHER::IdClusters(
  d.align,
  method = "complete",
  cutoff = 0.05)
colnames(cluster95) <- "cluster95"

head(viral_clusters <- cbind(cluster99, cluster98, cluster97, cluster96, cluster95))





dna <- DNAStringSet(c("ACTG", "ACCG"))
dna
DistanceMatrix(dna)
# changing the output type to "dist":
d <- DistanceMatrix(seqs, type="dist")
head(d)



