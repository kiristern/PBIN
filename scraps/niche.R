####### clusters ASVs to form new OTUs #######
#https://github.com/adriaaulaICM/bbmo_niche_sea/blob/master/src/analysis/nucdist_otuclustering.R
#https://github.com/benjjneb/dada2/issues/947
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


