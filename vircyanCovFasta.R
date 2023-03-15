#select only viral ASVs that have an associated cyanobacteria (to blast)
vircyan_cov <- read.csv("vircyan.cov.csv", header = T, row.names=1)
head(vircyan_cov)

vircyan.top06 <- vircyan_cov %>% 
  filter(weight > 0.6) 
head(vircyan.top06)

unique(vircyan_cov$from)
library(stringr)
vir2ID <- str_replace_all(unique(vircyan_cov$from), "vir_", "")

library("Biostrings")
vir.fasta <- readDNAStringSet("data/ASVs_UPDATED.fa")
seq_name <- names(vir.fasta)[names(vir.fasta) %in% vir2ID]
sequence <- paste(vir.fasta)[names(vir.fasta) %in% seq_name]

#export as fasta file
cat(file="vir2ID.fa", paste(paste0(">",seq_name), 
    sapply(sequence, paste, collapse=""), sep="\n"), sep="\n") 


##### get top 20 #####



# import blastn ID's from ncbi blastn nucleotide
virID <- read.table("data/7VZEVRXA013-Alignment-2.txt", fill=T, sep = "\n")
virname <- virID$V1[grep('^>', virID$V1)]
seqID <- virID$V1[grep('^Sequence ID:', virID$V1)]
length(virname)
length(seqID)
nameID <- as.data.frame(virname, seqID)

blastn.ncbi <- read.csv("7VZEVRXA013-Alignment-HitTable.csv", header = F)
head(blastn.ncbi)
dim(blastn.ncbi)
blastn.ncbi$name <- virname



