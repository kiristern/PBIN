#####   EVALUATE ACCURACY #####

#evaluating dada2 accuracy on mock community
rownames(seqtab.nochim)

unqs.mock <- seqtab.nochim["Mock",] #change "Mock" to file name with mock community
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")



##### PHYLOSEQ #####
library(phyloseq)
library(Biostrings)
library(ggplot2)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_names(sample_info_tab))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#visualize alpha-diversity
plot_richness(ps, x="Bloom", measures=c("Shannon", "Simpson")) #need bloom/no bloom data
