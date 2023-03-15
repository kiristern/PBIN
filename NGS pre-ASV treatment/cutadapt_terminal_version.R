# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2)
packageVersion("dada2") # 1.11.5 when this was put together

setwd("~/Desktop/Shapiro_lab/data/data_test_cutadapt/terminal")
setwd("~/Desktop/Shapiro_lab/data")


list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("samples", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_R1_trimmed.fastq.gz")
# and one with the reverse
reverse_reads <- paste0(samples, "_R2_trimmed.fastq.gz")

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_R1_filtered.fastq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fastq.gz")

plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)
plotQualityProfile(reverse_reads[2:5])

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads,
                              maxN=0, truncQ=2, minLen = 50,
                              rm.phix=FALSE, compress=TRUE, multithread=TRUE)

class(filtered_out) # matrix
dim(filtered_out) # 181 2

filtered_out

#count presence of primers in the sample (should all be 0)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtered_forward_reads), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtered_reverse_reads), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtered_forward_reads), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtered_reverse_reads))

