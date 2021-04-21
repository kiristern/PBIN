library(dada2)
packageVersion("dada2") # 1.11.5 when this was put together

setwd("~/Desktop/Shapiro_lab/data/data_test_cutadapt/nico")

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
                              rm.phix=FALSE, minLen=175, truncLen=c(210,200))

class(filtered_out) # matrix
dim(filtered_out) # 20 2

filtered_out
