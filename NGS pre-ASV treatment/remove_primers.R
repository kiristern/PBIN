##https://benjjneb.github.io/dada2/ITS_workflow.html 
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#load test data
path <- "~/Desktop/Shapiro_lab/data/data_test_cutadapt/R"
#run on all data
path <- "~/Desktop/Shapiro_lab/data"
list.files(path)

##get lists of the fwd and rev fastq files in matched order and parse out sample name
# Forward and reverse fastq filenames have format: XXX.001.SAMPLENAME.DATE_RX.fastq.gz 
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

#Inspect read quality profiles
#visualize forward reads
plotQualityProfile(fnFs[5:8])
#visualize rev reads
plotQualityProfile(fnRs[5:8])

##### IDENTIFY PRIMERS ####

# #primers
# g20_F-CS1 <- "ACACTGACGACATGGTTCTACAGTAGWATTTTCTACATTGAYGTTGG"
# g20_R-CS2 <- "TACGGTAGCAGAGACTTGGTCTGGTARCCAGAAATCYTCMAGCAT"
# #CS1-2 TAGs.
# CS1 <- "ACACTGACGACATGGTTCTACA"
# CS2 <- "TACGGTAGCAGAGACTTGGTCT"
#primers sequences without the CS1 tags
FWD <- "GTAGWATTTTCTACATTGAYGTTGG"
REV <- "GGTARCCAGAAATCYTCMAGCAT"

#verify presence and orientation of primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#remove sequences only with ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filteredN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
filteredN
#write.csv(filteredN, "filtN.csv")

#count number of times the primers appear in the fwd and rev read, while considering all orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN))

##look only at first sample -- assuming all files were created using the same library prep
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

##### REMOVE PRIMERS #####

cutadapt <- "/Users/kiristern/miniconda3/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R

#create new folder for trimmed files
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt ("cutadapt -h" for list of parameters)
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-j", 0, #work in parallel, 0 specifies on all available cores
                             # "-m", 200, #min length
                             # "-M", 300, #max length
                             "--discard-untrimmed",
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#count presence of primers in the first cuatadapt sample (should all be 0)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut))

#read names of cutadapt-ed FASTQ files and get matched lists of fwd and rev fastq files
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

# #Extract sample names, assuming filenames have format:
# get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
# sample.names <- unname(sapply(cutFs, get.sample.name))
# head(sample.names)
##better sample.names extract##
sample.names <- gsub(".*001.(.*)\\_R.*","\\1", cutFs, perl=T) 
  #.*001. get rid of everything before and incl 001
    #without * will only remove .001.
  #(.*) keep
  #\\_R.* remove everything after and incl _R
    #without * will just remove _R and keep fastq.gz
sample.names

full.sample.name <- gsub(".*cutadapt.(.*)_R.*","\\1", cutFs, perl=T)
samples <- full.sample.name
##alternative
#samples <- scan("samples", what="character") #from samples txt in working directory

##### INSPECT QUALITY READS #####

#visualize the quality profiles of the forward reads
plotQualityProfile(cutFs[2:5])
#visualize the quality profile of the reverse reads
plotQualityProfile(cutRs[2:5])



##### FILTER AND TRIM #####

#assign filenames for the output of the filtered reads to be stored as fastq.gz files
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

##parameters explained: https://astrobiomike.github.io/amplicon/dada2_workflow_ex 
  #phi.X=FALSE because using viral sequences
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN=0, truncQ=2, minLen=50, rm.phix=FALSE,
                     compress=TRUE, multithread=TRUE)
out
#write.csv(out, "reads.csv")

#count presence of primers in the sample (should all be 0)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs))

#visualize the quality profiles of the forward reads
plotQualityProfile(filtFs[2:5])
#visualize the quality profile of the reverse reads
plotQualityProfile(filtRs[2:5])

# # Place filtered files in filtered/ subdirectory ALTERNATE SAVING
# filtFs2 <- file.path(path.cut, "filtered2", paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs2 <- file.path(path.cut, "filtered2", paste0(sample.names, "_R_filt.fastq.gz"))
# names(filtFs2) <- sample.names
# names(filtRs2) <- sample.names
# out2 <- filterAndTrim(fnFs, filtFs2, fnRs, filtRs2, truncLen=c(210,200),
#                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=FALSE, minLen=175,
#                      compress=TRUE, multithread=TRUE)
# head(out2)



##### LEARN THE ERROR RATES #####

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#visualize the estimated error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)



#####   DEREPLICATE IDENTICAL READS #####

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##### DEREPLICATION SECOND METHOD #####
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# #process only one sequence and assign a number to it (how many times it shows up instead of processing all of them)
# derep_forward <- derepFastq(filtFs, verbose=TRUE)
# names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
# derep_reverse <- derepFastq(filtRs, verbose=TRUE)
# names(derep_reverse) <- samples


##### SAMPLE INFERENCE #####

#core sample inference algorithm (https://www.nature.com/articles/nmeth.3869#methods) is applied to derepl data
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#inspecting returned dada-class object
dadaFs[[1]]


##### MERGE PAIRED READS #####

#merge fwd and rev reads to obtain full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#inspect the merger data.frame from the first sample
head(mergers[[1]])


##### MERGING FWD AND REV READS #####
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
#                                derep_reverse, trimOverhang=TRUE, minOverlap=170)
# 
# # this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
# class(merged_amplicons) # list
# length(merged_amplicons) # 20 elements in this list, one for each of our samples
# names(merged_amplicons) # the names() function gives us the name of each element of the list 
# 
# class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe
# 
# names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names
# # "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"



##### CONSTRUCT ASV TABLE #####

#generating a count table
seqtab = makeSequenceTable(mergers)
class(seqtab)
dim(seqtab)

# Inspect distribution of sequence lengths
#rows = the samples
#cols = sequence variants (ASV)
table(nchar(getSequences(seqtab)))



##### REMOVE CHIMERAS #####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#frequencyof chimeric sequences (most reads should remain)
sum(seqtab.nochim)/sum(seqtab)

#inspect distributuion of sequence lengths
table(nchar(getSequences(seqtab.nochim)))



##### TRACK READS THRU THE PIPELINE #####

#inspect the number of reads that made it thru each step to verify everything worked
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), 
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#make table 2
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab



##### EXTRACTING FROM DADA2 #####

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# # tax table:
# asv_tax <- taxa
# row.names(asv_tax) <- sub(">", "", asv_headers)
# write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)



##### REMOVING LIKELY CONTAMINANTS #####
library(decontam)

colnames(asv_tab) # our blanks are the first 4 of 20 samples in this case
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 16))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 6 as contaminants

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

asv_tax[row.names(asv_tax) %in% contam_asvs, ]

###BLAST IN TERMINAL###

#remove from primary outputs and create new files
# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)



