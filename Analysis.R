library("phyloseq")
library("vegan")
library("DESeq2") 
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

setwd("~/Desktop/Shapiro_lab/data")
list.files()

rm(list=ls())

count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t") #[ , -c(1:4)]#[,-c(1:$)] removing "blank" samples (col 1:4) from decontaminated count table

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

#meta-data
sample_info_tab <- read.csv("mapping_bloom2_new_corrected.csv", header = TRUE,sep=",")

# Clean data: filter ASV less than 10 occurences
library(dplyr)
# create new col with ASV sum
new_count_tab <- mutate(count_tab, rowSums(count_tab))
# filter only ASVs with more than 10 
new_count_tab <- filter(new_count_tab, rowSums(count_tab)>10)
# remove rowSums col
new_count_tab <- select(new_count_tab, -"rowSums(count_tab)")
#rename sample names to just the date
names(new_count_tab) <- gsub(".*\\.", "", names(new_count_tab))

#reformat date: change sampleID "." to "-"
new_sample_info_tab <- sample_info_tab
rename_x.sampleID <- gsub('[.]', '-', new_sample_info_tab$X.SampleID)
#replace X.SampleID values with rename_x.sampleID
new_sample_info_tab[["X.SampleID"]] <- rename_x.sampleID

#view names
df1 <- as.data.frame(colnames(new_count_tab))
df2 <- as.data.frame(new_sample_info_tab$X.SampleID)

# #  set a color column to be of type "character", which helps for plotting later
# sample_info_tab$colour <- as.character(sample_info_tab$colour)
# head(sample_info_tab) # to take a peek

#see which data is the same 
library(arsenal)
summary(comparedf(df1, df2))

##### BETA DIVERSITY #####
# First generate exploratory visualizations (like ordinations and hierarchical clusterings) to see how samples relate and check for problems (ex. batch effect)
# Here, using Eclidean distances to generate exploratory visualization of samples
# must first normalize data using variance stabilizing transformation (recommended against turning counts into proportions or subsampling each sample to the lowest sample's depth)

# make a DESeq2 object. Need to incl "colData" and "design" for later processing
deseq_counts <- DESeqDataSetFromMatrix(new_count_tab, colData = new_sample_info_tab, design = ~type) 

```
count to phyloseq table
calculate Shannon for all samples, if it's below, then that's the cutoff for bloom no bloom
```

##### ALPHA DIVERSITY #####

# Rarefaction curves: only useful for visualization (cannot extrapolate total richness or anything else)
rarecurve(t(new_count_tab), step=100, lwd=2, ylab="ASVs", label=F) #t() transpose because vegan expects rows to be samples and observations (ASV) to be cols
# adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(count_tab)))))

