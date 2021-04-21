#https://microbiome.github.io/tutorials/

######## DIVERSITY ###########

library(microbiome)
library(knitr)
data(dietswap)
pseq <- dietswap

#A comprehensive list of global indicators of the ecosystem state can be obtained as follows. 
#This includes various measures of richness, evenness, diversity, dominance, and rarity with default parameters. 
#See the individual functions for more options regarding parameter tuning.
tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))

#Alpha diversity: This returns a table with selected diversity indicators.
tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))

#Richness: returns observed richness with given detection threshold(s)
tab <- richness(pseq)
kable(head(tab))

#The dominance index refers to the abundance of the most abundant species. 
#Various dominance indices are available (see the function help for a list of options).
# Absolute abundances for the single most abundant taxa in each sample
tab <- dominance(pseq, index = "all")
kable(head(tab))

#Also have a function to list the dominating (most abundant) taxa in each sample.
dominant(pseq)

#The rarity indices quantify the concentration of rare or low abundance taxa. 
#Various rarity indices are available (see the function help for a list of options).
tab <- rarity(pseq, index = "all")
kable(head(tab))

#The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5)
tab <- coverage(pseq, threshold = 0.5)
kable(head(tab))

#The core_abundance function refers to the relative proportion of the core species. 
#Non-core abundance provides the complement (1-x; see rare_abundance).
tab <- core_abundance(pseq, detection = .1/100, prevalence = 50/100)

#Gini index is a common measure for inequality in economical income. 
#The inverse gini index (1/x) can also be used as a community diversity measure.
tab <- inequality(pseq)

#Various evenness measures are also available.
tab <- evenness(pseq, "all")
kable(head(tab))


######## Beta diversity and microbiome divergence #########
data(peerj32)
pseq <- peerj32$phyloseq


#intra-individual divergence
library(tidyverse)
betas <- list()
groups <- as.character(unique(meta(pseq)$group))
for (g in groups) {
  
  df <- subset(meta(pseq), group == g)
  beta <- c()
  
  for (subj in df$subject) {
    # Pick the samples for this subject
    dfs <- subset(df, subject == subj)
    # Check that the subject has two time points
    if (nrow(dfs) == 2) {
      s <- as.character(dfs$sample)
      # Here with just two samples we can calculate the
      # beta diversity directly
      beta[[subj]] <- divergence(abundances(pseq)[, s[[1]]],abundances(pseq)[, s[[2]]],method = "bray")
    }
  }
  betas[[g]] <- beta
}
# boxplot
df <- as.data.frame(unlist(betas))
s<- rownames(df)
si<- as.data.frame(s)
si<- separate(si, s, into = c('names','s'))
df1<- bind_cols(df, si)
rownames(df1)<- df1$s ; df1$s<- NULL

p<- ggplot(df1, aes(x = names, y = `unlist(betas)`))+ geom_boxplot() + ylab('') + xlab('')

plot(p) 


# Divergence within individual over time
# Community divergence within individual often tends to increase over time with respect to the baseline sample.

library(MicrobeDS)
library(microbiome)
library(dplyr)
library(vegan)
data(MovingPictures)

# Pick the metadata for this subject and sort the
# samples by time

# Pick the data and modify variable names
pseq <- MovingPictures
s <- "F4" # Selected subject
b <- "UBERON:feces" # Selected body site

# Let us pick a subset
pseq <- subset_samples(MovingPictures, host_subject_id == s & body_site == b) 

# Rename variables
sample_data(pseq)$subject <- sample_data(pseq)$host_subject_id
sample_data(pseq)$sample <- sample_data(pseq)$X.SampleID

# Tidy up the time point information (convert from dates to days)
sample_data(pseq)$time <- as.numeric(as.Date(gsub(" 0:00", "", as.character(sample_data(pseq)$collection_timestamp)), "%m/%d/%Y") - as.Date("10/21/08", "%m/%d/%Y"))

# Order the entries by time
df <- meta(pseq) %>% arrange(time)

# Calculate the beta diversity between each time point and
# the baseline (first) time point
beta <- c() # Baseline similarity
s0 <- subset(df, time == 0)$sample
# Let us transform to relative abundance for Bray-Curtis calculations
a <-microbiome::abundances(microbiome::transform(pseq, "compositional")) 
for (tp in df$time[-1]) {
  # Pick the samples for this subject
  # If the same time point has more than one sample,
  # pick one at random
  st <- sample(subset(df, time == tp)$sample, 1)
  # Beta diversity between the current time point and baseline
  b <- vegdist(rbind(a[, s0], a[, st]), method = "bray")
  # Add to the list
  beta <- rbind(beta, c(tp, b))
}
colnames(beta) <- c("time", "beta")
beta <- as.data.frame(beta)

theme_set(theme_bw(20))
library(ggplot2)
p <- ggplot(beta, aes(x = time, y = beta)) +
  geom_point() +
  geom_line() +
  geom_smooth() +
  labs(x = "Time (Days)", y = "Beta diversity (Bray-Curtis)")
print(p)


#Inter-individual divergence / spread
#Divergence within a sample set quantifies the overall heterogeneity in community composition across samples or individuals. 
#This is sometimes quantified as the average dissimilarity of each sample from the group mean; the dissimilarity can be 
#quantified by beta diversity

#Calculate divergences within the LGG (probiotic) and Placebo groups with respect to the median profile within each group.

pseq <- peerj32$phyloseq

b.pla <- divergence(subset_samples(pseq, group == "Placebo"),
                    apply(abundances(subset_samples(pseq, group == "Placebo")), 1, median))

b.lgg <- divergence(subset_samples(pseq, group == "LGG"),
                    apply(abundances(subset_samples(pseq, group == "LGG")), 1, median))

#The group with larger values has a more heterogeneous community composition.
library(reshape)
l<- list(b.pla, b.lgg)
df<- melt(l)
df$L1[df$L1 == '1']<- 'placebo'
df$L1[df$L1 == '2']<- 'LGG'

df$L1<- factor(df$L1, levels = c('placebo','LGG'))

p<- ggplot(df, aes(x = L1, y = value)) + geom_boxplot()+ xlab('')

plot(p)



############# TEMPORAL MICROBIOTA TRAJECTORY ################
#plot trajactory of a microbiota through time

#### To install this pkg ####
# install.packages("devtools")
# devtools::install_github("microsud/jeevanuDB")
############################

library(microbiome)
library(jeevanuDB) # external database pkg for microbiome pkg with test data
library(dplyr)
library(ggplot2)
library(viridis)
library(knitr)
# Example data
data("moving_pictures")
# Rename
ps <- moving_pictures

#Check data for which and how many samples are present each subject.
kable(table(meta(ps)$host_subject_id, meta(ps)$sample_type))

#Use only the stool (gut) microbiota data.
ps.gut <- subset_samples(ps, sample_type == "stool")
taxa_names(ps.gut) <- paste0("ASV-", seq(ntaxa(ps.gut))) # rename sequences to ids

# remove asvs which are zero in all of these samples
ps.gut <- prune_taxa(taxa_sums(ps.gut) > 0, ps.gut)

# remove samples with less than 500 reads Note: this is user choice 
# here we just show example
ps.gut <- prune_samples(sample_sums(ps.gut) > 500, ps.gut)

# Covnert to relative abundances
ps.gut.rel <- microbiome::transform(ps.gut, "compositional")

#Using phyloseq we do ordination analysis.
# Ordination object
ps.ord <- ordinate(ps.gut.rel, "PCoA")

# Ordination object plus all metadata, note: we use justDF=T. This will not return a plot but a data.frame
ordip <- plot_ordination(ps.gut.rel, ps.ord, justDF = T)

#start to make a custom visualization
# Get axis 1 and 2 variation
evals1 <- round(ps.ord$values$Eigenvalues[1] / sum(ps.ord$values$Eigenvalues) * 100, 2)
evals2 <- round(ps.ord$values$Eigenvalues[2] / sum(ps.ord$values$Eigenvalues) * 100, 2)

# theme_set(theme_bw(14))
# set colors
subject_cols <- c(F4 = "#457b9d", M3 = "#e63946")

#Add trajectory for the subject of interest. Here, we randomnly choose subject F4
# choose data for subject F4
dfs <- subset(ordip, host_subject_id == "F4")

# arrange according to sampling time. Sort with increasing time
dfs <- dfs %>%
  arrange(days_since_experiment_start)


#Initiate plotting
# use the ordip
# first step get blank plot
p <- ggplot(ordip, aes(x = Axis.1, y = Axis.2))

# add path (lines) join only those samples that are from F4
p2 <- p +
  geom_path(
    data = dfs, alpha = 0.7,
    arrow = arrow(
      angle = 15, length = unit(0.1, "inches"),
      ends = "last", type = "closed"
    )
  ) +
  # now add a layer of points 
  geom_point(aes(color = host_subject_id), alpha = 0.6, size = 3) +
  scale_color_manual("Subject", values = subject_cols) +
  xlab(paste("PCoA 1 (", evals1, "%)", sep = "")) +
  ylab(paste("PCoA 2 (", evals2, "%)", sep = "")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
# coord_fixed(sqrt(evals[2] / evals[1]))

# Print figure
print(p2)

#Alternatively we can just focus on one subject.
# subset data for only M3
ps.gut.rel.m3 <- subset_samples(ps.gut.rel, host_subject_id == "M3")
# remove asvs which are zero in all of these samples
ps.gut.rel.m3 <- prune_taxa(taxa_sums(ps.gut.rel.m3) > 0, ps.gut.rel.m3)

ps.ord.m3 <- ordinate(ps.gut.rel.m3, "PCoA")

ordip.m3 <- plot_ordination(ps.gut.rel.m3, ps.ord.m3, justDF = T)

# Get axis 1 and 2 variation
evals1 <- round(ordip.m3$values$Eigenvalues[1] / sum(ordip.m3$values$Eigenvalues) * 100, 2)
evals2 <- round(ordip.m3$values$Eigenvalues[2] / sum(ordip.m3$values$Eigenvalues) * 100, 2)
# arrange according to sampling time
ordip.m3 <- ordip.m3 %>%
  arrange(days_since_experiment_start) # important to arrange the time

# Visualize
# blank plot initiate
p1 <- ggplot(ordip.m3, aes(x = Axis.1, y = Axis.2))
# add layers
p3 <- p1 +
  # add arrows with geom_path
  # geom_path(alpha = 0.5, arrow = arrow(
  #   angle = 30, length = unit(0.1, "inches"),
  #   ends = "last", type = "closed"
  # )) +
  # add points
  geom_point(aes(color = days_since_experiment_start), size = 3) +
  # add gradient colors
  scale_color_viridis("Days from first sampling") +
  # add x and y labels
  xlab(paste("PCoA 1 (", evals1, "%)", sep = "")) +
  ylab(paste("PCoA 2 (", evals2, "%)", sep = "")) +
  theme_bw() +
  # remove grids in the plot
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p3)




################# MICROBIOME COMPOSITION ###############
# Example data
library(microbiome)
library(dplyr)
data(dietswap)

# Make sure we use functions from correct package
transform <- microbiome::transform

# Merge rare taxa to speed up examples
pseq <- transform(dietswap, "compositional")
pseq <- aggregate_rare(pseq, level = "Genus", detection = 1/100, prevalence = 50/100)

# Pick sample subset
library(phyloseq)
pseq2 <- subset_samples(pseq, group == "DI" & nationality == "AFR" & timepoint.within.group == 1)

# Normal western adults
data(atlas1006)
pseq3 <- atlas1006 %>%
  subset_samples(DNA_extraction_method == "r") %>%
  aggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")


#Composition barplots


# Try another theme
# from https://github.com/hrbrmstr/hrbrthemes
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
theme_set(theme_bw(21))
p <- pseq3 %>%
  plot_composition(sample.sort = "Firmicutes", otu.sort = "abundance") +
  # Set custom colors
  scale_fill_manual(values = default_colors("Phylum")[taxa(pseq3)]) +
  scale_y_continuous(label = scales::percent)

print(p)


# Limit the analysis on core taxa and specific sample group
p <- plot_composition(pseq2,
                      taxonomic.level = "Genus",
                      sample.sort = "nationality",
                      x.label = "nationality") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data",
       subtitle = "Subtitle",
       caption = "Caption text.") + 
  theme_ipsum(grid="Y")
print(p)  

# Averaged by group
p <- plot_composition(pseq2,
                      average_by = "bmi_group", transform = "compositional")
print(p)



p <- NULL


# composition heatmaps
#Heatmap for CLR-transformed abundances, with samples and OTUs sorted with the neatmap method:
p <- plot_composition(microbiome::transform(pseq, "compositional"),
                      plot.type = "heatmap",
                      sample.sort = "neatmap", otu.sort = "neatmap")
print(p)



#plot taxa prevalence
#overview of OTU prevalences alongwith their taxonomic affiliations. 
#This will aid in checking if you filter OTUs based on prevalence, then what taxonomic affliations will be lost.

data(atlas1006)

# Use sample and taxa subset to speed up example
p0 <- subset_samples(atlas1006, DNA_extraction_method == "r")

# Define detection and prevalence thresholds to filter out rare taxa
p0 <- core(p0, detection = 0.1/100, prevalence = 1/100)

# For the available taxonomic levels
plot_taxa_prevalence(p0, "Phylum", detection = 0.1/100)




# Amplicon data

# Example data
library(microbiome)
# Try another theme
# from https://github.com/hrbrmstr/hrbrthemes
# you can install these if you don't have it already.
# devtools::install_github("hrbrmstr/hrbrthemes")
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(dplyr)
library(jeevanuDB)

ps1 <- emp_human
colnames(tax_table(ps1))


#can see the taxonomic classification is just lablled as “Rank1” … “Rank7”. We need to change this to proper designation and also do some formatting of the data. This can be a useful example for understanding simple file processing in R.
#In case you see the taxonomic classification is just lablled as “Rank1” … “Rank7” we can change it as follows
# First change the column names of the taxonomy table in phyloseq to following:

colnames(tax_table(ps1)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species" )

tax_table(ps1)[tax_table(ps1)[,"Domain"]== "NA", "Domain" ] <- "Unidentified_Domain"

tax_table(ps1)[tax_table(ps1)[,"Phylum"]== "p__", "Phylum" ] <- "p__Unidentified_Phylum"


#plot counts abundance
# try at Family level.
library(phyloseq)
# merge at family level.
# check how many samples are there

# Use only saliva samples 

ps1.saliva <- subset_samples(ps1, env_material == "saliva")

total_samples <- nsamples(ps1.saliva)
ps1.saliva.pruned <- prune_taxa(taxa_sums(ps1.saliva) >0, ps1.saliva)
# merge all taxa that are detected rare
pseq.fam <- aggregate_rare(ps1.saliva.pruned, level="Family", detection = 50, prevalence = 25/total_samples)

p.fam <- plot_composition(pseq.fam, sample.sort = NULL, 
                          otu.sort = NULL,
                          x.label = "empo_3", # sample type
                          plot.type = "barplot", 
                          verbose = FALSE) + 
  theme_bw() + scale_fill_brewer("Family", palette = "Paired")
# we can rotate x axis labels 
print(p.fam + theme(axis.text.x = element_text(angle = 90)))



#plot relative abundance

pseq.famrel <- microbiome::transform(pseq.fam, "compositional")

p.famrel <- plot_composition(pseq.famrel, sample.sort = NULL, otu.sort = NULL,
                             x.label = "empo_3", plot.type = "barplot", verbose = FALSE)

print(p.famrel)

# further improvements can be done as follows  

p.famrel <- plot_composition(pseq.famrel, 
                             sample.sort = NULL, 
                             otu.sort = NULL, 
                             x.label = "empo_3", 
                             plot.type = "barplot", 
                             verbose = FALSE) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") + 
  scale_fill_brewer("Family", palette = "Paired")

print(p.famrel)



# Average by group
# Use all samples 
ps1 <- emp_human
# get relative abudance
ps1.rel <- microbiome::transform(ps1, "compositional")
ps1.fam.rel <-aggregate_rare(ps1.rel, level = "Family", detection = 0.005, prevalence = 0.5)

p <- plot_composition(ps1.fam.rel,
                      average_by = "empo_3") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Samples", 
       y = "Relative abundance",
       title = "Relative abundance data", 
       subtitle = "Subtitle",
       caption = "Caption text.") 
print(p + scale_fill_brewer("Family", palette = "Paired") + theme_bw())


#heatmap composition

# Use all samples 
ps1 <- emp_human
ps1.rel <-aggregate_rare(ps1, level = "Family", detection = 10, prevalence = 0.5)

pseq.famlog <- microbiome::transform(ps1.rel, "log10")

p.famrel.heatmap <- plot_composition(pseq.famlog, 
                                     sample.sort = NULL, 
                                     otu.sort = NULL, 
                                     x.label = "empo_3", 
                                     plot.type = "heatmap", 
                                     verbose = FALSE)

print(p.famrel.heatmap)


#plot core taxa time trajectory

library(dplyr)
# select core
ps <- moving_pictures
table(meta(ps)$sample_type, meta(ps)$host_subject_id)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
# Filter the data to include only gut samples from M3 subject
ps.m3 <- subset_samples(ps, sample_type == "stool" & host_subject_id == "M3") 
#print(ps.m3)

ps.m3.rel <- microbiome::transform(ps.m3, "compositional")
pseq.core <- core(ps.m3.rel, detection = 0.001, prevalence = .95)

ps.stool.df <- psmelt(pseq.core)
#head(ps.stool.df)

# add genus name to ASVid
ps.stool.df <- ps.stool.df %>% 
  mutate(asv_gen= paste0(OTU, "-",Genus))

ps.stool.rel.plot <- ggplot(ps.stool.df) + 
  geom_line(aes(days_since_experiment_start, 
                Abundance, color = asv_gen)) +
  theme_bw() + 
  theme(legend.position="top") + 
  xlab("Days since experiment start") + 
  ylab("Relative abundance") + 
  scale_color_brewer("Core ASVs",palette = "Paired") +
  guides(col = guide_legend(ncol = 3, nrow = 3))

ps.stool.rel.plot 


#Highlight only one ASVs of interest.
ps.highlight.plot <- ggplot(ps.stool.df) + 
  geom_line(aes(days_since_experiment_start, 
                Abundance), color="grey80") 
# pick only data for ASV996-g__Faecalibacterium

asv996 <- subset(ps.stool.df, asv_gen =="ASV996-g__Faecalibacterium")

ps.highlight.plot <- ps.highlight.plot + 
  geom_line(data= asv996,aes(x=days_since_experiment_start, 
                             y=Abundance, color=asv_gen)) +
  theme_bw() + 
  theme(legend.position="top") + 
  xlab("Days since experiment start") + 
  ylab("Relative abundance") + 
  scale_color_manual("Core ASVs",values="brown3") +
  guides(col = guide_legend(ncol = 3, nrow = 3))

ps.highlight.plot 





############### CORE MICROBIOME ##############
#See also related functions for the analysis of rare and variable taxa:
#rare_members; rare_abundance; rare_members; rare_abundance; low_abundance.


#HITChip Data
# Load data
library(microbiome)
data(peerj32)

# Rename the data
pseq <- peerj32$phyloseq

# Calculate compositional version of the data
# (relative abundances)
pseq.rel <- microbiome::transform(pseq, "compositional")


#Prevalence of taxonomic groups
#Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE))

#Absolute population frequencies (sample count):
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE, count = TRUE))


#Core microbiota analysis
#If you only need the names of the core taxa, do as follows. This returns the taxa that exceed the given prevalence and detection thresholds.
core.taxa.standard <- core_members(pseq.rel, detection = 0, prevalence = 50/100)

#full phyloseq object of the core microbiota
pseq.core <- core(pseq.rel, detection = 0, prevalence = .5)

# collapse the rare taxa into an “Other” category
pseq.core2 <- aggregate_rare(pseq.rel, "Genus", detection = 0, prevalence = .5)

#Retrieving the core taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core)


#Core abundance and diversity
# Total core abundance in each sample (sum of abundances of the core members):
core.abundance <- sample_sums(core(pseq.rel, detection = .01, prevalence = .95))

#Core visualization
#Core line plots

#Determine core microbiota across various abundance/prevalence thresholds with the blanket analysis based on various signal and prevalences
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


plot_core(pseq.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")


#Core heatmaps

# Core with compositionals:
library(RColorBrewer)
library(reshape)
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))
#pseq.rel<- microbiome::transform(pseq, 'compositional')
p <- plot_core(pseq.rel, plot.type = "heatmap", colours = gray,
               prevalences = prevalences, detections = detections) +
  labs(x = "Detection Threshold (Relative Abundance (%))")

print(p)



# Core with absolute counts and horizontal view:
# and minimum population prevalence (given as percentage)
detections <- 10^seq(log10(1), log10(max(abundances(pseq))/10), length = 10)

p <- plot_core(pseq, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = TRUE) 
print(p)



# Core Microbiota using Amplicon data
#Make phyloseq object
#This tutorial is useful for analysis of output files from (Mothur), (QIIME or QIIME2) or any tool that gives a biom file as output. 
#There is also a simple way to read comma seperated (*.csv) files.


library(microbiome)

otu.file <-
  system.file("extdata/qiita1629_otu_table.csv",
              package='microbiome')

tax.file <- system.file("extdata/qiita1629_taxonomy_table.csv",
                        package='microbiome')

meta.file <- system.file("extdata/qiita1629_mapping_subset.csv",
                         package='microbiome')

pseq.csv <- read_phyloseq(
  otu.file=otu.file, 
  taxonomy.file=tax.file, 
  metadata.file=meta.file, type = "simple")



# Read the biom file
biom.file <- 
  system.file("extdata/qiita1629.biom", 
              package = "microbiome")

# Read the mapping/metadata file
meta.file <- 
  system.file("extdata/qiita1629_mapping.csv", 
              package = "microbiome")
# Make phyloseq object
pseq.biom <- read_phyloseq(otu.file = biom.file, 
                           metadata.file = meta.file, 
                           taxonomy.file = NULL, type = "biom")


#Mothur shared OTUs and Consensus Taxonomy:
otu.file <- system.file(
  "extdata/Baxter_FITs_Microbiome_2016_fit.final.tx.1.subsample.shared",
  package='microbiome')

tax.file <- system.file(
  "extdata/Baxter_FITs_Microbiome_2016_fit.final.tx.1.cons.taxonomy",
  package='microbiome')

meta.file <- system.file(
  "extdata/Baxter_FITs_Microbiome_2016_mapping.csv",
  package='microbiome')

pseq.mothur <- read_phyloseq(otu.file=otu.file,
                             taxonomy.file =tax.file,
                             metadata.file=meta.file, type = "mothur")
print(pseq.mothur)


# Core microbiota analysis

# check the data 
library(jeevanuDB)
ps <- moving_pictures
table(meta(ps)$sample_type, meta(ps)$host_subject_id)

# Filter the data to include only gut samples from M3 subject
ps.m3 <- subset_samples(ps, sample_type == "stool" & host_subject_id == "M3") 
print(ps.m3)

# keep only taxa with positive sums
ps.m3 <- prune_taxa(taxa_sums(ps.m3) > 0, ps.m3)
print(ps.m3)


# Calculate compositional version of the data
# (relative abundances)
ps.m3.rel <- microbiome::transform(ps.m3, "compositional")

#Output of deblur/dada2 will most likely have seqs as rownames instead of OTU ids or taxa names
taxa_names(ps.m3.rel)[1:3]

#We can change it to ASVIDs
ibrary(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(ps.m3.rel))
names(dna) <- taxa_names(ps.m3.rel)
ps.m3.rel <- merge_phyloseq(ps.m3.rel, dna)
taxa_names(ps.m3.rel) <- paste0("ASV", seq(ntaxa(ps.m3.rel)))
# now check again
taxa_names(ps.m3.rel)[1:3]



#If you only need the names of the core taxa, do as follows. 
#This returns the taxa that exceed the given prevalence and detection thresholds.
core.taxa.standard <- core_members(ps.m3.rel, detection = 0.0001, prevalence = 50/100)
core.taxa.standard

#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(ps.m3.rel, detection = 0.0001, prevalence = .5)

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core)
class(core.taxa)


# get the taxonomy data
tax.mat <- tax_table(pseq.core)
tax.df <- as.data.frame(tax.mat)

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))



#Core visualization
#Determine core microbiota across various abundance/prevalence thresholds with the blanket analysis based on various signal and prevalences.
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(ps.m3.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))
p1 <- plot_core(ps.m3.rel, 
                plot.type = "heatmap", 
                colours = gray,
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1

library(viridis)
print(p1 + scale_fill_viridis())
#As it can be seen, we see only OTu IDs and this may not be useful to interpret the data. 
#We need to repreoccess this figure to include taxonomic information. We can do this as follows:
library(RColorBrewer)
library(knitr)
# get the data used for plotting 
df <- p1$data 

# get the list of OTUs
list <- df$Taxa 

# check the OTU ids
# print(list) 

# get the taxonomy data
tax <- as.data.frame(tax_table(ps.m3.rel))

# add the ASVs to last column
tax$ASV <- rownames(tax)

# select taxonomy of only 
# those OTUs that are used in the plot
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 

# head(tax2)

# We will merege all the column into one except the Doamin as all is bacteria in this case
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), sep = "_;", remove = TRUE)

tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)

# add this new information into the plot data df

df$Taxa <- tax.unit$Taxa_level

# you can see now we have the taxonomic information
knitr::kable(head(df))

# replace the data in the plot object
p1$data <- df

plot(p1 + theme(axis.text.y = element_text(face="italic")))


#Genus level
ps.m3.rel.gen <- aggregate_taxa(ps.m3.rel, "Genus")

library(RColorBrewer)
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

p1 <- plot_core(ps.m3.rel.gen, 
                plot.type = "heatmap", 
                colours = rev(brewer.pal(5, "RdBu")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1




####### Dirichlet Multinomial Mixtures #############
#Community typing with Dirichlet Multinomial Mixtures
#DMM is a probabilistic method for community typing (or clustering) of microbial community profiling data. 
#It is an infinite mixture model, which means that the method can infer the optimal number of community types. 
#Note that the number of community types is likely to grow with data size.

library(microbiome)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DirichletMultinomial")
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
# Load example data
data(dietswap)
pseq <- dietswap

# To speed up, only consider the core taxa
# that are prevalent at 0.1% relative abundance in 50% of the samples
# (note that this is not strictly correct as information is
# being discarded; one alternative would be to aggregate rare taxa)
pseq.comp <- microbiome::transform(pseq, "compositional")
taxa <- core_members(pseq.comp, detection = 0.1/100, prevalence = 50/100)
pseq <- prune_taxa(taxa, pseq)

# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- abundances(pseq)
count <- as.matrix(t(dat))

#Fit the DMM model. Let us set the maximum allowed number of community types to 3 to speed up the example.
fit <- lapply(1:3, dmn, count = count, verbose=TRUE)

#Check model fit with different number of mixture components using standard information criteria
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
#plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)

#pick the optimal model
best <- fit[[which.min(unlist(lplc))]]

#Mixture parameters pi and theta
mixturewt(best)

#Sample-component assignments
ass <- apply(mixture(best), 1, which.max)

#Contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}




###### Dissimilarity-Overlap analysis ##########
library(microbiome)
data(atlas1006)


#Estimate the overlap and dissimilarity quantities:
# Dissimilarity
d <- phyloseq::distance(microbiome::transform(atlas1006, "compositional"), "jsd", parallel=TRUE)

# Overlap
o <- overlap(atlas1006, detection = 0.2/100)

# Compare
dvec <- d[lower.tri(d)]
ovec <- o[lower.tri(o)]

# Assess rough correlation
cc <- cor(dvec, ovec, method = "spearman", use = "pairwise.complete.obs")

# Scatterplot
plot(dvec, ovec, pch = ".", main = paste("Spearman rho", round(cc, 2)), las = 1, xlab = "Dissimilarity (Jensen-Shannon)", ylab = "Overlap")

#p <- data.frame(D = dvec, O = ovec) %>%
#  ggplot(aes(x = D, y = O)) +
#  geom_hex()




##########Microbiome Landscapes###########
#Microbiome Landscaping refers to the analysis and illustration of population frequencies. 
#Typically, these are wrappers based on standard ordination methods

#2D microbiome landscape
library(microbiome)
library(phyloseq)
library(ggplot2)

data(dietswap)
pseq <- dietswap

# Convert to compositional data
pseq.rel <- microbiome::transform(pseq, "compositional")

# Pick core taxa
pseq.core <- core(pseq.rel, detection = 5/100, prevalence = 5/100)
pseq.core <- subset_samples(pseq.core, sex == "female" &
                              bmi_group == "overweight")


#landscape figure

#PCA
# PCA with euclidean distance and CLR transformation
p <- plot_landscape(pseq, method = "PCA", transformation = "clr") +
  labs(title = paste("PCA / CLR"))
print(p)

#PCoA / MDS
# PCoA for compositional data with Bray-Curtis distances
p <- plot_landscape(microbiome::transform(pseq.core, "compositional"),
                    method = "PCoA", distance = "bray") +
  labs(title = paste("PCoA / Compositional / Bray-Curtis"))
print(p)

#t-SNE
p <- plot_landscape(pseq, "t-SNE",
                    distance = "euclidean", transformation = "hellinger") +
  labs(title = paste("t-SNE / Hellinger / Euclidean"))       
print(p)

#NMDS
# Landscape plot directly from phyloseq object
p <- plot_landscape(pseq.core, "NMDS", "bray", col = "nationality") +
  labs(title = paste("NMDS / Bray-Curtis"))       
#For direct access to the ordination coordinates, use the following:
# Project the samples with the given method and dissimilarity measure. 
# Ordinate the data; note that some ordinations are sensitive to random seed
# "quiet" is used to suppress intermediate outputs
set.seed(423542)
x <- pseq.core
quiet(x.ord <- ordinate(x, "NMDS", "bray"))
# Pick the projected data (first two columns + metadata)
proj <- phyloseq::plot_ordination(x, x.ord, justDF=TRUE)
# Rename the projection axes
names(proj)[1:2] <- paste("Comp", 1:2, sep=".")

# Same with a generic data.frame
# (note that random seed will affect the exact ordination)
p <- plot_landscape(proj[, 1:2], col = proj$nationality, legend = T)
print(p)

# Visualize sample names:
ax1 <- names(proj)[[1]]
ax2 <- names(proj)[[2]]
p <- ggplot(aes_string(x = ax1, y = ax2, label = "sample"), data = proj) +
  geom_text(size = 2)
print(p)


#Abundance histograms (one-dimensional landscapes)
#Population densities for Dialister:
# Load libraries
library(microbiome)
library(phyloseq)
pseq <- dietswap

# Visualize population densities for specific taxa
plot_density(pseq, "Dialister") + ggtitle("Absolute abundance")

# Same with log10 compositional abundances
x <- microbiome::transform(pseq, "compositional")
tax <- "Dialister"
plot_density(x, tax, log10 = TRUE) +
  ggtitle("Relative abundance") +
  xlab("Relative abundance (%)")



######## Microbiome stability analysis #######
# Load the example data
set.seed(134)
library(microbiome)
library(dplyr)
data(atlas1006)

# Rename the example data
pseq <- atlas1006

# Focus on specific subset
pseq <- pseq %>% subset_samples(DNA_extraction_method == "r")

# Use relative abundances
pseq <- microbiome::transform(pseq, "compositional")

# Merge rare taxa to speed up examples
pseq <- aggregate_rare(pseq, level = "Genus", detection = .1/100, prevalence = 10/100)

# For cross-sectional analysis, use only the baseline time point:
pseq0 <- baseline(pseq)



#Intermediate stability quantification
#It has been reported that certain microbial groups exhibit bi-stable abundance distributions with distinct peaks at low and high abundances, 
#and an instable intermediate abundance range. Instability at the intermediate abundance range is hence one indicator of bi-stability. 
#Lahti et al. 2014 used straightforward correlation analysis to quantify how the distance from the intermediate abundance region (50% quantile) 
#is associated with the observed shifts between consecutive time points.

# Here we use the (non-baseline) phyloseq object that contains time series.
intermediate.stability <- intermediate_stability(pseq, output = "scores")

#bimodality quantification
#Bimodality of the abundance distribution provides another (indirect) indicator of bistability, although other explanations such as sampling biases etc. should be controlled. 
#Multiple bimodality scores are available.

#Multimodality score using potential analysis with bootstrap
# Bimodality is better estimated from abundances at log scale (such as CLR)
pseq0.clr <- microbiome::transform(pseq0, "clr")

set.seed(4433)
# In practice, it is recommended to use more bootstrap iterations than in this example
bimodality.score <- bimodality(pseq0.clr, method = "potential_analysis",
                               bs.iter = 20, peak.threshold = 10,
                               min.density = 10)

#Comparing bimodality and intermediate stability
#The analysis suggests that bimodal population distribution across individuals is often associated with instable intermediate abundances within individuals. 
#The specific bi-stable groups in the upper left corner were suggested to constitute bistable tipping elements of the human intestinal microbiota

taxa <- taxa(pseq0)
df <- data.frame(group = taxa,
                 intermediate.stability = intermediate.stability[taxa],
                 bimodality = bimodality.score[taxa])

theme_set(theme_bw(20))
library(ggrepel)
p <- ggplot(df,
            aes(x = intermediate.stability, y = bimodality, label = group)) +
  geom_text_repel() +
  geom_point() 
print(p)


#tipping point detection
#Identify potential minima in cross-sectional population data as tipping point candidates.

# Log10 abundance for a selected taxonomic group
# Pick the most bimodal taxa as an example
tax  <- names(which.max(bimodality.score))

# Detect tipping points at log10 abundances
x <- abundances(microbiome::transform(pseq, "clr"))[tax,]

# Bootstrapped potential analysis to identify potential minima
# in practice, use more bootstrap iterations
potential.minima <- potential_analysis(x, bs.iter = 10)$minima

# Same with earlywarnings package (without bootstrap ie. less robust)
# library(earlywarnings)
# res <- livpotential_ews(x)$min.points

# Identify the potential minimum locations as tipping point candidates
tipping.samples <- sapply(potential.minima, function (m) {names(which.min(abs(sort(x) - m)))})
tipping.point <- abundances(pseq)[tax, tipping.samples]
print(tipping.point)

# Illustrate the detected tipping point
plot(density(x), main = tax)
abline(v = x[tipping.samples])




########### Heatmaps for microbiome analysis #############
library(microbiome) # Load libraries
library(phyloseq)
library(dplyr)
library(reshape2)

library(knitr)


data(peerj32)
pseq <- peerj32$phyloseq    # Rename data


# Pick data subset (DI samples from Phylum Bacteroidetes)
pseq2 <- pseq %>%
  subset_taxa(Phylum == "Bacteroidetes") %>%
  subset_samples(group == "LGG")


# Z transformed abundance data
pseqz <- microbiome::transform(pseq2, "Z")


#matrix heatmaps
#Visualize the Z-transformed abundance matrix

# Plot the abundances heatmap
dfm <- melt(abundances(pseqz))
colnames(dfm) <- c("Taxa", "Sample", "value")
heat(dfm, "Taxa", "Sample", "value")

#Find visually appealing order for rows and columns with the Neatmap approach:
# Sort the matrix rows and cols directly
xo <- neat(abundances(pseqz), method = "NMDS", distance = "euclidean") 

# Heatmap visualization
dfm <- melt(xo)
colnames(dfm) <- c("Taxa", "Sample", "value")
heat(dfm, "Taxa", "Sample", "value")

# or use a shortcut to sorting rows (or columns) if just the order was needed 
sorted.rows <- neatsort(abundances(pseqz), "rows", method = "NMDS", distance = "euclidean")
  


#Cross-correlating data sets
#The function returns correlations, raw p-values, and fdr estimates (not strictly proper as the comparisons are not independent). 
#Keep only those elements that have at least only one significant correlation (n.signif):

# Load example data 
otu <- peerj32$microbes 
lipids <- peerj32$lipids 

# Define data sets to cross-correlate
x <- log10(otu) # OTU Log10 (44 samples x 130 genera)
y <- as.matrix(lipids) # Lipids (44 samples x 389 lipids)

# Cross correlate data sets
correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)

# Or, alternatively, the same output is also available in a handy table format
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

kable(head(correlation.table))



#Association heatmaps
#Rearrange the data and plot the heatmap and mark significant correlations with stars to reproduce microbiota-lipidome heatmap
p <- heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 
print(p)


#heatmaps with ggplot2
#The above examples provide handy shortcuts for heatmap visualization. 
#You can also directly modify the ggplot2 routines. This time, let us set q-value threshold also for cell coloring:

# Order the rows and columns with levels argument if needed:
correlation.table$X1 <- factor(correlation.table$X1, levels = unique(as.character(correlation.table$X1)))
correlation.table$X2 <- factor(correlation.table$X2, levels = unique(as.character(correlation.table$X2)))

# Set black-and-white theme
library(ggplot2)
theme_set(theme_bw())

# Pick only the correlations with q<0.05
# Note: this will leave other cells empty
library(dplyr)
subtable <- filter(correlation.table, p.adj < 0.05)

# Arrange the figure
p <- ggplot(subtable, aes(x = X1, y = X2, fill = Correlation))
p <- p + geom_tile() 
p <- p + scale_fill_gradientn("Correlation", 
                              breaks = seq(from = -1, to = 1, by = 0.2), 
                              colours = c("darkblue", "blue", "white", "red", "darkred"), 
                              limits = c(-1,1)) 

# Polish texts
p <- p + theme(axis.text.x=element_text(angle = 90))
p <- p + xlab("") + ylab("")

# Mark the most significant cells with stars
p <- p + geom_text(data = subset(correlation.table, p.adj < 0.02), 
                   aes(x = X1, y = X2, label = "+"), col = "white", size = 5)

# Plot
print(p)


#heatmap with text
#For detailed information, might be handy to print the actual values on top of the heatmap:

theme_set(theme_bw(20))
df <- correlation.table
p <- ggplot(df, aes(X1, X2, group=X2)) 
p <- p + geom_tile(aes(fill = Correlation)) 
p <- p + geom_text(aes(fill = df$Correlation, label = round(df$Correlation, 1)), size = 2) 
p <- p + scale_fill_gradientn("Correlation", 
                              breaks = seq(from = -1, to = 1,  by = 0.25), 
                              colours = c("blue", "white", "red"), 
                              limits = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "", y = "")
print(p)



#ggcorr
library(GGally)
ggcorr(x[, 1:10], method = c("pairwise", "spearman"), nbreaks = 20, hjust = 0.75)
ggcorr(x[, 1:10], method = c("pairwise", "spearman"), nbreaks = 20, geom = "circle")
ggcorr(x[, 1:10], method = c("pairwise", "spearman"), nbreaks = 20, label = TRUE, label_alpha = TRUE)
ggcorr(data = NULL, cor_matrix = cor(x[, 1:10], use = "everything"), low = "steelblue", mid = "white", high = "darkred", midpoint = 0)




######## ORDINATION ANALYSIS ###########
library(microbiome)
library(phyloseq)
library(ggplot2)
data(dietswap)
pseq <- dietswap

# Convert to compositional data
pseq.rel <- microbiome::transform(pseq, "compositional")

# Pick core taxa with with the given prevalence and detection limits
pseq.core <- core(pseq.rel, detection = .1/100, prevalence = 90/100)

# Use relative abundances for the core
pseq.core <- microbiome::transform(pseq.core, "compositional")


#sample ordination
#Project the samples with the given method and dissimilarity measure.

# Ordinate the data
set.seed(4235421)
# proj <- get_ordination(pseq, "MDS", "bray")
ord <- ordinate(pseq, "MDS", "bray")


#Multidimensional scaling (MDS / PCoA)
plot_ordination(pseq, ord, color = "nationality") +
  geom_point(size = 5)

#Canonical correspondence analysis (CCA)
# With samples
pseq.cca <- ordinate(pseq, "CCA")
p <- plot_ordination(pseq, pseq.cca,
                     type = "samples", color = "nationality")
p <- p + geom_point(size = 4)
print(p)

# With taxa:
p <- plot_ordination(pseq, pseq.cca,
                     type = "taxa", color = "Phylum")
p <- p + geom_point(size = 4)
print(p)



#Split plot
plot_ordination(pseq, pseq.cca,
  type = "split", shape = "nationality", 
  color = "Phylum", label = "nationality")



#t-SNE
library(vegan)
library(microbiome)
library(Rtsne) # Load package
set.seed(423542)

method <- "tsne"
trans <- "hellinger"
distance <- "euclidean"

# Distance matrix for samples
ps <- microbiome::transform(pseq, trans)

# Calculate sample similarities
dm <- vegdist(otu_table(ps), distance)

# Run TSNE
tsne_out <- Rtsne(dm, dims = 2) 
proj <- tsne_out$Y
rownames(proj) <- rownames(otu_table(ps))

library(ggplot2)
p <- plot_landscape(proj, legend = T, size = 1) 
print(p)



########## Regression ############
library(microbiome)
data(atlas1006)
plot_regression(diversity ~ age, meta(atlas1006))




######## DIFFERENTIAL ABUNDANCE TESTING #########

#Differential abundance testing: univariate data
#This section covers basic univariate tests for two-group comparison, covering t-test, Wilcoxon test, and multiple testing.

#The following example compares the abundance of a selected bug between two conditions. Let us assume that the data is already properly normalized.

library(microbiome)
data(dietswap)
d <- dietswap

# Pick microbial abundances for a given taxonomic group
taxa <- "Dialister"

# Construct a data.frame with the selected
# taxonomic group and grouping
df <- data.frame(Abundance = abundances(d)[taxa,],
                 Group = meta(d)$nationality)


#Compare the groups visually using a boxplot (left). However, we observe that the abundances are in absolute scale and therefore the comparison is not clear. Let us try the log10 transformation. Now, the data contains many zeros and taking log10 will yield infinite values. Hence we choose the commonly used, although somewhat problematic, log10(1+x) transformation (right).
p1 <- ggplot(df, aes(x = Group, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Absolute abundances", y = "Abundance\n (read count)")+ theme(plot.title = element_text(size=18))

# Let us add the log10(1+x) version:
df$Log10_Abundance <- log10(1 + df$Abundance)
p2 <- ggplot(df, aes(x = Group, y = Log10_Abundance)) +
  geom_boxplot() +
  labs(title = "Log10 abundances", y = "Abundance\n (log10(1+x) read count)")+ theme(plot.title = element_text(size=18))    

library(patchwork)
print(p1 + p2)


#The groups seem to differ. Let us test the difference statistically. First, let us perform t-test, which is based on Gaussian assumptions. Each group is expected to follow Gaussian distribution.

#Significance p-value with t-test:
print(t.test(Log10_Abundance ~ Group, data = df)$p.value)

#Now let us investigate the Gaussian assumption in more detail. Boxplots may not show deviations from Gaussian assumptions very clearly Let us try another visualization; the density plot.
p <- ggplot(df, aes(fill = Group, x = Log10_Abundance)) +
  geom_density(alpha = 0.5)
print(p)


#Apparently, the data is not Gaussian distributed. In such cases, a common procedure is to use non-parametric tests. These do not make assumptions of the data distribution but instead compare the ordering of the samples.

#So, let us look at the significance p-value with Wilcoxon test (log10 data):
print(wilcox.test(Log10_Abundance ~ Group, data = df)$p.value)

#But since the test is non-parametric, we can as well use the original absolute abundances since the log transformation does not change sample ordering on which the Wilcoxon test is based.

#Let us verify that the absolute abundances yield the same p-value for Wilcoxon test:
print(wilcox.test(Abundance ~ Group, data = df)$p.value)



#Let us compare how much the results would differ in the whole data between t-test and Wilcoxon test. To remove non-varying taxa that would demand extra scripting, let us for demonstration purposes now focus on core taxa that are observed in more than 20% of the samples with more than 3 reads.
# Core taxa to be tested
test.taxa <- core_members(d, prevalence = 20/100, detection = 3)

# Calculate p-values with the two different methods for each taxonomic unit
pvalue.ttest <- c()
pvalue.wilcoxon <- c()
for (taxa in test.taxa) {
  # Create a new data frame for each taxonomic group
  df <- data.frame(Abundance = abundances(d)[taxa,],
                   Log10_Abundance = log10(1 + abundances(d)[taxa,]),  
                   Group = meta(d)$nationality)
  
  pvalue.ttest[[taxa]] <- t.test(Log10_Abundance ~ Group, data = df)$p.value
  pvalue.wilcoxon[[taxa]] <- wilcox.test(Abundance ~ Group, data = df)$p.value  
}
# Arrange the results in a data.frame
pvalues <- data.frame(taxon = test.taxa,
                      pvalue.ttest = pvalue.ttest,
                      pvalue.wilcoxon = pvalue.wilcoxon)

# Note that multiple testing occurs.
# We must correct the p-values.
# let us apply the standard Benjamini-Hochberg False Discovery Rate (FDR)
# correction
pvalues$pvalue.ttest.adjusted <- p.adjust(pvalue.ttest)
#pvalues$pvalue.ttest.adjusted <- p.adjust(pvalues$pvalue.ttest)
pvalues$pvalue.wilcoxon.adjusted <- p.adjust(pvalue.wilcoxon)


#Compare the p-value histograms between raw and adjusteed p-values.
library(reshape2)
library(tidyverse)
pvalues$pvalue.wilcoxon<- as.numeric(pvalue.wilcoxon)
p1 <- ggplot(pvalues, aes(x = pvalue.wilcoxon)) + geom_histogram(bins = 50, binwidth = .03)+
  labs(title = "Raw p-values") +
  ylim(c(0, 80))

p2 <- ggplot(pvalues, aes(x = pvalue.wilcoxon.adjusted)) +
  geom_histogram(bins = 50, binwidth = .03) +
  labs(title = "Adjusted p-values") +
  ylim(c(0, 80))  

print(p1 + p2)


#Now compare these adjusted p-values between t-test and Wilcoxon test. Let us also highlight the p = 0.05 intervals.
p <- ggplot(data = pvalues,
            aes(x = pvalue.ttest.adjusted,
                y = pvalue.wilcoxon.adjusted)) +
  geom_text(aes(label = taxon)) + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_hline(aes(yintercept = 0.05), shape = 2) +
  geom_vline(aes(xintercept = 0.05), shape = 2)
print(p)


#Linear models: the role of covariates
#This section provides a brief hands-on introduction to the practical motivations and use linear (and generalized linear) models.

#Let us compare two groups with a linear model. We use Log10 abundances since this is closer to the Gaussian assumptions than the absolute count data. We can fit a linear model with Gaussian variation as follows:
res <- glm(Log10_Abundance ~ Group, data = df, family = "gaussian")

#Let us investigate model coefficients
kable(summary(res)$coefficients, digits = 5)

#The intercept equals to the mean in the first group:
print(mean(subset(df, Group == "AAM")$Log10_Abundance))

#The group term equals to the difference between group means:
print(mean(subset(df, Group == "AFR")$Log10_Abundance) -
        mean(subset(df, Group == "AAM")$Log10_Abundance))


#Note that the linear model (default) significance equals to t-test assuming equal variances.
print(t.test(Log10_Abundance ~ Group, data = df, var.equal=TRUE)$p.value)

#An important advantage of linear models, compared to plain t-test is that they allow incorporating additional variables, such as potential confounders (age, BMI, gender..):
# Add a covariate:
df$sex <- meta(d)$sex

# Fit the model:
res <- glm(Log10_Abundance ~ Group + sex, data = df, family = "gaussian")



#Generalized linear models: a very brief overview
#Let us briefly discuss the ideas underlying generalized linear models before proceeding to the next section.

#The Generalized linear model (GLM) allows a richer family of probability distributions to describe the data. Intuitively speaking, GLMs allow the modeling of nonlinear, nonsymmetric, and nongaussian associations. GLMs consist of three elements: - A probability distribution (from exponential family) - A linear predictor η = Xβ . - A link function g such that E(Y)=μ=g−1(η)

#We use Poisson with (its natural) log-link. Fit abundance (read counts) assuming that the data is Poisson distributed, and the logarithm of its mean, or expectation, is obtained with a linear model.
res <- glm(Abundance ~ Group, data = df, family = "poisson")
#Investigate the model output:
kable(summary(res)$coefficients, digits = 5)


#Advanced models of differential abundance
#GLMs are the basis for advanced testing of differential abundance in sequencing data. This is necessary, as the sequencing data sets deviate from symmetric, continuous, Gaussian assumptions in many ways.

#Sequencing data consists of discrete counts:
print(abundances(d)[1:5,1:3])

#the data is sparse:
hist(log10(1 + abundances(d)),100)

#Long tails of rare taxa:
medians <- apply(abundances(d),1,median)/1e3
library(reshape2)
A <- melt(microbiome::abundances(d))
A$Var1 <- factor(A$X1, levels = rev(names(sort(medians))))
p <- ggplot(A, aes(x = Var1, y = value)) +
  geom_boxplot() +
  labs(y = "Abundance (reads)", x = "Taxonomic Group") +
  scale_y_log10()

print(p)


#Overdispersion (variance exceeds the mean):
means <- apply(abundances(d),1,mean)
variances <- apply(abundances(d),1,var)

# Calculate mean and variance over samples for each taxon
library(reshape2)
library(dplyr)
df <- melt(abundances(d))
names(df) <- c("Taxon", "Sample", "Reads")
df <- df %>% group_by(Taxon) %>%
  summarise(mean = mean(Reads),
            variance = var(Reads))

# Illustrate overdispersion
library(scales)
p <- ggplot(df, aes(x = mean, y = variance)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  labs(title = "Overdispersion (variance > mean)")
print(p)



#DESeq2: differential abundance testing for sequencing data
#DESeq2 analysis can accommodate those particular assumptions about sequencing data.
# Start by converting phyloseq object to deseq2 format
library(DESeq2)
ds2 <- phyloseq_to_deseq2(d, ~ group + nationality)

# Run DESeq2 analysis (all taxa at once!)
dds <- DESeq(ds2)

# Investigate results
res <- results(dds)
deseq.results <- as.data.frame(res)
df <- deseq.results
df$taxon <- rownames(df)
df <- df %>% arrange(log2FoldChange, padj)

# Print the results; flitered and sorted by pvalue and effectsize
library(knitr)
df <- df %>% filter(pvalue < 0.05 & log2FoldChange > 1.5) %>%
  arrange(pvalue, log2FoldChange)
kable(df, digits = 5)



#For comparison purposes, assess significances and effect sizes based on Wilcoxon test.
test.taxa <- taxa(d)
pvalue.wilcoxon <- c()
foldchange <- c()
for (taxa in test.taxa) {
  # Create a new data frame for each taxonomic group
  df <- data.frame(Abundance = abundances(d)[taxa,],
                   Log10_Abundance = log10(1 + abundances(d)[taxa,]),
                   Group = meta(d)$nationality)
  # Calculate pvalue and effect size (difference beween log means)       
  pvalue.wilcoxon[[taxa]] <- wilcox.test(Abundance ~ Group, data = df)$p.value
  foldchange[[taxa]] <- coef(lm(Log10_Abundance ~ Group, data = df))[[2]]
}
# Correct p-values for multiple testing
pvalue.wilcoxon.adjusted <- p.adjust(pvalue.wilcoxon)


par(mfrow = c(1,2))
plot(deseq.results$padj, pvalue.wilcoxon.adjusted,
     xlab = "DESeq2 adjusted p-value",
     ylab = "Wilcoxon adjusted p-value",
     main = "P-value comparison")
abline(v = 0.05, h = 0.05, lty = 2)

plot(deseq.results$log2FoldChange, foldchange, 
     xlab = "DESeq2",
     ylab = "Linear model",
     main = "Effect size comparison")
abline(0,1)
















