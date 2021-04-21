#### An analysis of cyanobacteria 
#### in Lake Champlain

setwd("/Users/adw96/Dropbox/LakeChamplain/")
require("biom")
require(devtools)
devtools::install_github("adw96/breakaway")
library(breakaway)
require(ggplot2)
require(reshape2)
require(dplyr)

###########################################################
################## READ IN DATA ###########################
###########################################################

############# OTUs ########################################
otus_500 <-  read.table(file = "MED/filtered_MATRIX-COUNT_M500.txt", header = T)

meds <-  as.character(otus_500[,1])
set.seed(1)
otus <- "otus_500" 
### create otu and taxonomy
my_table <- get(otus)
my_taxonomy <- my_table[, (names(my_table) %in% c("taxonomy", "taxa"))]
my_table <- my_table[, !(names(my_table) %in% c("Node_ID", "Nodes", "taxa", "taxonomy"))]
assign(paste(otus, "_otus", sep=""), as.matrix(my_table))
assign(paste(otus, "_taxonomy", sep=""), as.character(my_taxonomy))
rownames(my_table) <- meds


## look at cyanobacteria
my_total <- apply(as.matrix(my_table), 2, sum)
assign(paste(otus, "_count", sep=""), my_total)

## consider the rows with taxonomy 
my_microcystis <- apply(my_table[grepl("Microcystis", my_taxonomy), ], 2, sum)
my_dolicho <- apply(my_table[grepl("Dolichospermum", my_taxonomy), ], 2, sum)
assign(paste(otus, "_microcystis", sep=""), my_microcystis)
assign(paste(otus, "_dolichospermum", sep=""), my_dolicho)

## diversity indices
shannon_resamples <- t(replicate(20, apply(get(paste(otus, "_otus", sep="")), 2,  breakaway::resample_estimate, breakaway::shannon, my_total)))
shannon_e_resamples <- t(replicate(20, apply(get(paste(otus, "_otus", sep="")), 2,  breakaway::resample_estimate, breakaway::shannon_e, my_total)))
simpson_resamples <- t(replicate(20, apply(get(paste(otus, "_otus", sep="")), 2,  breakaway::resample_estimate, breakaway::simpson, my_total)))
assign(paste(otus, "_shannon", sep=""), shannon_resamples)
assign(paste(otus, "_shannon_e", sep=""), shannon_e_resamples)
assign(paste(otus, "_simpson", sep=""), simpson_resamples)

## create frequency tables
frequencytablelist <- lapply(apply(as.matrix(my_table),2,table),as.data.frame)
frequencytablelist <- lapply(frequencytablelist,function(x) x[x[,1]!=0,])

assign(paste(otus, "_frequency_tables", sep=""), frequencytablelist)

### create otu and taxonomy
my_breakaway <- matrix(NA, nrow=length(otus_500_frequency_tables), ncol=2)
my_objective <- matrix(NA, nrow=length(otus_500_frequency_tables), ncol=2)
for (i in 1:length(otus_500_frequency_tables)) {
  ## originally I tried to do this with lapply but occasionally  breaks
  temp <- try(breakaway::breakaway(otus_500_frequency_tables[[i]], print = F, answers = T, plot = F), silent = TRUE)
  #my_breakaway[i, ] <- ifelse(rep(class(temp) != "NULL", 2), c(NA, NA), c(temp$est, temp$seest))
  try(my_breakaway[i, ] <- c(temp$est, temp$seest), silent=T)
}
# FIX!
for (i in 1){#:length(otus_500_frequency_tables)) {
  temp2 <- try(breakaway::objective_bayes_negbin(data = otus_500_frequency_tables[[i]], plot = F, print = F, answers = T)$results, silent = TRUE)
  #my_objective[i, ] <- ifelse(class(temp2) == rep("try-error", 2), c(NA, NA), c(temp2["mean.N", ], temp2["stddev.C", ]))
  try(my_objective[i, ] <- c(temp2["mean.N", ], temp2["stddev.C", ]), silent=T)
  
}
objective_bayes_negbin(data = otus_500_frequency_tables[[i]], plot = T, print = F, answers = T)

assign(paste(otus, "_breakaway", sep=""), my_breakaway)
assign(paste(otus, "_objective", sep=""), my_objective)


### within cyanobacteria diversity
my_microcystis_counts <- my_table[grepl("Microcystis", my_taxonomy), ]
my_dolicho_counts <- my_table[grepl("Dolichospermum", my_taxonomy), ]
assign(paste(otus, "_microcystis_counts", sep=""), my_microcystis_counts)
assign(paste(otus, "_dolichospermum_counts", sep=""), my_dolicho_counts)

## diversity indices
shannon_e_resamples_microcystis <- t(replicate(20, apply(my_microcystis_counts, 2, resample_estimate, shannon_e, my_total)))
simpson_resamples_microcystis <- t(replicate(20, apply(my_microcystis_counts, 2, resample_estimate, simpson, my_total)))
assign(paste(otus, "_shannon_e_microcystis", sep=""), shannon_e_resamples_microcystis)
assign(paste(otus, "_simpson_microcystis", sep=""), simpson_resamples_microcystis)

shannon_e_resamples_dolicho <- t(replicate(20, apply(my_dolicho_counts, 2, resample_estimate, shannon_e, my_total)))
simpson_resamples_dolicho <- t(replicate(20, apply(my_dolicho_counts, 2, resample_estimate, simpson, my_total)))
assign(paste(otus, "_shannon_e_dolicho", sep=""), shannon_e_resamples_dolicho)
assign(paste(otus, "_simpson_dolicho", sep=""), simpson_resamples_dolicho)

###########################################################
################## READ IN METADATA #######################
###########################################################

meta <- read.table(file = "MED/mapping_file_environm.txt", header = T)
rownames(meta) <- paste("X", meta[,1], sep="")
usable <- intersect(rownames(meta), colnames(shannon_e_resamples_dolicho))

Year <- meta$Years
names(Year) <- rownames(meta)
Year <- Year[colnames(my_microcystis_counts)]
Month <- meta$Months; names(Month) <- rownames(meta)
Month <- Month[colnames(my_microcystis_counts)]

Date <- as.Date(colnames(my_microcystis_counts), "X%d.%m.%Y")
tmp <- as.POSIXlt(Date[(Date < "1900-01-01") ])
tmp$year <- tmp$year+2000
Date[(Date < "1900-01-01") ] <- format(tmp)


#save.image("/Users/adw96/Dropbox/LakeChamplain/AmyRScripts/workspace_161222.RData")

############################################################
############ FANCY BARPLOTS ################################
############################################################

###### ###### ###### ###### ###### ###### ###### ###### 
###### Stacked bar plot with MED proportions
###### ###### ###### ###### ###### ###### ###### ###### 


#rownames(my_microcystis_counts) <- paste("M", 1:dim(my_microcystis_counts)[1], sep="")
#rownames(my_dolichospermum_counts) <- paste("D", 1:dim(my_dolichospermum_counts)[1], sep="")
#my_cyanobacteria_matrix <- (data.frame( rbind(my_microcystis_counts, my_dolichospermum_counts), Type = c(rep("M", dim(my_microcystis_counts)[1]), rep("D", dim(my_dolichospermum_counts)[1]))))
# my_cyanobacteria_matrix <- data.frame(Year = Year,
#                                       Month = Month,
#                                       Sample=colnames(my_microcystis_counts),
#                                       cbind(t(my_microcystis_counts), t(my_dolichospermum_counts)))#, Type = c(rep("M", dim(my_microcystis_counts)[1]), rep("D", dim(my_dolichospermum_counts)[1]))))
my_microcystis_counts_named <- my_microcystis_counts
rownames(my_microcystis_counts_named) <- paste("M", rownames(my_microcystis_counts_named), sep="")
my_dolicho_counts_named <- my_dolicho_counts
rownames(my_dolicho_counts_named) <- paste("D", rownames(my_dolicho_counts_named), sep="")
my_cyanobacteria_matrix <- 
  data.frame(Sample= colnames(my_microcystis_counts), 
             cbind(t(my_microcystis_counts_named), t(my_dolicho_counts_named)))#, Type = c(rep("M", dim(my_microcystis_counts)[1]), rep("D", dim(my_dolichospermum_counts)[1]))))


my_cyanobacteria_matrix <- my_cyanobacteria_matrix[order(Date), ]

colnames(my_microcystis_counts)
Year = Year[colnames(my_microcystis_counts)]
colnames(my_microcystis_counts) ## 150
rownames(meta) ## 143
meta$Years[colnames(my_microcystis_counts)]
colnames(my_microcystis_counts)

###### ###### ###### ###### 
###### everyone
###### ###### ###### ###### 
my_width <- 20
#### set up basic df and plot
my_cyanobacteria_matrix_renamed <- my_cyanobacteria_matrix
cbind(rownames(my_cyanobacteria_matrix), as.character(Date[order(Date)])) ## yes, matching
df <- melt(my_cyanobacteria_matrix_renamed, id = "Sample")

## create date column
get_date <- function(rowname) {
  as.Date(ifelse(grepl(".20", rowname), as.Date(rowname, "X%d.%m.%Y"),  as.Date(rowname, "X%d.%m.%y")), origin="1970-01-01")
}
data.frame("changed"=get_date(df$Sample), 
           "original"=df$Sample)
df$Date <- get_date(df$Sample)

df <- df[order(df$Date, decreasing = F), ]

g <- ggplot(df, aes(x=Sample, y = value, fill = variable)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_x_discrete()
g + geom_bar(stat='identity', position = "fill")

#### incorporate bloom info
## start by creating a dates and bloom data frame
bloom_df <- data.frame("ID"=meta$SampleID, "bloom"=meta$bloom2, "XID"=paste("X", meta$SampleID, sep=""))
intersect(bloom_df$XID, my_cyanobacteria_matrix_renamed$Sample)
strsplit(intersect(bloom_df$XID, my_cyanobacteria_matrix_renamed$Sample), split = "X")

bloom_df_useable <- data.frame("Sample"=my_cyanobacteria_matrix_renamed$Sample, "bloom"=NA)
for (i in 1:length(my_cyanobacteria_matrix_renamed$Sample)) {
  j <- my_cyanobacteria_matrix_renamed$Sample[i]
  if (j %in% intersect(bloom_df$XID, my_cyanobacteria_matrix_renamed$Sample)) {
    bloom_df_useable$bloom[i] <- meta$bloom2[which(as.character(j) == bloom_df$XID)] == "yes"
  }
}
bloom_df_useable
bloom_df_useable$Sample <- as.character(bloom_df_useable$Sample)
rownames(bloom_df_useable) <- bloom_df_useable$Sample

###### All nodes with bloom info
df$bloom <- bloom_df_useable[df$Sample, "bloom"]
g_all <- ggplot(df, aes(x=Date, y = value, fill = variable, bloom)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_alpha_discrete(range=c(0.5, 0.9)) +
  theme_bw()

pdf("AmyPlots/all-barplot.pdf", width = 20, height=5)
g_all + geom_bar(stat='identity', position = "fill", aes(alpha = bloom), width=my_width) +
  scale_fill_manual(values = c(
    '#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d', # reds
    '#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d', # purples
    '#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b', # blues
    '#c7e9c0','#74c476','#41ab5d','#238b45','#006d2c','#00441b', # greens
    '#bdbdbd','#969696','#737373','#525252','#252525','#000000'))   # greys
dev.off()

###### M nodes
df_m <- na.omit(df[substr(df$variable,1,1) == "M", ])
g_m <- ggplot(df_m, 
              aes(x=Date, y = value, fill = variable, bloom)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_alpha_discrete(range=c(0.6, 0.9)) +
  theme_bw()


pdf("AmyPlots/m-barplot.pdf", width = 20, height=5)
g_m + geom_bar(stat='identity', position = "fill", aes(alpha = bloom), width=my_width) + 
  scale_fill_discrete(h=c(0,80))
dev.off()

###### D nodes
df_d <- na.omit(df[substr(df$variable,1,1) == "D", ])
g_d <- ggplot(df_d, aes(x=Date, y = value, fill = variable, bloom)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_alpha_discrete(range=c(0.6, 0.9)) + 
  theme_bw()

pdf("AmyPlots/d-barplot.pdf", width = 20, height=5)
g_d + geom_bar(stat='identity', position = "fill", aes(alpha = bloom), width = my_width) + scale_fill_discrete(h=c(220,280))
dev.off()

###### M vs D nodes
df_md_too_big <- data.frame(df, "cyano"=substr(df$variable,1,1))
df_md <- na.omit(df_md_too_big %>% group_by(Sample, cyano, bloom, Date) %>% summarise(total = sum(value)))
g_md <- ggplot(df_md, aes(x=Date, y = total, fill = cyano, bloom)) +
  scale_alpha_discrete(range=c(0.6, 0.9))
pdf("AmyPlots/md-barplot.pdf", width = 20, height=5)
g_md + geom_bar(stat='identity', position = "fill", aes(alpha = bloom), width = my_width) + scale_fill_discrete(h=c(250, 0))
dev.off()

df_md_collapse <- na.omit(df_md_too_big %>% group_by(Date, cyano, bloom) %>% summarise(total = sum(value)))
g_md_collapse <- ggplot(df_md_collapse, aes(x=Date, y = total, fill = cyano, bloom)) +
  scale_alpha_discrete(range=c(0.6, 0.9))
pdf("AmyPlots/md-barplot-collapse-same-day.pdf", width = 20, height=5)
g_md_collapse + geom_bar(stat='identity', position = "fill", aes(alpha = bloom), width = my_width) + scale_fill_discrete(h=c(250, 0))
dev.off()

####### Irrelevant: Nico wants whole community
# ## new coloured line plot
# df_md
# df_wide <- data.frame(acast(df_md, Sample~cyano, value.var="total"))
# names(df_wide) <- c("Sample", "D", "M")
# ggplot(df_wide, aes(x=Sample, ymin = D, ymax = M)) +
#   geom_line() +
#   geom_ribbon()
# reshape(df_md, timevar = c("cyano"), idvar = c("Sample"), direction = "wide")
# reshape(df_md, timevar = c("cyano"), idvar = c("bloom", "Sample", "Date"), direction = "wide")
# df_md_collapse$cyano == "D"

#### Ribbon plot for Nico
df_ribbon <- data.frame("Date" = get_date(names(my_total)), 
                        "M" = my_microcystis/my_total, 
                        "D" = my_dolicho/my_total)
rownames(df_ribbon) <- names(my_total)
df_ribbon <- df_ribbon[order(get_date(rownames(df_ribbon)), decreasing = F), ]
df_ribbon <- cbind(df_ribbon, "min_line"=pmin(df_ribbon$M, df_ribbon$D))
df_ribbon <- melt(df_ribbon, id.vars=c("Date","min_line"), variable.name="Cyano", value.name="Proportions")
names(df_ribbon)
head(df_ribbon)
head(df)

pdf("AmyPlots/md-ribbon.pdf")
ggplot(data=df_ribbon, aes(x=Date, fill=Cyano)) + 
  geom_ribbon(aes(ymax=Proportions, ymin=min_line)) +
  scale_fill_manual(values=c(M="darkred", D="darkblue")) +
  labs(y = "Relative abundance") +
  theme_bw()
dev.off()
###
#save.image("/Users/adw96/Dropbox/LakeChamplain/AmyRScripts/workspace_170108.RData")

#### End Jan 8 2016 ########
##############################


############################################################
############################################################
############################################################

## what type of singletons filtering is present?
head(sort(apply(otus_500_otus, 1, sum)), 50)

##### make some pretty pictures

pdf("AmyPlots/Shannon.pdf")
par(mfrow=c(3,3))
### microcystis, shannon
for (otus in c("otus_500")) {
  ### create otu and taxonomy
  my_shannon <- get(paste(otus, "_shannon", sep=""))
  my_simpson <- get(paste(otus, "_simpson", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  breakaway::betta_pic(apply(my_shannon, 2, mean), apply(my_shannon, 2, sd), my_microcystis, mymain = otus, myy = "Shannon index")
  axis(1, at=0.35, tick = F, labels = "Microcystis proportion", pos = -1)
}
## d, shannon
for (otus in c("otus_500")) {
  ### create otu and taxonomy
  my_shannon <- get(paste(otus, "_shannon", sep=""))
  my_simpson <- get(paste(otus, "_simpson", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  breakaway::betta_pic(apply(my_shannon, 2, mean), apply(my_shannon, 2, sd), my_dolicho, mymain = otus, myy = "Shannon index")
  axis(1, at=0.35, tick = F, labels = "Dolichospermum proportion", pos = -1)
}
for (otus in c("otus_500")) {
  ### create otu and taxonomy
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  plot(my_microcystis, my_dolicho, main=otus)
}
dev.off()

pdf("AmyPlots/Simpson.pdf")
par(mfrow=c(2,3))
par(xpd=FALSE)
# micro,  simpson
for (otus in c("otus_500")) {
  ### create otu and taxonomy
  my_simpson <- get(paste(otus, "_simpson", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  plot(my_microcystis, apply(my_simpson, 2, mean), type="n", main = otus, ylab = "Simpson index", xlab = "Microcystis proportion")
  betta_points(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_microcystis)
  #breakaway::betta_pic(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_microcystis,)
}
## d, simpson
for (otus in c("otus_500")) {
  ### create otu and taxonomy
  my_simpson <- get(paste(otus, "_simpson", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  plot(my_dolicho, apply(my_simpson, 2, mean), type="n", main = otus, ylab = "Simpson index", xlab = "Dolichospermum proportion")
  betta_points(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_dolicho)
  #breakaway::betta_pic(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_dolicho, mymain = otus, myy = "Simpson index")
}
dev.off()

pdf("AmyPlots/Simpson_e.pdf")
par(mfrow=c(2,3))
par(xpd=FALSE)
# micro,  simpson
for (otus in c("otus_500")) {
  ### create otu and taxonomy
  my_simpson <- get(paste(otus, "_shannon_e", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  plot(my_microcystis, apply(my_simpson, 2, mean), type="n", main = otus, ylab = "Shannon/log(sample richness)", xlab = "Microcystis proportion")
  betta_points(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_microcystis)
  #breakaway::betta_pic(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_microcystis,)
}
## d, simpson
for (otus in c("otus_500")) {
  ### create otu and taxonomy
  my_simpson <- get(paste(otus, "_shannon_e", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  plot(my_dolicho, apply(my_simpson, 2, mean), type="n", main = otus, ylab = "Shannon/log(sample richness)", xlab = "Dolichospermum proportion")
  betta_points(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_dolicho)
  #breakaway::betta_pic(apply(my_simpson, 2, mean), apply(my_simpson, 2, sd), my_dolicho, mymain = otus, myy = "Simpson index")
}
dev.off()

pdf("AmyPlots/Richness_breakaway.pdf")
par(mfrow=c(2,3))
par(xpd=F)
# micro,  richness
for (otus in c("otus_500")) {
  my_breakaway <- get(paste(otus, "_breakaway", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(my_breakaway[,1], my_breakaway[,2], my_microcystis, ylimu = 3000, myy = "Estimated species richness (breakaway)", mymain = "Species richness with M")
}
## d, breakaway
for (otus in c("otus_500")) {
  my_breakaway <- get(paste(otus, "_breakaway", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(my_breakaway[,1], my_breakaway[,2], my_dolicho, ylimu = 3000, myy = "Estimated species richness (breakaway)", mymain = "Species richness with D")
}
dev.off()

pdf("AmyPlots/Richness_objective_bayes.pdf")
par(mfrow=c(2,3))
par(xpd=F)
# micro,  richness
for (otus in c("otus_500")) {
  my_objective <- get(paste(otus, "_objective", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(my_objective[,1], my_objective[,2], my_microcystis, ylimu = 1000, myy = "Estimated species richness (objective Bayes)", mymain = "Species richness with M")
}
## d, objective
for (otus in c("otus_500")) {
  my_objective <- get(paste(otus, "_objective", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(my_objective[,1], my_objective[,2], my_dolicho, ylimu = 1000, myy = "Estimated species richness (objective Bayes)", mymain = "Species richness with D")
}
dev.off()

###

corrected_shannon_breakaway <- function(my_data) {
  x <- as.data.frame((table(my_data)))
  x <- x[x[,1]!=0,]
  y <- try(shannon(x)/log(breakaway::breakaway(x, plot = F, print = F, answers = T)$est), silent = T)
  ifelse(class(y) == "numeric", y, NA)
}

corrected_shannon_objective <- function(my_data) {
  x <- as.data.frame((table(my_data)))
  x <- x[x[,1]!=0,]
  y <- try(shannon(x)/log(breakaway::objective_bayes_negbin(x, plot = F, print = F, answers = T)$results["mean.N", ]), silent = T)
}

##### run overnight
set.seed(2)
for (otus in c("otus_500")) {
  number_iterations <- 20
  corrected_shannon_breakaway_resamples <- t(replicate(number_iterations, unlist(apply(get(paste(otus, "_otus", sep="")), 2, resample_estimate, corrected_shannon_breakaway, my_total))))
  corrected_shannon_objective_resamples <- t(replicate(number_iterations, unlist(apply(get(paste(otus, "_otus", sep="")), 2, resample_estimate, corrected_shannon_objective, my_total))))
  corrected_shannon_objective_resamples1 <- matrix(as.numeric(corrected_shannon_objective_resamples), nrow=number_iterations)
  assign(paste(otus, "_corrected_shannon_breakaway", sep=""), corrected_shannon_breakaway_resamples )
  assign(paste(otus, "_corrected_shannon_objective", sep=""), corrected_shannon_objective_resamples1 )
}

########### evenness adjusted for richness
pdf("AmyPlots/Shannon_e_breakaway.pdf")
par(mfrow=c(2,3))
par(xpd=F)
# micro,  shannon e
for (otus in c("otus_500")) {
  my_objective <- get(paste(otus, "_corrected_shannon_breakaway", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(apply(my_objective, 2, mean, na.rm=T), apply(my_objective, 2, sd, na.rm=T), my_microcystis, myy = "Estimated Shannon E (breakaway)", mymain = "Shannon E with M")
}
## d,  shannon e
for (otus in c("otus_500")) {
  my_objective <- get(paste(otus, "_corrected_shannon_breakaway", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(apply(my_objective, 2, mean, na.rm=T), apply(my_objective, 2, sd, na.rm=T), my_dolicho, myy = "Estimated Shannon E (breakaway)", mymain = "Shannon E with D")
}
dev.off()
pdf("AmyPlots/Shannon_e_objective.pdf")
par(mfrow=c(2,3))
par(xpd=F)
# micro,  shannon e
for (otus in c("otus_500")) {
  my_objective <- get(paste(otus, "_corrected_shannon_objective", sep=""))
  my_microcystis <- get(paste(otus, "_microcystis", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(apply(my_objective, 2, mean, na.rm=T), apply(my_objective, 2, sd, na.rm=T), my_microcystis, myy = "Estimated Shannon E (breakaway)", mymain = "Shannon E with M")
}
## d,  shannon e
for (otus in c("otus_500")) {
  my_objective <- get(paste(otus, "_corrected_shannon_objective", sep=""))
  my_dolicho <- get(paste(otus, "_dolichospermum", sep=""))/get(paste(otus, "_count", sep=""))
  betta_pic(apply(my_objective, 2, mean, na.rm=T), apply(my_objective, 2, sd, na.rm=T), my_dolicho, myy = "Estimated Shannon E (objective Bayes)", mymain = "Shannon E with D")
}
dev.off()

### within cyanobacteria diversity
for (otus in c("otus_500")) {
  my_taxonomy <- get(paste(otus, "_taxonomy", sep=""))
  my_table <- get(paste(otus, "_otus", sep=""))
  
  ## consider the rows with taxonomy 
  my_microcystis_counts <- my_table[grepl("Microcystis", my_taxonomy), ]
  my_dolicho_counts <- my_table[grepl("Dolichospermum", my_taxonomy), ]
  assign(paste(otus, "_microcystis_counts", sep=""), my_microcystis_counts)
  assign(paste(otus, "_dolichospermum_counts", sep=""), my_dolicho_counts)
  
  ## diversity indices
  shannon_e_resamples_microcystis <- t(replicate(20, apply(my_microcystis_counts, 2, resample_estimate, breakaway::shannon_e, my_total)))
  simpson_resamples_microcystis <- t(replicate(20, apply(my_microcystis_counts, 2, resample_estimate, breakaway::simpson, my_total)))
  assign(paste(otus, "_shannon_e_microcystis", sep=""), shannon_e_resamples_microcystis)
  assign(paste(otus, "_simpson_microcystis", sep=""), simpson_resamples_microcystis)
  
  shannon_e_resamples_dolicho <- t(replicate(20, apply(my_dolicho_counts, 2, resample_estimate, breakaway::shannon_e, my_total)))
  simpson_resamples_dolicho <- t(replicate(20, apply(my_dolicho_counts, 2, resample_estimate, breakaway::simpson, my_total)))
  assign(paste(otus, "_shannon_e_dolicho", sep=""), shannon_e_resamples_dolicho)
  assign(paste(otus, "_simpson_dolicho", sep=""), simpson_resamples_dolicho)
}


############# Metadata ####################################

meta <- read.table(file = "MED/mapping_file_environm.txt", header = T)
rownames(meta) <- paste("X", meta[,1], sep="")
usable <- intersect(rownames(meta), colnames(shannon_e_resamples_dolicho))

new_uncertainty_figure <- function(resamples, covariate, ylab ="", main = "") {
  means <- apply(resamples, 2, mean)
  sds <- apply(resamples, 2, sd)
  means_subset <- means[usable] 
  sds_subset <- sds[usable]
  names(covariate) <- rownames(meta)
  covariate_subset <- covariate[usable]
  betta_pic( means_subset[sds_subset>0], sds_subset[sds_subset>0], covariate_subset[sds_subset>0], myy = ylab, mymain = main)
}

### want to plot evenness with m and d as a function of temperature, TN, TP, DN, DP
pdf("AmyPlots/Within-cyano-diversity-nutrients.pdf")
par(mfrow=c(2,3))
variables <- c("TP", "TN", "DP", "DN")
for (covariate1 in 1:length(variables)) {
  covariate  <- meta[,variables[covariate1]]
  for (otus in c("otus_500")) {
    shannon_e_resamples_microcystis <- get(paste(otus, "_shannon_e_microcystis", sep=""))
    new_uncertainty_figure(shannon_e_resamples_microcystis, log(covariate), ylab = "Microcystis Shannon e", main = paste("Shannon e of M with log of", variables[covariate1]))
  }
  for (otus in c("otus_500")) {
    shannon_e_resamples_dolicho <- get(paste(otus, "_shannon_e_dolicho", sep=""))
    new_uncertainty_figure(shannon_e_resamples_dolicho, log(covariate), ylab = "Dolicho Shannon e", main = paste("Shannon e of D with log of", variables[covariate1]))
  }
  for (otus in c("otus_500")) {
    simpson_resamples_microcystis <- get(paste(otus, "_simpson_microcystis", sep=""))
    new_uncertainty_figure(simpson_resamples_microcystis, log(covariate), ylab = "Microcystis Simpson", main = paste("Simpson of M with log of", variables[covariate1]))
  }
  for (otus in c("otus_500")) {
    simpson_resamples_dolicho <- get(paste(otus, "_simpson_dolicho", sep=""))
    new_uncertainty_figure(simpson_resamples_dolicho, log(covariate), ylab = "Dolicho Simpson", main = paste("Simpson of D with log of", variables[covariate1]))
  }
}
dev.off()
pdf("AmyPlots/Within-cyano-diversity-temp.pdf")
par(mfrow=c(2,3))
variables <- c("Temperature_Water_Celsius")
covariate  <- meta[,variables[covariate1]]
for (otus in c("otus_500")) {
  shannon_e_resamples_microcystis <- get(paste(otus, "_shannon_e_microcystis", sep=""))
  new_uncertainty_figure(shannon_e_resamples_microcystis, covariate, ylab = "Microcystis Shannon e", main = paste("Shannon e of M with Water Temp"))
}
for (otus in c("otus_500")) {
  shannon_e_resamples_dolicho <- get(paste(otus, "_shannon_e_dolicho", sep=""))
  new_uncertainty_figure(shannon_e_resamples_dolicho, covariate, ylab = "Dolicho Shannon e", main = paste("Shannon e of D with Water Temp"))
}
for (otus in c("otus_500")) {
  simpson_resamples_microcystis <- get(paste(otus, "_simpson_microcystis", sep=""))
  new_uncertainty_figure(simpson_resamples_microcystis, covariate, ylab = "Microcystis Simpson", main = paste("Simpson of M with Water Temp"))
}
for (otus in c("otus_500")) {
  simpson_resamples_dolicho <- get(paste(otus, "_simpson_dolicho", sep=""))
  new_uncertainty_figure(simpson_resamples_dolicho, covariate, ylab = "Dolicho Simpson", main = paste("Simpson of D with Water Temp"))
}
dev.off()

#save.image("/Users/adw96/Dropbox/LakeChamplain/AmyRScripts/workspace_161013.RData")

#####
#load("/Users/adw96/Dropbox/LakeChamplain/AmyRScripts/workspace_161013.RData")
#####
#### look at cyanobacteria richness
rownames(otus_500_microcystis_counts)
my_taxonomy[grepl("Microcystis", my_taxonomy)]

for (otus in c("otus_500")) {
  my_taxonomy <- get(paste(otus, "_taxonomy", sep=""))
  my_table <- get(paste(otus, "_otus", sep=""))
  my_microcystis_counts <- my_table[grepl("Microcystis", my_taxonomy), ]
  rownames(my_microcystis_counts) <- my_taxonomy[grepl("Microcystis", my_taxonomy)]
  my_dolicho_counts <- my_table[grepl("Dolichospermum", my_taxonomy), ]
  rownames(my_dolicho_counts) <- my_taxonomy[grepl("Dolichospermum", my_taxonomy)]
  assign(paste(otus, "_microcystis_counts", sep=""), my_microcystis_counts)
  assign(paste(otus, "_dolichospermum_counts", sep=""), my_dolicho_counts)
  
  my_microcystis_richness <- apply(my_microcystis_counts, 2, function(x) sum(x>0))
  my_dolicho_richness <- apply(my_dolicho_counts, 2, function(x) sum(x>0))
  assign(paste(otus, "_microcystis_richness", sep=""), my_microcystis_richness)
  assign(paste(otus, "_dolichospermum_richness", sep=""), my_dolicho_richness)
}

pdf("AmyPlots/within-cyano-richness.pdf")
par(mfrow=c(2,3))
variables <- c("TP", "TN", "DP", "DN", "Temperature_Water_Celsius")
names(my_microcystis_richness)
for (covariate1 in 1:length(variables)) {
  covariate  <- meta[,variables[covariate1]]
  names(covariate) <- rownames(meta)
  this_xlab <- variables[covariate1]
  if (covariate1 < 5) { covariate <- log(covariate); this_xlab <- paste("Log of", this_xlab) }
  for (otus in c("otus_500")) {
    my_microcystis_richness <- get(paste(otus, "_microcystis_richness", sep=""))
    plot(covariate[usable], my_microcystis_richness[usable], xlab=this_xlab, ylab="Number of different Microcystis MEDs present", main=otus)
  }
  for (otus in c("otus_500")) {
    my_dolicho_richness <- get(paste(otus, "_dolichospermum_richness", sep=""))
    plot(covariate[usable], my_dolicho_richness[usable], xlab=this_xlab, ylab="Number of different Dolicho MEDs present", main=otus)
  }
}
dev.off()

#####
#load("/Users/adw96/Dropbox/LakeChamplain/AmyRScripts/workspace_161020.RData")
#####

####### OLD STUFF
############# Oligotype ###################################
otus_500_taxonomy
otus_500_taxonomy[grepl("Cyano", otus_500_taxonomy)]
otus_500_taxonomy[grepl("Microcystis", otus_500_taxonomy)]
otus_500_taxonomy[grepl("Dolichospermum", otus_500_taxonomy)]

oligotypes <- read.table("MATRIX-COUNT_oligotypes_Microcystis.txt", header = T)
colnames(oligotypes) <- substring(colnames(oligotypes), 2, 100000)
sort(colnames(oligotypes))
sort(names(richness_estimates))
intersect(sort(colnames(oligotypes)), sort(names(richness_estimates)))
oligotypes <- oligotypes[, intersect(sort(colnames(oligotypes)), sort(names(richness_estimates)))]

############# Metadata ####################################

meta <- read.table("mapping_bloom2.txt", header = T)
rownames(meta) <- paste("X", rownames(meta), sep="")

### check it matches otus
usable <- intersect(rownames(covariates), names(otus))
covariates <- meta[usable, ]
otus <- otus[, usable]


###########################################################
# ANALYSIS
###########################################################

## investigate how total richness is affected by the oligotypes

## first construct frequency tables

frequency <- lapply(apply(otus,2,table),as.data.frame) 
frequency <- lapply(frequency,function(x) x[x[,1]!=0,])
head(frequency)

## definitely breakaway structures
blank_template <- rep(NA, length(frequency)); names(blank_template) <- names(frequency)
richness_estimates <- blank_template
richness_errors <- blank_template
richness_models <- blank_template
for (i in 132:length(frequency)) {
  breakaway_results <- breakaway(frequency[[i]], print = FALSE, answers = TRUE, plot = FALSE)
  if (!is.null(breakaway_results)) {
    richness_estimates[i] <- breakaway_results$est
    richness_errors[i] <- breakaway_results$seest
    richness_models[i] <- breakaway_results$name    
  }
}
cbind(richness_estimates, richness_errors)
table(richness_models)



betta_pic(richness_estimates, richness_errors, ylimu = 20000)

richness_estimates[richness_errors < 5] <- NA

betta_pic(richness_estimates, richness_errors, x = covariates[, "Temperature_Water_Celsius"], ylimu = 20000)
betta_pic(richness_estimates, richness_errors, x = covariates[, "Mean_temperature_t0_t7"], ylimu = 20000)

global <- betta(richness_estimates, richness_errors)$table[1,1]

#pdf("Figures/first_oligotypes_look.pdf")
par(mfrow=c(1,2))
for (i in 1:4) {
  betta_pic(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], x = oligotypes[i, ], ylimu = 20000, mymain = rownames(oligotypes)[i], myy = "Est'd species richness"); abline(h=global)
  betta_pic(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], x = oligotypes[i, ]/apply(oligotypes, 2, sum), ylimu = 20000, mymain = rownames(oligotypes)[i], myy = "Est'd species richness"); abline(h=global)
}
par(mfrow=c(1,1))
betta_pic(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], x = apply(oligotypes, 2, sum), ylimu = 20000, mymain = "Total number of observed Microcystis", myy = "Est'd species richness"); abline(h=global)
dev.off()
betta_pic(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], x = oligotypes[1, ], ylimu = 20000, mymain = rownames(oligotypes)[1]); abline(h=global)
betta_pic(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], x = oligotypes[2, ], ylimu = 20000); abline(h=global)
betta_pic(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], x = oligotypes[3, ], ylimu = 20000); abline(h=global)
betta_pic(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], x = oligotypes[4, ], ylimu = 20000); abline(h=global)

betta(richness_estimates[names(oligotypes)], richness_errors[names(oligotypes)], X = cbind(1, c(oligotypes[4, ]/apply(oligotypes, 2, sum))))

#################### Hill numbers ####################
hill <- function(data, alpha) {
  if (mean(data[data>0]) > 1) {
    data <- data/sum(data)
  }
  data <- data[data>0]
  if (length(alpha) == 1) {
    ifelse(alpha == 0, length(data), 
           ifelse(alpha == 1, exp(-sum(data*log(data))), 
                  (sum(data^alpha))^(1 /(1-alpha))))
  } else {
    blank <- rep(NA, length(alpha))
    for (i in 1:length (alpha)) blank[ i] <- hill(data, alpha[i])
    blank
  }
}

alpha_collection <- seq(length.out = 20, from = 0, to = 2)
hill_numbers <- sapply(X = otus, FUN = hill, alpha = alpha_collection)
sapply(X = otus, FUN = hill, alpha = 1.0001)

plot(alpha_collection, hill_numbers[, 1], type="n", ylim=c(0, 500))
for (i in 1:dim(hill_numbers)[2]) {
  points(x=alpha_collection, hill_numbers[, i], type='l', col=ifelse(covariates$bloom2[i] == "yes", "red", "blue"))
}

### richness doesn't change but even this does?
#pdf("Hill_by_bloom.pdf", height=5,  width=7)
threshold <- 0.25
plot(alpha_collection, hill_numbers[, 1], type="n", ylim=c(0, 500), ylab="Sample Hill number(q)", xlab="q: 0 = Richness, 1 = Simpson, 2 = Shannon", main=paste("Hill numbers by bloom status: ", 100*(1-threshold), "% intervals", sep=""))
points(alpha_collection, apply(hill_numbers[, covariates$bloom2 == "yes"], 1, quantile, threshold), type="l", col="red", lty=2)
points(alpha_collection, apply(hill_numbers[, covariates$bloom2 == "yes"], 1, quantile, 0.5), type="l", col="red")
points(alpha_collection, apply(hill_numbers[, covariates$bloom2 == "yes"], 1, quantile, 1-threshold), type="l", col="red", lty=2)
points(alpha_collection, apply(hill_numbers[, covariates$bloom2 == "no"], 1, quantile, threshold), type="l", col="blue", lty=2)
points(alpha_collection, apply(hill_numbers[, covariates$bloom2 == "no"], 1, quantile, 0.5), type="l", col="blue")
points(alpha_collection, apply(hill_numbers[, covariates$bloom2 == "no"], 1, quantile, 1-threshold), type="l", col="blue", lty=2)
legend("top", c("Bloom", "No bloom"), col=c("red", "blue"), lty=1, bty="n")
#dev.off()

### richness doesn't change but even this does?
plot(alpha_collection, hill_numbers[, 1], type="n", ylim=c(0, 500), ylab="Hill number(q)", xlab="q", main="Hill numbers by ")
threshold <- 0.25
hill_quintiles <- apply(hill_numbers, 2, quantile, c(threshold, 0.5, 1-threshold))


plot(alpha_collection, hill_numbers[, 1], type="n", ylim=c(0, 500))
hill_plotting_function <- function(y, name) {
  print(rownames(y))
  points(x=alpha_collection, y, type='l', col=ifelse(covariates$bloom2[names(y)] == "yes", "red", "blue"))
}

plot(alpha_collection, hill_numbers[, 1], type="n", ylim=c(0, 500))
apply(X = hill_numbers, 2, points, x=alpha_collection, type='l')

names(metadata)
metadata[, 7]

#save.image()

################## Dec 6 ########################## 

###### ###### ###### ###### ###### ###### ###### ###### 
###### Stacked bar plot with oligotype proportions
###### ###### ###### ###### ###### ###### ###### ###### 

#for (otus in c("otus_500")) {
Year <- meta$Years; names(Year) <- rownames(meta)
Year <- Year[colnames(my_microcystis_counts)]; names(Year) <- rownames(meta)
Month <- meta$Months; names(Month) <- rownames(meta)
Month <- Month[colnames(my_microcystis_counts)]; names(Month) <- rownames(meta)

my_microcystis_counts <- get(paste(otus, "_microcystis_counts", sep=""))
rownames(my_microcystis_counts) <- paste("M", 1:dim(my_microcystis_counts)[1], sep="")
my_dolichospermum_counts <- get(paste(otus, "_dolichospermum_counts", sep=""))
rownames(my_dolichospermum_counts) <- paste("D", 1:dim(my_dolichospermum_counts)[1], sep="")
#my_cyanobacteria_matrix <- (data.frame( rbind(my_microcystis_counts, my_dolichospermum_counts), Type = c(rep("M", dim(my_microcystis_counts)[1]), rep("D", dim(my_dolichospermum_counts)[1]))))
my_cyanobacteria_matrix <- data.frame(Sample= colnames(my_microcystis_counts), cbind(t(my_microcystis_counts), t(my_dolichospermum_counts)))#, Type = c(rep("M", dim(my_microcystis_counts)[1]), rep("D", dim(my_dolichospermum_counts)[1]))))
# my_cyanobacteria_matrix <- data.frame(Year = Year,
#                                       Month = Month,
#                                       Sample=colnames(my_microcystis_counts),
#                                       cbind(t(my_microcystis_counts), t(my_dolichospermum_counts)))#, Type = c(rep("M", dim(my_microcystis_counts)[1]), rep("D", dim(my_dolichospermum_counts)[1]))))

Date <- as.Date(colnames(my_microcystis_counts), "X%d.%m.%Y")
tmp <- as.POSIXlt(Date[(Date < "1900-01-01") ]); tmp$year <- tmp$year+2000; Date[(Date < "1900-01-01") ] <- format(tmp)

my_cyanobacteria_matrix <- my_cyanobacteria_matrix[order(Date), ]

colnames(my_microcystis_counts)
Year = Year[colnames(my_microcystis_counts)]
length()
colnames(my_microcystis_counts) ## 150
rownames(meta) ## 143
meta$Years[colnames(my_microcystis_counts)]
colnames(my_microcystis_counts)

###### ###### ###### ###### 
###### everyone
###### ###### ###### ###### 

#### set up basic df and plot
my_cyanobacteria_matrix_renamed <- my_cyanobacteria_matrix
cbind(rownames(my_cyanobacteria_matrix), as.character(Date[order(Date)])) ## yes, matching
df <- melt(my_cyanobacteria_matrix_renamed, id = "Sample")

g <- ggplot(df, aes(x=Sample, y = value, fill = variable)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_x_discrete(labels=as.character(Date[order(Date)]))
g + geom_bar(stat='identity', position = "fill")

#### incorporate bloom info
my_cyanobacteria_matrix_renamed$Sample
bloom_df <- data.frame("ID"=meta$SampleID, "bloom"=meta$bloom2, "XID"=paste("X", meta$SampleID, sep=""))
intersect(bloom_df$XID, my_cyanobacteria_matrix_renamed$Sample)
strsplit(intersect(bloom_df$XID, my_cyanobacteria_matrix_renamed$Sample), split = "X")

bloom_df_useable <- data.frame("Sample"=my_cyanobacteria_matrix_renamed$Sample, "bloom"=NA)
for (i in 1:length(my_cyanobacteria_matrix_renamed$Sample)) {
  j <- my_cyanobacteria_matrix_renamed$Sample[i]
  if (j %in% intersect(bloom_df$XID, my_cyanobacteria_matrix_renamed$Sample)) {
    bloom_df_useable$bloom[i] <- meta$bloom2[which(as.character(j) == bloom_df$XID)] == "yes"
  }
}
bloom_df_useable

###### All nodes with bloom info
df_bloom_all <- data.frame(df, "bloom"=bloom_df_useable$bloom)
g_all <- ggplot(df_bloom_all[!is.na(df_bloom_all$bloom), ], aes(x=Sample, y = value, fill = variable, bloom)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_x_discrete(labels=as.character(Date[order(Date)])) +
  scale_alpha_discrete(range=c(0.5, 0.9))

pdf("AmyPlots/all-proportions-barplot.pdf", width = 20, height=5)
g_all + geom_bar(stat='identity', position = "fill", aes(alpha = bloom)) +
  scale_fill_manual(values = c("#990000", "#FF6633", "#FF6666", "#FF3333", "#CC6666", #reds 
                               "#3399FF", "#0066CC", "#003366", "#3399CC", "#FF66FF", 
                               "#0033FF", "#3366FF", "#CC66FF", "#3333FF", "#6666FF", 
                               "#CC0099", "#0000CC", "#6600CC", "#3300FF", "#0099FF"))
dev.off()

###### M nodes
df_m <- df_bloom_all[substr(df_bloom_all$variable,1,1) == "M", ]
g_m <- ggplot(na.omit(df_m), aes(x=Sample, y = value, fill = variable, bloom)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_x_discrete(labels=as.character(Date[order(Date)])) +
  scale_alpha_discrete(range=c(0.6, 0.9))

pdf("AmyPlots/m-proportions-barplot.pdf", width = 20, height=5)
g_m + geom_bar(stat='identity', position = "fill", aes(alpha = bloom)) + scale_fill_discrete(h=c(0,80))
dev.off()

###### D nodes
df_d <- df_bloom_all[substr(df_bloom_all$variable,1,1) == "D", ]
g_d <- ggplot(na.omit(df_d), aes(x=Sample, y = value, fill = variable, bloom)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_x_discrete(labels=as.character(Date[order(Date)])) +
  scale_alpha_discrete(range=c(0.6, 0.9))

pdf("AmyPlots/d-proportions-barplot.pdf", width = 20, height=5)
g_d + geom_bar(stat='identity', position = "fill", aes(alpha = bloom)) + scale_fill_discrete(h=c(220,280))
dev.off()

###### M vs D nodes
df <- melt(my_cyanobacteria_matrix_renamed, id = "Sample")
df_bloom_all <- data.frame(df, "bloom"=bloom_df_useable$bloom)
df_bloom_all_combine <- data.frame(df_bloom_all, "cyano"=substr(df_bloom_all$variable,1,1))
df_md <- df_bloom_all_combine %>% group_by(Sample, cyano, bloom) %>% summarise(total = sum(value))

g_md <- ggplot(na.omit(df_md), aes(x=Sample, y = total, fill = cyano, bloom)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_x_discrete(labels=as.character(Date[order(Date)])) +
  scale_alpha_discrete(range=c(0.6, 0.9))

pdf("AmyPlots/md-proportions-barplot.pdf", width = 20, height=5)
g_md + geom_bar(stat='identity', position = "fill", aes(alpha = bloom)) + scale_fill_discrete(h=c(250, 0))
dev.off()

## diversity analysis
i=1
df.richness  <- data.frame("sample"=my_cyanobacteria_matrix$Sample,"m"=NA, "mse"=NA,"d"=NA,"dse"=NA)
df.richness
for (i in 1:dim(my_cyanobacteria_matrix)[1]) {
  # m diversity
  ms <- my_cyanobacteria_matrix[i, 2:7]
  ds <- my_cyanobacteria_matrix[i, 8:32]
  if (any(ms > 0)) {
    freq_m <- make_frequency_count_table(ms)
    df.richness$m[i] <- sum(freq_m[,2])
    df.richness$mse[i] <- 1e-3
    
  } else{
    df.richness$m[i] <- 0
    df.richness$mse[i] <- 1e-3
    
  }
  if (any(ds >0)) {
    freq_d <- make_frequency_count_table(ds)
    df.richness$d[i] <- sum(freq_d[,2])
    df.richness$dse[i] <- 1e-3
    
  } else{
    df.richness$d[i] <- 0
    df.richness$dse[i] <- 1e-3
    
  }
  # if (length(freq_m)[1] > 5) {
  #   tmp <- freq_m ## catchall is really fussy! Need numbers not strings/factors
  #   write.table(tmp,"tmp.csv", row.names = FALSE,col.names = FALSE)
  #   system("./runcatchall.sh")
  # } else {
  # }
}

df.richness  <- data.frame("richness"=rep(NA,2*dim(my_cyanobacteria_matrix)[1]), "se"=NA, "type"=NA)
df.richness
j=1
for (i in 1:dim(my_cyanobacteria_matrix)[1]) {
  # m diversity
  ms <- my_cyanobacteria_matrix[i, 2:7]
  ds <- my_cyanobacteria_matrix[i, 8:32]
  
  df.richness$type[j] <- "M"
  if (any(ms > 0)) {
    freq_m <- make_frequency_count_table(ms)
    df.richness$richness[j] <- sum(freq_m[,2])
    df.richness$se[j] <- 1e-3
    j <- j + 1
  } else {
    df.richness$richness[j] <- 0
    df.richness$se[j] <- 1e-3
    j <- j + 1
  }
  
  df.richness$type[j] <- "D"
  if (any(ds >0)) {
    freq_d <- make_frequency_count_table(ds)
    df.richness$richness[j] <- sum(freq_d[,2])
    df.richness$se[j] <- 1e-3
    j <- j + 1
  } else{
    df.richness$richness[j] <- 0
    df.richness$se[j] <- 1e-3
    j <- j + 1
  }
}
df.richness


betta(df.richness$richness, df.richness$se, cbind(1,ifelse(df.richness[,3]=="M",1,0)))$table


dim(my_cyanobacteria_matrix)

apply(my_cyanobacteria_matrix[,2:32], 1, function(x) sum(x>0))

#save.image("/Users/adw96/Dropbox/LakeChamplain/AmyRScripts/workspace_161206.RData")

## 
load("/Users/adw96/Dropbox/LakeChamplain/AmyRScripts/workspace_161222.RData")
getwd()
Date
max(Date)
