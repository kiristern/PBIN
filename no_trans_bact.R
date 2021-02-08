#upload cyano ASV data
cyano_counts <- read.table("cyano/Champ_ASVs_counts.txt", header = TRUE, row.names = 1)
head(cyano_counts)
cyano_taxa <- read.csv("cyano/ASVs_taxonomy_Champ_Greengenes.csv", header = T, row.names = 1, fill=T)
head(cyano_taxa)

nrow(meta)
length(cyano_counts)

cyano_taxa$Species <- rownames(cyano_taxa)

colnames(cyano_counts)
#remove X at beginning of date
colnames(cyano_counts)[1:135] <- substring(colnames(cyano_counts)[1:135], 2)
#select dates cols only
cyano_counts <- cyano_counts[1:135]

#match sample dates
meta2 <- meta

#make sure meta matches cyano samples
nrow(meta2)
length(cyano_counts)

#select cols that match dates
bact_counts <- cyano_counts[,(colnames(cyano_counts) %in% meta2$description)]
length(bact_counts)

meta2 <- meta2[meta2$description %in% colnames(bact_counts),]
nrow(meta2)
rownames(meta2) <- meta2$description


#Phyloseq
library(phyloseq)

bac_count <- otu_table(bact_counts, taxa_are_rows = T)
dim(bac_count)
bact_tax_tab <- tax_table(cyano_taxa)
dim(cyano_taxa_ps)
rownames(bact_tax_tab) <- rownames(cyano_taxa)

#add to phyloseq object
sample_info_cyano <- sample_data(meta2)
dim(sample_info_cyano)

bact_physeq <- phyloseq(bac_count, bact_tax_tab, sample_info_cyano)
print(bact_physeq)

#rename cols
colnames(tax_table(bact_physeq)) <- c("Kingdom", "Phylum", "Class",
                                      "Order", "Family", "Genus", "ASV")
#quick check
bact_physeq %>% tax_table %>% head()
bps <- bact_physeq %>% tax_table
bps[260,]

bact_abun <- bact_physeq %>% otu_table()

#reorder phyloseq by chronological date
all(colnames(bact_counts) %in% rownames(meta2))
map <- meta2[order(meta2$Date),]
toorder <- rownames(map)
otu_table(bact_physeq) <- otu_table(bact_physeq)[,toorder]
bact_physeq %>% otu_table()

cyano_ps <- subset_taxa(bact_physeq, Phylum == "p__Cyanobacteria")


# ALPHA DIV #
#richness by year
ba_bact <- breakaway(bact_physeq)
ba_bact

ba_bact_df = data.frame("richness" = (ba_bact %>% summary)$estimate,
                        "sample" = (ba_bact %>% summary)$sample_names,
                        "error" = (ba_bact %>% summary)$error,
                     "Years" = bact_physeq %>% sample_data %>% get_variable("Years"),
                     "Upper" = (ba_bact %>% summary)$upper,
                     "Lower" = (ba_bact %>% summary)$lower)
head(ba_bact_df)
str(ba_bact_df)
ggplot(ba_bact_df, aes(x = fct_inorder(sample), y = richness, color = Years))+ #fct_inorder ensures plotting in order of sample date
        geom_point()+
        geom_pointrange(aes(ymin=richness-abs(richness-Lower), ymax=richness+abs(richness-Upper)))+ #vs. geom_errorbar
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
        ggtitle("Breakaway richness by sample")+
        scale_x_discrete(labels = bact_physeq %>% sample_data %>% get_variable("Months"), name="Month")#change x-axis sample name to Month

# # alpha diversity for (filtered) cyano only
# ba_cyano <- breakaway(cyano_ps)
# ba_cyano
# 
# ba_cyano_df = data.frame("richness" = (ba_cyano %>% summary)$estimate,
#                         "sample" = (ba_cyano %>% summary)$sample_names,
#                         "error" = (ba_cyano %>% summary)$error,
#                         "Years" = cyano_ps %>% sample_data %>% get_variable("Years"),
#                         "Upper" = (ba_cyano %>% summary)$upper,
#                         "Lower" = (ba_cyano %>% summary)$lower)
# head(ba_cyano_df)
# str(ba_cyano_df)
# ggplot(ba_cyano_df, aes(x = fct_inorder(sample), y = richness, color = Years))+ #fct_inorder ensures plotting in order of sample date
#   geom_point()+
#   geom_pointrange(aes(ymin=richness-abs(richness-Lower), ymax=richness+abs(richness-Upper)))+ #vs. geom_errorbar
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
#         plot.title = element_text(hjust = 0.5))+ #center title
#   ggtitle("Breakaway richness by sample")+
#   scale_x_discrete(labels = cyano_ps %>% sample_data %>% get_variable("Months"), name="Month")#change x-axis sample name to Month

  

#boxplot years
(ba_plot <-  ggplot(ba_bact_df, aes(x = Years, y = richness))+
    geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                                 geom="crossbar", width=0.5) + theme_minimal()+
  ggtitle("Observed richness by year")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_y_continuous(name = "Observed richness")





### Divnet ###
## R crashes when using viral_physeq (data too large??)
#check all variables in filtered phyloseq object
cyano_ps <- subset_taxa(bact_physeq, Phylum == "p__Cyanobacteria")
filt_cyaps <- filter_taxa(cyano_ps, function(x) sum(x > 1) > (0.10*length(x)), TRUE) ##remove taxa not seen more than 1 times in at least 10% of the samples. This protects against ASV with small mean & trivially large coef of var

toptaxa <- find.top.taxa(filt_cyaps,"ASV")
head(toptaxa)
tt <- toptaxa %>% select(c("ASV"))
#see which taxa comes up most
tt <- as.data.frame(table(tt))
tt <- tt %>% arrange(desc(Freq))
head(tt)



#TODO: make new for bact and viral divnets with the same samples/dates (taken from mantel)




sample_variables(filt_cyaps)

ASV11_years_cyan <- filt_cyaps %>%
  divnet(X = "Years", ncores = 4, 
         base = "ASV_11")
ASV11_years_cyan

ASV11_years_cyan$shannon %>% head

#to test if alpha-diversity (by default, Shannon) is equal across the values of the covariate X:
testDiversity(ASV11_years_cyan)

filtcyan <- filt_cyaps %>% otu_table()
filtcyan <- t(filtcyan)
#isolate for date only
df.filtcyan <- as.data.frame(row.names(filtcyan))
df.filtcyan$date <- gsub("^([^_]*_[^_]*_[^_]*_[^_]*)_.*$", "\\1", df.filtcyan[,1]) #removes everything after last _
df.filtcyan$date <- sub("*._*._*._*._*._*._*._","", df.filtcyan$date) #remove sample ID at beginning

#compare the plug-in Shannon with divnet estimates
library(ggplot2)
ASV11_years_cyan$shannon %>%
  plot(filt_virseq, color = "Years") +
  scale_x_discrete(labels = df.filtcyan$date, name="Sample date")+ #change x-axis sample name to date
  ylab("Shannon diversity estimate\n(ASV1 level)")+
  ggtitle("Shannon diversity estimate between years\n(base ASV1)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), #rotate axis labels
        plot.title = element_text(hjust = 0.5)) #center title
#only a single DivNet estimate for each year (along with error bars). For characteristics for which many samples were observed, there are smaller error bars than for samples for which there was only one sample (seems reasonable -- we had less data).




#distribution of Bray-Curtis distances between the samples
simplifyBeta(ASV11_years_cyan, filt_cyaps, "bray-curtis", "Years") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est,
             col = interaction(Covar1, Covar2))) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

merge_samples(ASV11_years_cyan, "Years") %>%
  sample_shannon %>%
  plot()

#Shannon index using breakaway
estimates <- ASV11_years_cyan$shannon %>% summary %>% select("estimate")
ses <- sqrt(ASV11_years_cyan$`shannon-variance`)
X <- breakaway::make_design_matrix(filt_cyaps, "description")
(ba_shannon <- betta(estimates, ses, X)$table)

#write.csv(ba_shannon, "ba_shannonDescription_vir.csv")

