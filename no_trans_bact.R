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


#Phyloseq
library(phyloseq)

bac_count <- otu_table(bact_counts, taxa_are_rows = T)
cyano_taxa_ps <- tax_table(cyano_taxa)
rownames(cyano_taxa_ps) <- rownames(cyano_taxa)

#add to phyloseq object
sample_info_cyano <- sample_data(meta2)
bact_physeq <- phyloseq(bac_count, cyano_taxa_ps, sample_info_cyano)
print(bact_physeq)

#rename cols
colnames(tax_table(bact_physeq)) <- c("Kingdom", "Phylum", "Class",
                                      "Order", "Family", "Genus", "ASV")
#quick check
bact_physeq %>% tax_table %>% head(n=50)


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

