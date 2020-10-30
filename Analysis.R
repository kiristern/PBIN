setwd("~/Documents/GitHub/PBIN/data")

#####Abundance#####
# TODO: RUN ABUNDANCE_CLEANING.R SCRIPT or
# df <- read.csv("FINAL_ASV_TAXA_ABUN.csv")
 # df2 <- read.csv("asvID_filt_rem_rel_abun.csv", row.names = 1, header = T)
 # head(df2)
#TODO: make sure col name is "abundance"
head(df.l)
head(df.p)
# df %>% group_by(samples, ASV_ID) %>% summarise(Relative_Abundance = sum(abundance)) %>%
#   ggplot(aes(x = samples, y = Relative_Abundance, fill = ASV_ID)) + 
#   geom_bar(stat = "identity", show.legend = T) +
#   theme_minimal()

ggplot(df.l, aes(x = samples.l, y = abundance, fill = ASV_ID.l))+
  geom_bar(stat="identity", show.legend = T)+
  theme_minimal()

ggplot(df.p, aes(x = samples.p, y = abundance, fill = ASV_ID.p))+
  geom_bar(stat="identity", show.legend = T)+
  theme_minimal()


# asv_rel_abun
# library(reshape2)
# df_long <- melt(asv_rel_abun, id.vars=rownames(asv_rel_abun), variable.name = colnames(asv_rel_abun))




#Convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(filt_virseq, function(x){x / sum(x)})
phyloseq::otu_table(filt_virseq)[1:5, 1:5]
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5]
#Plot
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Status, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())







##ordination
##https://joey711.github.io/phyloseq/plot_ordination-examples.html

#just plotting asv
GP.ord <- ordinate(viral_physeq, "NMDS", "bray")
p1 = plot_ordination(viral_physeq, GP.ord)
print(p1)






##using breakaway to explore alpha diversity
#species richness estimate
ba <- breakaway(filt_virseq)
ba
plot(ba, filt_virseq, color="Years") 
#y = 0.9103x + 188.7873

#get sample names
x <- viral_physeq %>% sample_data %>% rownames()


#table of sample names, richness, and error 
 alpha_err <- data.frame(x,
                        summary(ba)$estimate, #richness
                        summary(ba)$error, #stnd dev
                        make_design_matrix(viral_physeq, "Years"))
alpha_err$count <- 1:nrow(alpha_err)
 
plot(ba)


#for variables Years and bloom2
#boxplot bloom
ba_bloom = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
                       "Bloom" = filt_virseq %>% sample_data %>% get_variable("bloom2"))
 
(ba_plot <-  ggplot(ba_bloom, aes(x = Bloom, y = ba_observed_richness))+
   geom_point()) + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                             geom="crossbar", width=0.5) + theme_minimal()

#boxplot years
ba_years = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
           "Years" = filt_virseq %>% sample_data %>% get_variable("Years"))
#write.table(ba_years, "ba_years_filt_virseq.txt")
ba_years <- read.table("ba_years_filt_virseq.txt", header = T)

#finding linear regression:
#http://r-statistics.co/Linear-Regression.html 
fit <- lm(ba_observed_richness ~ Years, data = ba_years)
coefs <- coef(fit)
summary(fit)

#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coef[2,4])) #wrong! coef[2,4] != p-val

(ba_plot <-  ggplot(ba_years, aes(x = Years, y = ba_observed_richness))+
  geom_point()+
    geom_abline(intercept = coefs[1], slope = coefs[2])+
    labs(title = paste("Adj R2 = ", r2,
                       #"Intercept =", signif(coefs[1]),
                       #"Slope =", signif(coefs[2]),
                       "p-value =", pval)))
#geom_crossbar()
ba_plot + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                       geom="crossbar", width=0.5) + theme_minimal()




#look at just one sample
tr <- viral_physeq %>% subset_samples(Date == "2006-5-30")
tr
#look at the structure of this dataset
freq_count <- tr %>% otu_table %>% make_frequency_count_table
#this is the frequency count table for this dataset
freq_count %>% head(10) #14 ASVs observed only twice  
freq_count %>% tail() #1 ASV observed 293 times

#fit breakaway to the sample
ba_tr <- breakaway(freq_count)
ba_tr #this is an alpha diversity estimate -- a special class for alpha div. estimates

#check which model breakaway picked
ba_tr$model

#take estimates and turn into df 
# summary(ba) %>%
#   add_column("SampleNames" = viral_physeq %>% otu_table %>% sample_names)

#chose a different species richness estimate
chao_est <- viral_physeq %>%
  chao1 %>%
  plot(viral_physeq, color="Months")



#betta() works like a regression model but accounts for the uncertainty in estimating diversity
# test hypothesis that different types of water systems have the same microbial div (not in this case)
bt <- betta(summary(ba)$estimate,
            summary(ba)$error,
            make_design_matrix(viral_physeq, "Site"))
bt$table
#betta() estimates that the mean ASV-level diversity in ... is 368?
#estimates the the diversity in March is significantly lower (~222) than all others but that they're all below 0??
#betta accounts for the error bars in diversity when doing hypothesis testing

#Shannon index: downweights the importance of rare taxa (reflects low abundance of a taxon)
#DivNet adjusts for different sequencing depths so don't need to rarefy (throw away data)

#Run in parallel
dv_viral_ps <- divnet(viral_physeq, ncores = 4)




######## Beta Diversity ##########
# #subset data between bloom and non bloom samples
# data_bloom <- subset_samples(viral_physeq, bloom2=="yes")
# data_no_bloom <- subset_samples(viral_physeq, bloom2=="no")

# PCoA plot using the unweighted UniFrac and JSD as distance
unifrac_dist <- phyloseq::distance(viral_physeq, method="unifrac", weighted=T)
jsd_dist <- sqrt(phyloseq::distance(viral_physeq, "jsd")) #jsd is a semi-metric
#PCoA bloom
unifrac_dist_bloom <- phyloseq::distance(data_bloom, method="unifrac", weighted=T)
jsd_dist_bloom <- sqrt(phyloseq::distance(data_bloom, "jsd")) #jsd is a semi-metric
#PCoA no bloom
unifrac_dist_nob <- phyloseq::distance(data_no_bloom, method="unifrac", weighted=T)
jsd_dist_nob <- sqrt(phyloseq::distance(data_no_bloom, "jsd")) #jsd is a semi-metric


#change distance=unifrac_dist to jsd_dist
ordination <- ordinate(viral_physeq, method="PCoA", distance=jsd_dist)
plot_ordination(viral_physeq, ordination, color="Years") + 
  #stat_ellipse(type = "norm")+#, linetype = 2) +
  #stat_ellipse(type = "t") +
  theme(aspect.ratio=1) + 
  theme_classic() #+
#theme_bw()
#plot PCoA bloom
ordination_bloom <- ordinate(data_bloom, method="PCoA", distance=jsd_dist_bloom)
plot_ordination(data_bloom, ordination, color="Years") + 
  #stat_ellipse(type = "norm")+#, linetype = 2) +
  #stat_ellipse(type = "t") +
  theme(aspect.ratio=1) + 
  theme_classic() #+
   #theme_bw()
#Plot PCoA no bloom
ordination_nob <- ordinate(data_no_bloom, method="PCoA", distance=jsd_dist_nob)
plot_ordination(data_no_bloom, ordination, color="Years") + 
  #stat_ellipse(type = "norm")+#, linetype = 2) +
  #stat_ellipse(type = "t") +
  theme(aspect.ratio=1) + 
  theme_classic() #+
#theme_bw()



#Test whether the Years/Months differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:
adonis(unifrac_dist ~ sample_data(viral_physeq)$Years)
adonis(unifrac_dist ~ sample_data(viral_physeq)$Months)

#Permanova test using adonis function
adonis(dist ~ variable, as(sample_data(data), "data.frame"))


#####Phyloseq analysis#####

#Transform to relative abundance
SmileT<- transform_sample_counts(filt_vir, function(x) x / sum(x) )

# PCoA plot using the unweighted UniFrac and JSD as distance
unifrac_dist <- phyloseq::distance(filt_vir_seq, method="unifrac", weighted=T)
jsd_dist <- sqrt(phyloseq::distance(filt_vir_seq, "jsd")) #jsd is a semi-metric

#ordination_betadiversity_PCOA
pcoa=ordinate(filt_vir_seq, "PCoA", distance=jsd_dist)
nmds=ordinate(filt_vir_seq,"NMDS",distance=jsd_dist)

#Permanova test using adonis function
adonis(jsd_dist ~ Years, as(sample_data(filt_vir_seq), "data.frame"))



#Dispersion analysis
betatax=betadisper(jsd_dist,data.frame(sample_data(filt_vir_seq))$Years)
p=permutest(betatax)
p$tab

#ANOSIM test
variable_group = get_variable(filt_vir_seq, "Years")
variable_group = anosim(jsd_dist, variable_group)





