setwd("~/Documents/GitHub/PBIN/data")

#####Abundance#####
# df <- read.csv("FINAL_ASV_TAXA_ABUN.csv")
#TODO: make sure col name is "abundance"

# TODO: RUN ABUNDANCE_CLEANING.R SCRIPT
df %>% group_by(samples, ASV_ID) %>% summarise(Relative_Abundance = sum(abundance)) %>%
  ggplot(aes(x = samples, y = Relative_Abundance, fill = ASV_ID)) + 
  geom_bar(stat = "identity", width=1, show.legend = F) +
  theme_minimal()

# df_unk <- read.csv("rel_ab_top20_unknown_FINAL.csv")
# df_unk %>% group_by(samples, ASV_ID) %>% summarise(Relative_Abundance = sum(abundance)) %>%
#   ggplot(aes(x = samples, y = Relative_Abundance, fill = ASV_ID)) + 
#   geom_bar(stat = "identity", show.legend = F)

cya <- read.csv("cyano/cyano_samples.csv")
mic <- read.csv("cyano/micro_samples.csv")
dol <- read.csv("cyano/doli_samples.csv")







##ordination
##https://joey711.github.io/phyloseq/plot_ordination-examples.html

#just plotting asv
GP.ord <- ordinate(viral_physeq, "NMDS", "bray")
p1 = plot_ordination(viral_physeq, GP.ord)
print(p1)


##using breakaway to explore alpha diversity
##https://github.com/adw96/stamps2018/blob/master/estimation/diversity-lab.R

viral_physeq
viral_physeq %>% sample_data

#how to extract only values corresponding to period "spring" "summer" "fall"
# season <- viral_physeq %>% 
#   subset_samples(Period %in% c("Spring", "Summer", "Fall"))

#look at observed richness plot colour by month
observed <- sample_richness(viral_physeq)
summary(observed)
plot(observed, viral_physeq, color="Years")

#get linear regression
alpha_df <- estimate_richness(viral_physeq)
alpha_df <- as.data.frame(alpha_df$Observed)
alpha_df <- tibble::rowid_to_column(alpha_df, "Sample")
names(alpha_df)[2] <- "richness"
plot(alpha_df)
abline(lm(richness ~ Sample, data=alpha_df))
#get adjusted R2
alphadf_summary <- summary(lm(richness~Sample, data = alpha_df))
alpha_r2 <- alphadf_summary$adj.r.squared
#get p-value
alpha_p <- alphadf_summary$coefficients[2,4]
#add to graph
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(alpha_r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(alpha_p, digits = 2)))[2]
legend("topleft", legend = rp, bty = "n")

#find linear regression formula
lm(richness ~ Idx, data = alpha_df)
#y = 0.9103x + 188.7873


#observed richness 
#Sequencing depth by Year
# data.frame("observed_richness" = (observed %>% summary)$estimate,
#            "Depth" = phyloseq::sample_sums(viral_physeq), # sequence depth
#            "type" = viral_physeq %>% sample_data %>% get_variable("Years")) %>%
#   ggplot(aes(x = Depth, y = observed_richness, color = type)) +
#   geom_point()
# 
# #Year by month
# data.frame("observed_richness" = (observed %>% summary)$estimate,
#            "Years" = viral_physeq %>% sample_data %>% get_variable("Years"),
#            "type" = viral_physeq %>% sample_data %>% get_variable("Months")) %>%
#   ggplot(aes(x = Years, y = observed_richness, color = type)) +
#   geom_point()
# 
# #Month by year
# data.frame("observed_richness" = (observed %>% summary)$estimate,
#            "Months" = viral_physeq %>% sample_data %>% get_variable("Months"),
#            "type" = viral_physeq %>% sample_data %>% get_variable("Years")) %>%
#   ggplot(aes(x = Months, y = observed_richness, color = type)) +
#   geom_point()

#estimate number of missing species using a species richness estimate
ba <- breakaway(viral_physeq)
ba
plot(ba, viral_physeq, color="Years") 
#y = 0.9103x + 188.7873

#get sample names
# x <- viral_physeq %>% sample_data %>% rownames()

#table of sample names, richness, and error 
# alpha_err <- data.frame(x,
                        # summary(ba)$estimate,
                        # summary(ba)$error,
                        # make_design_matrix(viral_physeq, "Years"))

# richness <- summary(ba)$estimate
# std_dev <- summary(ba)$error
# 
# plot(x, richness,
#      ylim = range(c(richness-std_dev, richness+std_dev)),
#      xlim = 
#      pch=19, xlab = "Sample", ylab = "Richness +/- Std Dev",
#      main = "Alpha div with std.dev error bars"
# )
# #hack: draw arrow heads as lines
# arrows(x, richness-std_dev, x, richness+std_dev, length=0.05, angle = 90, code=3)

#for variables Years and bloom2
ba_years = data.frame("ba_observed_richness" = (ba %>% summary)$estimate,
           "Years" = viral_physeq %>% sample_data %>% get_variable("Years"))
#write.table(ba_years, "ba_years.txt")
ba_years <- read.table("ba_years.txt", header = T)
fit <- lm(ba_observed_richness ~ Years, data = ba_years)
coefs <- coef(fit)
#get r2
(r2 = signif(summary(fit)$adj.r.squared))
#get p-value
(pval <- signif(summary(fit)$coef[2,4]))

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


# #geom_errorbar()
# ba_plot + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#                geom="errorbar", color="red", width=0.2) +
#   stat_summary(fun.y=mean, geom="point", color="red")


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
#subset data between bloom and non bloom samples
data_bloom <- subset_samples(viral_physeq, bloom2=="yes")
data_no_bloom <- subset_samples(viral_physeq, bloom2=="no")


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






