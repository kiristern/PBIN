#stats test between box plots

#significance between bloom and no bloom richness
ba_alpha

#remove NA samples
bloom_nobloom <- ba_alpha[complete.cases(ba_alpha), ]

#compute summary statistics by groups
group_by(ba_alpha, Bloom) %>%
  summarise(
    count = n(),
    mean = mean(ba_observed_richness, na.rm = TRUE),
    sd = sd(ba_observed_richness, na.rm = TRUE)
  )

#visualize using box plots
#plot richness by bloom and colour by bloom
library(ggpubr)
ggboxplot(bloom_nobloom, x = "Bloom", y = "ba_observed_richness",
          color = "Bloom", palette = c("#00AFBB", "#E7B800"),
          ylab = "Observed Richness", xlab = "Bloom")


#tests to check independent t-test assumptions
#Assumption 1: are the two samples independent? Yes. not taken from the same time 
#Assumption 2: does the data from each of the 2 groups follow a normal distirbution?
  #Use Shapiro-Wilk normaility test
    #Null hypothesis: the data are normally distributed
    #Alternative hypothesis: the data are not normally distributed

# Shapiro-Wilk normality test for Bloom
with(bloom_nobloom, shapiro.test(ba_observed_richness[Bloom == "yes"])) # p = 0.2278
# Shapiro-Wilk normality test for No Bloom
with(bloom_nobloom, shapiro.test(ba_observed_richness[Bloom == "no"])) # p = 0.9858
#the two p-values are greater than the significance level 0.05 implying that the distribution of the 
#data are not significantly different from the normal distribution. Ie, we can assume the normality.
  #if the data are not normally distributed, it’s recommended to use the non parametric two-samples Wilcoxon rank test.

#Assumption 3: Do the two populations have the same variances?
  #use F-test to test for homogeneity in variances

(res.ftest <- var.test(ba_observed_richness ~ Bloom, data = bloom_nobloom))
  #p-value: 0.2611. It’s greater than the significance level alpha = 0.05. There is no significant 
  #difference between the variances of the two sets of data. Therefore, we can use the classic t-test witch assume equality of the two variances.

#Compute unpaired two-samples t-test
#Is there any significant difference between Bloom and no-bloom richness?
# Compute t-test
(res <- t.test(ba_observed_richness ~ Bloom, data = bloom_nobloom, var.equal = TRUE))
#p-value = 0.4782 > 0.05 therefore fail to reject null hypothesis. ie. richness between bloom/no bloom is the same
#aka no significant difference between groups






