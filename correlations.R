#Correlation Species-Env
spec <- abund_clean_env_hel[,2:68]
env <- abund_clean_env_hel[,71:86]
# Define data sets to cross-correlate
x <- as.matrix(spec)
y <- as.matrix(env)
# Cross correlate data sets
correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
# Or, alternatively, the same output is also available in a handy table format
correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
kable(head(correlation.table))


#### Hellinger Transformed & filtered ####
#transformed viral and cyano table with hellinger:
tcyano_helli_filt <- t(cyano_helli_filt)
dim(tcyano_helli_filt)
tvir_helli_filt <- t(vir_helli_filt)
dim(tvir_helli_filt)

#set to same dims
vir_hel_filt<- tvir_helli_filt[rownames(tvir_helli_filt) %in% rownames(tcyano_helli_filt),]
dim(vir_hel_filt)
#vir_hel_filt <- read.csv("vir_hel_filt.csv", header = T, row.names = 1)
colnames(vir_hel_filt) <- paste0("vir_", colnames(vir_hel_filt))

cyano_hel_filt <- tcyano_helli_filt[rownames(tcyano_helli_filt) %in% rownames(tvir_helli_filt),]
dim(cyano_hel_filt)
#cyano_hel_filt <- read.csv("cyano_hel_filt.csv", header=T, row.names = 1)

#https://microbiome.github.io/tutorials/Heatmap.html
helli_filt_corr <- associate(vir_hel_filt, cyano_hel_filt, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
head(helli_filt_corr)
heat(helli_filt_corr)

corr.filt <- helli_filt_corr %>% filter(Correlation > 0.6)
corr.filt[order(corr.filt$X2),]

unique(corr.filt$X1)
unique(corr.filt$X2)

### Graph temporal series ###
X1 <- corr.filt %>% filter(X2 == "ASV_260") %>%
  select(X1) 
X1rm <- str_remove(X1$X1, "vir_")
X1rm

corr.ts[,!(names(corr.ts) %in% rm.col)]

head(vir.corr <- ASV_count[(rownames(ASV_count) %in% X1rm),])
tvir.corr <- t(samp_date(vir.corr)) #apply samp_date custom function to change sample col names to dates

cyano.corr <- cyano_counts[c("ASV_260"),]
head(tcyano.corr <- t(cyano.corr))

head(corr.ts <- merge(tcyano.corr, tvir.corr, by="row.names"))
rownames(corr.ts) <- corr.ts[,1] #set col1 as rownames
orgDate <- corr.ts[,1]
orgDate2 <- sub("^(.*)[.].*", "\\1", orgDate[-c(1,3)]) #remove everything after last period. 1st and 3rd entries don't have same dims so omit
orgDate[-c(1,3)] <- orgDate2 
corr.ts[,1] <- orgDate

#break up into own cols
for (i in 1:nrow(corr.ts)){
  corr.ts$day[i] <- str_extract_all(corr.ts$Row.names, "[^.]+")[[i]][[1]]
  corr.ts$month[i] <- str_extract_all(corr.ts$Row.names, "[^.]+")[[i]][[2]]
  corr.ts$year[i] <- str_extract_all(corr.ts$Row.names, "[^.]+")[[i]][[3]]
}

corr.ts[,1] <- NULL #remove col1 can also call "timeseriesdf$Row.names"

corr.ts$date <- as.Date(with(corr.ts, paste(year, month, day, sep="-")), "%Y-%m-%d")
corr.ts <- corr.ts[order(as.Date(corr.ts$date, format="%Y-%m-%d")),] #order by date
rm.col <- c("day", "month", "year", "date")
# corr.ts <- corr.ts[,!(names(corr.ts) %in% rm.col)] #keep only ASV cols
ts.corr <- log(corr.ts[,!(names(corr.ts) %in% rm.col)])
head(ts.corr)
dim(ts.corr)

names(ts.corr)[names(ts.corr) == "ASV_260"] <- "Cyanobacteria (ASV_260)"

#plot
ts.corr.plot <- ts.corr %>% 
  rownames_to_column() %>% 
  gather(key = key, value = value, "Cyanobacteria (ASV_260)":ASV_447) %>% 
  mutate(rowname = factor(rowname))

head(ts.corr.plot)

#organize doliASV50 last and viral ASV of interest second to last, so it's plotted line is brought to the front on graph
ts.corr.plot$key <- factor(ts.corr.plot$key, c(
                                               "ASV_447", 
                                               "Cyanobacteria (ASV_260)"))
ts.corr.plot %>%
  ggplot(aes(x = as.numeric(rowname), y = value, color = key)) + 
  geom_point() +
  geom_line() +
  ggtitle("Timeseries: Cyanobacteria (ASV_260) and viral ASV_443")+
  scale_x_discrete(labels = corr.ts$date, name="Date")+#change x-axis sample name to Month
  scale_y_continuous(name = "log(abondance)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  scale_color_manual(values=c(
                              "ASV_447" = "blue",  
                              "Cyanobacteria (ASV_260)" = "red"),
                     breaks = c("ASV_447",
                       "Cyanobacteria (ASV_260)"
                     ))+ #ensures legend stays in same order
  theme_bw()


# #organize doliASV50 last and viral ASV of interest second to last, so it's plotted line is brought to the front on graph
# ts.corr.plot$key <- factor(ts.corr.plot$key, c("ASV_380",
#                                                "ASV_636", 
#                                                "ASV_347", 
#                                                "ASV_289", 
#                                                "ASV_346",
#                                                "ASV_261",
#                                                "ASV_146", 
#                                                "Dolichospermum (ASV_50)"))
# ts.corr.plot %>%
#   ggplot(aes(x = as.numeric(rowname), y = value, color = key)) + 
#   geom_point() +
#   geom_line() +
#   ggtitle("Timeseries: Dolichospermum ASV_50 and viral ASV_146")+
#   scale_x_discrete(labels = corr.ts$date, name="Date")+#change x-axis sample name to Month
#   scale_y_continuous(name = "log(abondance)")+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #rotate axis labels
#        plot.title = element_text(hjust = 0.5))+ #center title
#   scale_color_manual(values=c("ASV_636" = "lightgrey", 
#                               "ASV_347"= "lightgrey", 
#                               "ASV_289"= "lightgrey",
#                               "ASV_346"= "lightgrey", 
#                               "ASV_146"= "blue", 
#                               "ASV_261"= "lightgrey", 
#                               "ASV_380"= "lightgrey",
#                               "Dolichospermum (ASV_50)"= "red"),
#                      breaks = c(
#                        "ASV_146", 
#                        "ASV_261", 
#                        "ASV_289",
#                        "ASV_346",  
#                        "ASV_347", 
#                        "ASV_380",
#                        "ASV_636",
#                        "Dolichospermum (ASV_50)"
#                      ))+ #ensures legend stays in same order
#   theme_bw()








### Compare alpha diversities of cyano and phage ###
dim(ba_bact_df)
dim(ba_vir_df)
head(ba_bact_df)
head(ba_vir_df)

#merge dfs to ensure same dims
df <- merge(ba_bact_df,ba_vir_df, by="sample")
dim(df)
head(df)

# #check that ba_bact_df is .x
# which(ba_bact_df$sample == "01.06.2008") 
# ba_bact_df[37,]

head(alpha_bact <- df[,1:2])
head(alpha_vir <- df[,c(1, 7)])

#cross correlation and lagged regressions
alpha_corr <- ccf(alpha_bact$richness.x, alpha_vir$richness.y)
alpha_corr

alphadiv.corr <- df[,c("richness.x", "richness.y")]
rownames(alphadiv.corr) <- df$sample

names(alphadiv.corr)[names(alphadiv.corr) == "richness.x"] <- "bact.rich"
names(alphadiv.corr)[names(alphadiv.corr) == "richness.y"] <- "vir.rich"
head(alphadiv.corr)

#http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
library("ggpubr")

### Is the covariation linear? 
ggscatter(alphadiv.corr, x = "bact.rich", y = "vir.rich", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Bacterial richness", ylab = "Viral richness", title = "Correlation between viral and bacterial alpha diversity")
#Yes, from the plot the relationship is linear. 
#situations where the scatter plots show curved patterns, dealing with nonlinear association b/w the 2 variables.

### Are the data from each of the 2 variables (x, y) follow a normal distribution?

##Shapiro-Wilk test can be performed as follow:
#Null hypothesis: the data are normally distributed
#Alternative hypothesis: the data are not normally distributed

# Shapiro-Wilk normality test for bact
shapiro.test(alphadiv.corr$bact.rich) # => p = 0.1102
# Shapiro-Wilk normality test for vir
shapiro.test(alphadiv.corr$vir.rich) # => p = 0.9555
#p-values >0.05 implying that the distribution of the data are not significantly different from normal distribution. 
#In other words, we can assume the normality.

#Visual inspection of the data normality 
ggqqplot(alphadiv.corr$bact.rich, ylab = "bact.rich")
ggqqplot(alphadiv.corr$vir.rich, ylab = "vir.rich")

richness.corr <- cor.test(alphadiv.corr$bact.rich, alphadiv.corr$vir.rich, method = "pearson")
richness.corr
# p-val = 0.6889 > 0.05, therefore bact rich and vir rich are not signif correlated with a corr coeff of -0.0363

