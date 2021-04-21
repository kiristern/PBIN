library(stringr)

### Graph temporal series ###
X1 <- corr.filt %>% filter(X2 == "ASV_50") %>%
  select(X1) 
X1rm <- str_remove(X1$X1, "vir_")
X1rm

head(vir.corr <- ASV_count[(rownames(ASV_count) %in% X1rm),])
tvir.corr <- t(colsamp2date(vir.corr)) #apply colsamp2date custom function (no_transf.R) to change sample col names to dates

cyano.corr <- cyano_counts[c("ASV_50"),]
head(tcyano.corr <- t(cyano.corr))

head(corr.ts <- merge(tcyano.corr, tvir.corr, by="row.names"))
rownames(corr.ts) <- corr.ts[,1] #set col1 as rownames
corr.ts$Row.names <- sub("^([^.]*.[^.]*.[^.]*).*$",'\\1', corr.ts$Row.names) #rm everything after third period: 3x ".[^.]*" specifies 3 periods
names(corr.ts)[names(corr.ts) == "Row.names"] <- "Date"
corr.ts$Date <- as.Date(as.character(corr.ts$Date), "%d.%m.%Y") 
### add all dates
timeser <- merge(data.frame(Date = as.Date(min(corr.ts$Date):max(corr.ts$Date), "1970-1-1")), corr.ts, by = "Date", all = T)
head(timeser)
tail(timeser)

names(timeser)[names(timeser) == "ASV_50"] <- "Dolichospermum (ASV_50)"

timeser.keep <- timeser[, c("Dolichospermum (ASV_50)", "ASV_636")]
timeser.keep$Date <- timeser$Date
write.csv(timeser.keep, "LVdf.csv")

timeser.log <- log(timeser.keep)

#plot
ts.corr.plot <- timeser.log %>% 
  rownames_to_column() %>% 
  gather(key = key, value = value, "Dolichospermum (ASV_50)":ASV_636) %>% 
  mutate(rowname = factor(rowname))

head(ts.corr.plot)

ts.corr.plot$key <- factor(ts.corr.plot$key, c(
  "ASV_636", 
  "Dolichospermum (ASV_50)"))
ts.corr.plot %>%
  ggplot(aes(x = as.numeric(rowname), y = value, color = key)) + 
  geom_point() +
  geom_line() +
  ggtitle("Timeseries: Dolichospermum (ASV_50) and viral ASV_636")+
  #scale_x_discrete(labels = timeser$Date, name="Date")+#change x-axis sample name to Month
  scale_y_continuous(name = "log(abondance)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #rotate axis labels
        plot.title = element_text(hjust = 0.5))+ #center title
  scale_color_manual(values=c(
    "ASV_636" = "blue",  
    "Dolichospermum (ASV_50)" = "red"),
    breaks = c("ASV_636",
               "Dolichospermum (ASV_50)"))+ #ensures legend stays in same order
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