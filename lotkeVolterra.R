library(stringr)

### Graph temporal series ###
X1 <- corr.filt %>% filter(X2 == "ASV_260") %>%
  select(X1) 
X1rm <- str_remove(X1$X1, "vir_")
X1rm

corr.ts[,!(names(corr.ts) %in% rm.col)]

head(vir.corr <- ASV_count[(rownames(ASV_count) %in% X1rm),])
tvir.corr <- as.data.frame(t(colsamp2date(vir.corr))) #apply colsamp2date custom function (no_transf.R) to change sample col names to dates

cyano.corr <- cyano_counts[c("ASV_260"),]
head(tcyano.corr <- t(cyano.corr))

head(corr.ts <- merge(tcyano.corr, tvir.corr, by="row.names"))
rownames(corr.ts) <- corr.ts[,1] #set col1 as rownames
corr.ts$Row.names <- sub("^([^.]*.[^.]*.[^.]*).*$",'\\1', corr.ts$Row.names) #rm everything after third period: 3x ".[^.]*" specifies 3 periods
names(corr.ts)[names(corr.ts) == "Row.names"] <- "Date"
corr.ts$Date <- as.Date(as.character(corr.ts$Date), "%d.%m.%Y") 
### add all dates
timeser <- merge(data.frame(Date = as.Date(min(corr.ts$Date):max(corr.ts$Date), "1970-1-1")), corr.ts, by = "Date", all = T)



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