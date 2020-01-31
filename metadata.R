library("dplyr")
library("tidyverse")

samples <- read.table("data/samples.txt")


#assign new col name
colnames(samples)[1] <- "SampleID"

#add cols for rest of meta data
samples[c("Day", "Month", "Year", "Date", "Site", "Period", "bloom2", "Total_Phosphorus_ug", 
        "Phosph_Range", "Nitrog_Range", "Total_Nitrogen_mg", "Temperature_Water_Celsius",
        "Dissolved_P", "Dissolved_N", "Cumulative_precipitation_t1_t7_mm", 
        "Profondeur_Secchi_cm", "Mean_temperature_t0_t7", "Microcystin_ug_L", "Description")] <- NaN

#fill Description col with values from SampleID col
samples["Description"] <- samples["SampleID"]

#remove everything before 4th period in sampleID (just to keep date)
samples$SampleID <- gsub(".*\\.","\\1", samples$SampleID, perl=T) 

#extract proper year
samples$Year <- gsub(".*-(.*)\\-.*", "\\1", samples$SampleID, perl=T)

#extract proper month
samples$Month <- gsub(".-*.-(.*)\\-.*-.*", "\\1", samples$SampleID, perl=T)

#extract day
samples$Day <- gsub("-.*", "\\1", samples$SampleID, perl=T)

#reorganize date format
samples$Date <- as.Date(with(samples, paste(Year, Month, Day, sep="-")), "%Y-%m-%d")

#fix dates that didn't format correctly
samples$Date[47] <- "2008-06-01"
samples$Date[48] <- "2007-06-02"
samples$Date[172] <- "2016-07-19"
samples$Date[173] <- "2016-07-20"
samples$Date[174] <- "2016-07-25"
samples$Date[175] <- "2016-07-27"
samples$Date[176] <- "2016-08-01"
samples$Date[177] <- "2016-08-09"
samples$Date[178] <- "2016-08-18"
samples$Date[179] <- "2016-08-23"
samples$Date[180] <- "2016-09-03"
samples$Date[181] <- "2016-09-15"
samples$Date[182] <- "2016-09-22"
samples$Date[183] <- "2016-10-25"


#change month from numeric to title
samples$Month <- gsub("03", "March", samples$Month)
samples$Month <- gsub("04", "April", samples$Month)
samples$Month <- gsub("05", "May", samples$Month)
samples$Month <- gsub("06", "June", samples$Month)
samples$Month <- gsub("07", "July", samples$Month)
samples$Month <- gsub("08", "August", samples$Month)
samples$Month <- gsub("09", "September", samples$Month)
samples$Month <- gsub("10", "October", samples$Month)

write.csv(samples, "metadata2.csv")

################# Using new generated meta-data2 table
weather_2006 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2006_P1D.csv")
weather_2007 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2007_P1D.csv")
weather_2008 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2008_P1D.csv")
weather_2009 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2009_P1D.csv")
weather_2010 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2010_P1D.csv")
weather_2011 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2011_P1D.csv")
weather_2012 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2012_P1D.csv")
weather_2013 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2013_P1D.csv")
weather_2014 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2014_P1D.csv")
weather_2015 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2015_P1D.csv")
weather_2016 <- read.csv("data/fr_climat_quotidiennes_QC_7022579_2016_P1D.csv")

#merge all weather dataframes
weather <- Reduce(function(x,y) merge(x, y, all=TRUE), list(weather_2006, weather_2007, weather_2008,
                                                            weather_2009, weather_2010, weather_2011,
                                                            weather_2012, weather_2013, weather_2014,
                                                            weather_2015, weather_2016))

#format as.Date
weather$Date.Heure <- as.Date(weather$Date.Heure)

#rename cols & keep only what is useful to me
weather <- weather %>% rename(Date = Date.Heure, Temp.moy = Temp.moy...C., Precip.tot = Précip..tot...mm.) %>%
                        select(Date, Temp.moy, Precip.tot) 

write.csv(weather, "weather.csv")

#rename col name
new_samples <- samples %>% select(Date)
new_weather <- weather 

temperature_selection <- inner_join(new_samples, new_weather, by = "Date")

#get mean temp of 7 days leading to date
library(lubridate)

#function to select temp from range of 7 days leading up to date x
get_temprange <- function(x){
  weather[weather$Date >= as.Date(x) - 7 & weather$Date <= as.Date(x),]$Temp.moy
}
get_temprange(weather$Date[200])

#function to get avg temp
getmeantemp <- function(x) {
  mean(c(get_temprange(weather$Date[x])))
}
getmeantemp(weather$Date[110])

mean(c(get_temprange(weather$Date[110])))
mean(c(weather$Temp.moy[104:105]))
mean(weather$Temp.moy[104:107])
mean(c(2,5,9,8))

weather[104:105,]

for(x in samples$Date) {
  samples[samples$Date==x,"Temp.moy"] <- func.widthmeans(prefix=x,target.df="weather")
}
rm(x)
df1
weather$Date
weather$Precip.tot
weather$Temp.moy[1]
