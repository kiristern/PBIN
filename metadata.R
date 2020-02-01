library("dplyr")
library("tidyverse")
library(tidyr)
library(lubridate)
library(naniar)

samples <- read.table("/Users/kiristern/Desktop/Shapiro_lab/data/samples.txt")

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

# write.csv(samples, "metadata2.csv")

################# Using new generated meta-data2 table
weather_2006 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2006_P1D.csv")
weather_2007 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2007_P1D.csv")
weather_2008 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2008_P1D.csv")
weather_2009 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2009_P1D.csv")
weather_2010 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2010_P1D.csv")
weather_2011 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2011_P1D.csv")
weather_2012 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2012_P1D.csv")
weather_2013 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2013_P1D.csv")
weather_2014 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2014_P1D.csv")
weather_2015 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2015_P1D.csv")
weather_2016 <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/data/fr_climat_quotidiennes_QC_7022579_2016_P1D.csv")

#merge all weather dataframes
weather <- Reduce(function(x,y) merge(x, y, all=TRUE), list(weather_2006, weather_2007, weather_2008,
                                                            weather_2009, weather_2010, weather_2011,
                                                            weather_2012, weather_2013, weather_2014,
                                                            weather_2015, weather_2016))
#replace missing values with NA
weather <-replace_na(weather)

#format as.Date
weather$Date.Heure <- as.Date(weather$Date.Heure)

#rename cols & keep only what is useful to me
weather <- weather %>% rename(Date = Date.Heure, Temp.moy = Temp.moy...C., Precip.tot = Précip..tot...mm.) %>%
                        select(Date, Temp.moy, Precip.tot) 

#convert commas to periods & keep values as.numeric
weather$Temp.moy <- as.numeric(gsub(",", ".", weather$Temp.moy))

#convert commas to periods & keep values as.numeric
weather$Precip.tot <- as.numeric(gsub(",", ".", weather$Precip.tot))

# write.csv(weather, "weather.csv")

##get mean temp of 7 days leading to date
#function to select range of t-7:t
get_date_range <- function(x){
  weather[weather$Date >= as.Date(x) - 7 & weather$Date <= as.Date(x),]
}
get_date_range(samples$Date[100])

#function to get mean temp from t-7:t
get_mean_temp <- function(x){
  y = get_date_range(x)
  return(mean(y$Temp.moy))
}

#function to get mean precipitation from t-7:t
get_mean_prec <- function(x){
  y = get_date_range(x)
  return(mean(y$Precip.tot))
}

#get mean temp for each sample
for (i in 1:length(samples$Date)){
  samples$Mean_temperature_t0_t7[i] <- get_mean_temp(samples$Date[i])
}

#get mean prec for each sample
for (i in 1:length(samples$Date)){
  samples$Cumulative_precipitation_t1_t7_mm[i] <- get_mean_prec(samples$Date[i])
}


##########
meta <- read.csv("/Users/kiristern/Desktop/Shapiro_lab/mapping_bloom2_new_corrected2.csv")
meta2 <- meta
#View(meta2)

#replace missing values with NA
meta2 <- meta2 %>% mutate_all(~replace(., . ==0, NA))

#rename
meta2<- rename(meta2, "Date" = "Sample")

#remove everything before 4th period in sampleID (just to keep date)
meta2$Date <- gsub("*........(.*)\\-.*","\\1", meta2$SampleID)

#extract proper year
meta2$Years <- gsub(".*-(.*)\\-.*", "\\1", meta2$SampleID, perl=T)

#extract proper month
#meta2$Months <- gsub(".*-.(.*)\\-.*-.*", "\\1", meta2$SampleID, perl=T) #"." before ( indicates to remove the first character between xx-Xx-xxxx
meta2$Months <- gsub(".*-(.*)\\-.*-.*", "\\1", meta2$SampleID, perl=T)

#change month from numeric to title
meta2$Months <- gsub("03", "March", meta2$Months)
meta2$Months <- gsub("04", "April", meta2$Months)
meta2$Months <- gsub("05", "May", meta2$Months)
meta2$Months <- gsub("06", "June", meta2$Months)
meta2$Months <- gsub("07", "July", meta2$Months)
meta2$Months <- gsub("08", "August", meta2$Months)
meta2$Months <- gsub("09", "September", meta2$Months)
meta2$Months <- gsub("10", "October", meta2$Months)

#assign season depending on date
getSeason <- function(DATES) {
  WS <- as.Date("15-12-2012", format = "%d-%m-%Y") # Winter Solstice
  SE <- as.Date("15-3-2012",  format = "%d-%m-%Y") # Spring Equinox
  SS <- as.Date("15-6-2012",  format = "%d-%m-%Y") # Summer Solstice
  FE <- as.Date("15-9-2012",  format = "%d-%m-%Y") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Fall")))
}
#format as.Date
meta2$Date <- as.Date(meta2$Date, "%d-%m-%Y")
#assign season depending on date
meta2$Period <- getSeason(meta2$Date)



