#https://blogs.oracle.com/datascience/introduction-to-forecasting-with-arima-in-r

library('ggplot2')
library('forecast')
library('tseries')

#example using: number of bicycles checkouts from a bike sharing service
daily_data = read.csv('/Users/kiristern/Downloads/Bike-Sharing-Dataset/day.csv', header=TRUE, stringsAsFactors=FALSE)

daily_data$Date = as.Date(daily_data$dteday)

#plot the series and visually examine it for any outliers, volatility, or irregularities
ggplot(daily_data, aes(Date, cnt)) + geom_line() + scale_x_date('month')  + ylab("Daily Bike Checkouts") +
  xlab("")
  #bicycle checkouts are showing a lot of fluctuations from one day to another. However, even with this volatility present, we already see some patterns emerge. 
  #For example, lower usage of bicycles occurs in the winter months and higher checkout numbers are observed in the summer months


