library(ncdf4)
library(tidyr)
library(zoo)
library(lubridate)
library(forecast)
library(prophet)
library(sqldf)

# read_cru_hemi.r
#
# Reads a CRU-format hemispheric average file, as provided at
# http://www.cru.uea.ac.uk/cru/data/temperature
#
# Format has two lines for each year.
#  1) monthly mean anomalies plus an annual mean
#  2) coverage percentages
#
# Returns a data frame with columns:
#  year (1850 to final year)
#  annual (mean annual anomaly)
#  month.1 ... month.12 (mean monthly anomaly)
#  cover.1 ... cover.12 (percentage coverage)
#
read_cru_hemi <- function(filename) {
  # read in whole file as table
  tab <- read.table(filename,fill=TRUE)
  nrows <- nrow(tab)
  # create frame
  hemi <- data.frame(
    year=tab[seq(1,nrows,2),1],
    annual=tab[seq(1,nrows,2),14],
    month=array(tab[seq(1,nrows,2),2:13]),
    cover=array(tab[seq(2,nrows,2),2:13])
  )
  # mask out months with 0 coverage
  hemi$month.1 [which(hemi$cover.1 ==0)] <- NA
  hemi$month.2 [which(hemi$cover.2 ==0)] <- NA
  hemi$month.3 [which(hemi$cover.3 ==0)] <- NA
  hemi$month.4 [which(hemi$cover.4 ==0)] <- NA
  hemi$month.5 [which(hemi$cover.5 ==0)] <- NA
  hemi$month.6 [which(hemi$cover.6 ==0)] <- NA
  hemi$month.7 [which(hemi$cover.7 ==0)] <- NA
  hemi$month.8 [which(hemi$cover.8 ==0)] <- NA
  hemi$month.9 [which(hemi$cover.9 ==0)] <- NA
  hemi$month.10[which(hemi$cover.10==0)] <- NA
  hemi$month.11[which(hemi$cover.11==0)] <- NA
  hemi$month.12[which(hemi$cover.12==0)] <- NA
  #
  return(hemi)
}


#####____________________________Prepare MetOffice Data__________________####################3

MetOffice_Data <- read_cru_hemi("C:\\Users\\panes\\Desktop\\Jaspal\\Queens MMA\\MMA 867 - Predictive Modelling\\A2\\HadCRUT4-gl.dat")
col <- c(1,3:14)
MetOffice_Data <- MetOffice_Data[,col]
MetOffice_Data <- gather(MetOffice_Data,Month,AnomalyTemp, month.1:month.12)
MetOffice_Data$Month <- as.integer(gsub('month.', '', MetOffice_Data$Month))
#MetOffice_Data <- mutate(MetOffice_Data, YearMonth = paste(MetOffice_Data$year, "-", MetOffice_Data$Month)) 
#MetOffice_Data$YearMonth <- as.yearmon(MetOffice_Data$YearMonth)
MetOffice_Data <- mutate(MetOffice_Data, YearMonth= make_datetime(year, Month, 1) ) 
MetOffice_Data$YearMonth <- as.Date(MetOffice_Data$YearMonth)
#MetOffice_Data <- MetOffice_Data[order(MetOffice_Data$YearMonth)]
MetOffice_Data <- arrange(MetOffice_Data, YearMonth)
MetOffice_Data <- MetOffice_Data[!is.na(MetOffice_Data$AnomalyTemp),]
MetOffice_Data <- mutate(MetOffice_Data, AnomalyTempC= AnomalyTemp+14 ) 
MetOffice_Data <- mutate(MetOffice_Data, AnomalyTempF= AnomalyTempC ) 
MetOffice_Data <- mutate(MetOffice_Data, AnomalyTempF= ((AnomalyTempC*9)/5)+32) 

#####____________________________Prepare Nasa Data__________________####################3


NASA_data<-read.csv("C:\\Users\\panes\\Desktop\\Jaspal\\Queens MMA\\MMA 867 - Predictive Modelling\\A2\\GLB.Ts+dSST.csv", header=TRUE, sep=",", skip=1)
NASA_data<- NASA_data[,1:13]
NASA_data <- gather(NASA_data,Month,AnomalyTemp, Jan:Dec)
NASA_data$Month <- match(NASA_data$Month, month.abb)
NASA_data$AnomalyTemp <-as.numeric(NASA_data$AnomalyTemp)
#NASA_data <- mutate(NASA_data, YearMonth = paste(NASA_data$Year, "-",NASA_data$Month))
#NASA_data$YearMonth <- as.yearmon(NASA_data$YearMonth)
NASA_data <- mutate(NASA_data, YearMonth= make_datetime(Year, Month, 1) ) 
NASA_data$YearMonth <- as.Date(NASA_data$YearMonth)
#NASA_data <- NASA_data[order(NASA_data$YearMonth)]
NASA_data <- arrange(NASA_data, YearMonth)
NASA_data <- NASA_data[!is.na(NASA_data$AnomalyTemp),]
NASA_data <- mutate(NASA_data, AnomalyTempC= AnomalyTemp+14 ) 
NASA_data <- mutate(NASA_data, AnomalyTempF= ((AnomalyTempC*9)/5)+32) 



#########_____________________Explore Data____________________________________________________
str(NASA_data)
str(MetOffice_Data)


## subset data
#MetOffice_Data <- filter(MetOffice_Data, MetOffice_Data$year>=1950)
#NASA_data <- filter(MetOffice_Data, MetOffice_Data$year>=1950)

MetOffice_DataFB <- sqldf('SELECT YearMonth, AnomalyTempC FROM MetOffice_Data')
names(MetOffice_DataFB)[1] = "ds"
names(MetOffice_DataFB)[2] = "y"
NASA_dataFB <- sqldf('SELECT YearMonth, AnomalyTempC FROM NASA_data')
names(NASA_dataFB)[1] = "ds"
names(NASA_dataFB)[2] = "y"
#####____________________________Take Subsets & Create various TS__________________####################3


NASA_data_ts <- ts(NASA_data$AnomalyTempF,start=1880, frequency=12) # ts function defines the dataset as timeseries starting Jan 2004 and having seasonality of frequency 12 (monthly) 

MetOffice_data_ts <- ts(MetOffice_Data$AnomalyTempF,start=1850, frequency=12) # ts function defines the dataset as timeseries starting Jan 2004 and having seasonality of frequency 12 (monthly) 

#NASA_data_ts <- ts(NASA_data$AnomalyTempC,start=1950, frequency=12) # ts function defines the dataset as timeseries starting Jan 2004 and having seasonality of frequency 12 (monthly) 

#MetOffice_data_ts <- ts(MetOffice_Data$AnomalyTempC,start=1950, frequency=12) # ts function defines the dataset as timeseries starting Jan 2004 and having seasonality of frequency 12 (monthly) 


plot(NASA_data_ts)
plot(MetOffice_data_ts)

fit <- decompose(NASA_data_ts, type="multiplicative") #decompose using "classical" method, multiplicative form
plot(fit)

fit <- decompose(MetOffice_data_ts, type="multiplicative") #decompose using "classical" method, multiplicative form
plot(fit)




#######################___________Modeling______________________


NAS_AAN <- ets(NASA_data_ts, model="AAN")
NAS_AAZ <- ets(NASA_data_ts, model="AAZ",damped=FALSE)
NAS_MMN <- ets(NASA_data_ts, model="MMN",damped=FALSE)
NAS_MMZ <- ets(NASA_data_ts, model="MMZ",damped=FALSE)

# Create their prediction "cones" for 360 months (30 years) into the future with quintile confidence intervals
NASA_data_AAN_pred <- forecast(NAS_AAN, h=960, level=c(0.8, 0.95))
NASA_data_AAZ_pred <- forecast(NAS_AAZ, h=960, level=c(0.8, 0.95))
NASA_data_MMN_pred <- forecast(NAS_MMN, h=960, level=c(0.8, 0.95))
NASA_data_MMZ_pred <- forecast(NAS_MMZ, h=960, level=c(0.8, 0.95))


par(mfrow=c(1,4)) # This command sets the plot window to show 1 row of 4 plots
plot(NASA_data_AAN_pred, xlab="Year", ylab="Predicted Temperature")
plot(NASA_data_AAZ_pred, xlab="Year", ylab="Predicted Temperature")
plot(NASA_data_MMN_pred, xlab="Year", ylab="Predicted Temperature")
plot(NASA_data_MMZ_pred, xlab="Year", ylab="Predicted Temperature")





#Create a trigonometric box-cox autoregressive trend seasonality (TBATS) model
NASA_data_tbats <- tbats(NASA_data_ts)
NASA_tbats_pred <-forecast(NASA_data_tbats, h=960, level=c(0.8, 0.95))
par(mfrow=c(1,1))
plot(NASA_tbats_pred, xlab="Year", ylab="Predicted Temperature")



###__________auto-arima___________________


fit <- auto.arima(NASA_data_ts,seasonal=TRUE)
fit

par(mfrow=c(1,1))
Acf(residuals(fit))
plot(forecast(fit,960)) #60 stands for 60 months = 5 years






# Met_Data

MET_AAN <- ets(MetOffice_data_ts, model="AAN")
MET_AAZ <- ets(MetOffice_data_ts, model="AAZ",damped=FALSE)
MET_MMN <- ets(MetOffice_data_ts, model="MMN",damped=FALSE)
MET_MMZ <- ets(MetOffice_data_ts, model="MMZ",damped=FALSE)

# Create their prediction "cones" for 360 months (30 years) into the future with quintile confidence intervals
MET_data_AAN_pred <- forecast(MET_AAN, h=960, level=c(0.8, 0.95))
MET_data_AAZ_pred <- forecast(MET_AAZ, h=960, level=c(0.8, 0.95))
MET_data_MMN_pred <- forecast(MET_MMN, h=960, level=c(0.8, 0.95))
MET_data_MMZ_pred <- forecast(MET_MMZ, h=960, level=c(0.8, 0.95))


par(mfrow=c(1,4)) # This command sets the plot window to show 1 row of 4 plots
plot(MET_data_AAN_pred, xlab="Year", ylab="Predicted Temperature")
plot(MET_data_AAZ_pred, xlab="Year", ylab="Predicted Temperature")
plot(MET_data_MMN_pred, xlab="Year", ylab="Predicted Temperature")
plot(MET_data_MMZ_pred, xlab="Year", ylab="Predicted Temperature")





#Create a trigonometric box-cox autoregressive trend seasonality (TBATS) model
MET_data_tbats <- tbats(MetOffice_data_ts)
MET_tbats_pred <-forecast(MET_data_tbats, h=960, level=c(0.8, 0.95))
par(mfrow=c(1,1))
plot(MET_tbats_pred, xlab="Year", ylab="Predicted Temperature")



###__________auto-arima___________________

fit <- auto.arima(MetOffice_data_ts,seasonal=TRUE)
fit

par(mfrow=c(1,1))
Acf(residuals(fit))
plot(forecast(fit,960)) #60 stands for 60 months = 5 years






#### fb

m <- prophet(MetOffice_DataFB)
future <- make_future_dataframe(m, periods = 960)
forecast <- predict(m, future)

plot(m,forecast)