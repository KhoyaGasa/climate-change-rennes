#Setting working directory
setwd("...../Econometrie_Appliqué/Devoirs/Evolution_de_la_temp")
##Loading necessary libraries
library(dplyr)
library(questionr)
library(ade4)
library(dynlm)
library(sandwich)
library(ggplot2)
library(tidyverse)
library(gtsummary)


## Importing the dataset
temps <- read.csv("rennes_temp.csv", sep = ";", header = T)
temp <- temps[,c(1,3,4)]
temp$AVG <- ((temp$TMAX..Degrees.Fahrenheit.+temp$TMIN..Degrees.Fahrenheit.)/2)
temp_clean = `colnames<-`(temp, c("Date","tmax","tmin","AVG"))

## writing a new clean csv file for use
#write.csv(temp_clean, file="temperatures_rennes.csv")

## converting from Fahrenheit to celsius
conv_celsius <- function (a){
  temp_cel = round(((a - 32)* (5/9)),2)
  temp_cel
  
}
temp_clean[,c(2,3,4)] <- lapply(temp_clean[,c(2,3,4)], conv_celsius) # applying the function to all the columns
temp_clean$DTR = temp_clean$tmax - temp_clean$tmin ## Creation of the Daily temperature range(DTR) column
## Setting the date sequence
inds <- seq.Date(from =as.Date("1945-01-01"), to = as.Date("2022-10-27"), by= "day")
##Creating the time series object for the average temperature
AVG <- ts(temp_clean$AVG, start = c(1945,as.numeric(format(inds[1],"%j"))),
          frequency = 365)

##Creating the linear trend
time.pts = c(1:length(AVG))

reg1 <- lm(AVG~time.pts)
summary(reg1)
coeftest(reg1, vcov = NeweyWest(reg1,  lag = 18,
                                prewhite = F)) ## Correction of autocorrelation and for heteroskedasticity
est_values <- fitted(reg1)

## creation of time series object for the fitted values
est_values_reg1 <- ts(est_values, start = c(1945,as.numeric(format(inds[1],"%j"))),
                      frequency = 365)
plot(est_values_reg1) #>= 1.5 degrees celsius

plot(AVG,ylab="", col="gray",sub="Rennes",xlim =c(1945,2022), main="Evolution de la temperature moyenne")
lines(est_values_reg1,lwd=2,col="blue")

## Analysis of the DTR
reg2 <- lm(DTR ~time.pts)
summary(reg2)
coeftest(reg2, vcov = NeweyWest(reg1,  lag = 18,
                                prewhite = F))## Correction of autocorrelation and for heteroskedasticity
                                ## Trend is not signicative in the results, why ??
est_vals_dtr <- fitted(reg2)
plot(est_vals_dtr) ##only 0.25 degrees
## creation of time series object for the fitted values of DTR
est_vals_dtr_reg <- ts(est_vals_dtr, start = c(1945,as.numeric(format(inds[1],"%j"))),
                       frequency = 365)
plot(DTR,ylab="", col="gray",sub="Rennes",xlim =c(1945,2022), main="Evolution de l'�tendue de temp�rature `\n` (DTR) ")
lines(est_vals_dtr,lwd=2,col="blue") ## Warming happening more reliably, gap between max and min dropping.
                                    ## Min is increasing faster than max temperatures.

## Seasonality analysis
##Construction of of the months binary dummies
temp_clean$D1 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "01", 1,0)
temp_clean$D2 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "02", 1,0)
temp_clean$D3 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "03", 1,0)
temp_clean$D4 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "04", 1,0)
temp_clean$D5 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "05", 1,0)
temp_clean$D6 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "06", 1,0)
temp_clean$D7 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "07", 1,0)
temp_clean$D8 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "08", 1,0)
temp_clean$D9 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "09", 1,0)
temp_clean$D10 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "10", 1,0)
temp_clean$D11 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "11", 1,0)
temp_clean$D12 <- ifelse(format((as.Date(temp_clean$Date, format = "%d/%m/%Y")), "%m")== "12", 1,0)


## Extraction of the seasonality effects in the series of AVG
dat=decompose(na.approx(AVG)) #replacing NAs by means of interpolation
plot(dat)
plot(dat$trend)
detrended_data <- na.approx(AVG) - dat$trend ## Detrending the time series

## Inclusion de saisonalité avec les dummy
mod6 <- dynlm(detrended_data ~ -1 + temp_clean$D1+temp_clean$D2+temp_clean$D3+
                temp_clean$D4+
                temp_clean$D5+temp_clean$D6+temp_clean$D7+temp_clean$D8 +
                temp_clean$D9+temp_clean$D10+temp_clean$D11+temp_clean$D12)
summary(mod6)

temp.fit.lm <- ts(fitted(mod6),     # random data
                   start = c(1945, as.numeric(format(inds[1], "%j"))),
                   frequency = 365) ## time series object of the estimated values with months in account
plot(detrended_data,ylab="",sub="Rennes", col="gray")
lines(temp.fit.lm,lwd=2,col="blue")
st1 = coef(mod6) # Extraction of the coefficients of each month

plot(1:12,st1,lwd=2,type="s",main="AVG",sub="Rennes",ylab="", xaxt = "n",
     col="blue",xlab='Time')
axis(1, at=1:12,
     labels=c("J","F","M","A","M","J","J","A",
              "S","O","N","D"))

## Analysis of seasonality effects on the DTR
dtr_data = decompose(na.approx(DTR)) # Breaking down the time series object and replacing any NAs
plot(dtr_data)
stationary_dtr <- na.approx(DTR) - dtr_data$trend ## Detrended/Stationary DTR

## regression analysis of the stationary dtr data
dtr.mod <-dynlm(stationary_dtr ~ -1 + temp_clean$D1+temp_clean$D2+temp_clean$D3+
                  temp_clean$D4+
                  temp_clean$D5+temp_clean$D6+temp_clean$D7+temp_clean$D8 +
                  temp_clean$D9+temp_clean$D10+temp_clean$D11+temp_clean$D12)
summary(dtr.mod)
dtr.fit.lm <- ts(fitted(dtr.mod),     # random data
                 start = c(1945, as.numeric(format(inds[1], "%j"))),
                 frequency = 365) ## time series object of the estimated values with months in account
plot(stationary_dtr, ylab="", sub= "Rennes", col= "gray")
lines(dtr.fit.lm, lwd=2, col="red") ## graphs of DTR against the estimated effects of the months
st2 = coef(dtr.mod)
plot(1:12,st2,lwd=2,type="s",main="DTR",sub="Rennes",ylab="", xaxt = "n",
     col="red",xlab='Time')
axis(1, at=1:12,
     labels=c("J","F","M","A","M","J","J","A",
              "S","O","N","D"))

