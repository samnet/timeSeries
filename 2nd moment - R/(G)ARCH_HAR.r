
### Set Up ####

## libraries

library(xts)
library(tseries)
library(highfrequency)
library(zoo)
library(its)
library(car)
library(FinTS)
library(fGarch)
library(highfrequency)
library(timeDate)

library(stargazer)
#################################################

## Data cleaning ##

rearrangeIt <- function(data, newname) { # this function:
    # 1). selects appropriate column, rearrange them,
    # 2). convert the file to the xts format
    # 3). does the house cleaning
    # running time: app. 2 minutes
    raw <- read.csv(file=data, header = T, sep = ',')
    raw <- subset(raw, select=c(2:4))
    write.csv(raw, file = newname, row.names=FALSE)
    newname <- read.csv(file=newname, sep =',', header=T)
    dayTime <- timeDate(paste(newname[,1], newname[,2]), format = "%m/%d/%Y %H:%M:%S")
    newname <- xts(newname[,3],dayTime)
    names(newname)<-'price'
    return(newname)
}

rearrangeIt('SP.csv', 'SPProcessed.csv') -> prices


##################################################



## Signature Plot ##

# Parameters #

rangeFreq <- 30

# Implementation #

indMat <- function(mat, rge){
    rawMat<- matrix(0,length(mat),rge)
    for (i in 1:(length(mat))){
        for (j in (1:rge)){
            freq <-j*4
            if ((i%%freq)==0){
                rawMat[i,j] <- mat[i]
            }
        }
        if (i%%10000==0){print(i)} # to see whether it is alive

    }
    return(rawMat)
}

gridPrice <- indMat(prices, rangeFreq)

gridToList<- function(mat){
    raw <- mat
    out <- list()
    for (i in 1:(dim(mat)[2])){
        out[[i]]<-raw[,i]
        print(i)
    }
    return(out)
}

gridList <- gridToList(gridPrice)

NzList <-function(lst){
    raw<-list()
    for (i in 1:length(lst)){
        tempV <- (lst[[i]])
        raw[[i]] <- tempV[tempV>0]
    }
    return(raw)
}

gridListNz <- NzList(gridList)

priceToLogM <- function(mat){
    raw <- diff(log(mat))
    raw <- raw[-1,] # drops first value (NA)
    return(raw)
    # lag's default behavior matches the common time-series interpretation of that operator --- specifically that a value at time 't' should be the value at time 't-1' for a positive lag.
}

priceToLogV <- function(v){
    raw <- diff(log(v))
    raw <- raw[-1] # drops first value (NA)
    return(raw)
    # lag's default behavior matches the common time-series interpretation of that operator --- specifically that a value at time 't' should be the value at time 't-1' for a positive lag.
}

gridRts <- lapply(gridListNz, priceToLogV)

sumSq <- function(v){
    raw <- v
    raw <- raw^2
    out <- sum(raw)
    return(out)
}

sampleVar <- lapply(gridRts, sumSq)
png('volSign.png')
plot(1:rangeFreq, sampleVar, type='lines', xlab='Interval between measurements (minutes)', ylab='Sample Variance', main='sample volatility as a function of measurement interval')
dev.off()

###########################################

## RV ##

# Parameters / constants #

dayLenM <- 6.5*60 # in minutes
optInt <- 5 # in minutes
obsPerDay <- dayLenM/optInt

optSubs <-gridRts[[optInt]]
optSubsSq <- (gridRts[[optInt]])^2
obsNum <- length(optSubsSq)
days <- round(obsNum/obsPerDay)

# Implementation #

if (length(prices)-length(optSubsSq)*(optInt*4) > 25){
    print('PROBLEM HERE')
}

dayRts <- function(data, periods, obserPerPeriod){
    rt <- c()
    for (i in 1:periods){
        dayObs <- data[((i-1)*obserPerPeriod+1):(i*obserPerPeriod)]
        rt[i] <- sum(dayObs)
    }
    return(rt)
}

dailyRts <- dayRts(optSubs, days, obsPerDay)

dayRV <- function(data, periods, obserPerPeriod){
    rv <- c()
    for (i in 1:periods){
        dayObs <- data[((i-1)*obserPerPeriod+1):(i*obserPerPeriod)]
        rv[i] <- sum(dayObs)
    }
    if (length(rv) != periods){
        print('PROBLEM HERE!!!')
    }
    return(rv)
}

dailyRV <- dayRV(optSubsSq, days, obsPerDay)

png('RealizedVol.png')
plot(1:days, dailyRV, type='line', main= 'Daily Realized Volatility', xlab = 'days (starting 02.01.2004)', ylab = '')
dev.off()
densRV <- density(dailyRV)

png('Density.png')
plot(densRV, col = "Red", main='Density of the Realized Volatility', xlab='Value')
dev.off()

###########################################

## HAR ##

# Parameters / Regressors #

windL <- 21
windM <- 5
windS <- 1

RVW = NULL
RVM = NULL
RVYes = NULL
signYes = NULL

trainRV <- dailyRV[(windL+1):length(dailyRV)]
trainRts <- dailyRts[(windL+1):length(dailyRV)]

seriesLength <- length(trainRV)

movingAv <- function(lags, series, len){
    inp <- series
    indb <- windL -lags
    inde <- windL-1
    raw  <- c()
    for (i in 1:len){
        sum <-sum(series[(i+indb):(i+inde)])/lags
        raw[i] <- sum
    }
    return(raw)
}

signY <- function(series, len){
    si  <- c()
    for (i in 1:len){
        val <-series[(windL-1+i)]
        si[i] <-sign(val)
    }
    return(si)
}

# Implementation #

RVW <- movingAv(5, dailyRV, seriesLength)
RVM <- movingAv(21, dailyRV, seriesLength)
RVYes <-movingAv(1, dailyRV, seriesLength)
signYes = signY(dailyRts, seriesLength)

length(RVM)
length(RVW)
length(signYes)
length(trainRV)

trainset <-list(trainRV, RVYes, RVW, RVM, signYes)
trainset <-as.data.frame(trainset)
attach(trainset)
HAR = lm(trainRV~RVYes+RVW+RVM, data=trainset); summary(HAR)
HARAsy <-lm(trainRV~RVYes+RVW+RVM+signYes, data=trainset); summary(HARAsy)
HARPred = predict(HAR)
png('QQplot.png')
plot(HARPred, trainRV, xlab= 'HAR', ylab = 'Training Data', main='QQ-Plot HAR vs. Training Data')
dev.off()
png('RVHAR.png')
plot(1:length(HARPred), trainRV, type= 'lines', main= 'RV (black) and HAR (red)', xlab = 'days (starting 02.01.2004)', ylab = '')
points(HARPred, type='lines', col="red")
dev.off()
