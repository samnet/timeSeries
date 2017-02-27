#delete data
rm(list = ls(all = TRUE))

#Install libraries

#install.packages("xts", dependencies = TRUE)
#install.packages("tseries", dependencies = TRUE)
#install.packages("highfrequency", dependencies = TRUE)
#install.packages("zoo", dependencies = TRUE)
#install.packages("its", dependencies = TRUE)
#install.packages("car", dependencies = TRUE)
#install.packages("FinTS", dependencies = TRUE)
#install.packages("fGarch", dependencies = TRUE)
#install.packages("rugarch", dependencies = TRUE)

# load libraries

library(xts)
library(tseries)
library(highfrequency)
library(zoo)
library(its)
library(car)
library(FinTS)
library(fGarch)
#library(rugarch)

# Data handling

# We create a new csv that contains solely the needed data for our analysis, change the path to your SP file.
#path = "SP.csv" #change the path to your SP file.
#SP = read.csv(file=path, , sep = ",", header=T)
#Delete first, second and last row and write a new file
#SP[5] <- NULL
#SP[1] <- NULL
#write.csv(SP, file = "SP2.csv", row.names=FALSE) #change the path to store new SP file.

#read data from csv
path = "SP2.csv" #change the path to your new SP file
SP2 = read.csv(file=path, , sep = ",", header=T)

#Prices as extended time series (xts) - full
timedate <- timeDate(paste(SP2[,1], SP2[,2]), format = "%m/%d/%Y %H:%M:%S")
sp.c     <- xts(SP2[,3],timedate)

#Prices as extended time series (xts) - 1 day
day = '2004-01-05'
sp_day.c <- xts(subset(SP2[,3],substr(timedate,0,10)==day),subset(timedate,substr(timedate,0,10)==day))

#compute as log returns, drop the first value (N/A)
sp.r   <- diff(log(sp.c))
sp_day.r   <- diff(log(sp_day.c))
sp.r <- sp.r[-1,]
sp_day.r  <- sp_day.r [-1,]

# Descriptive statistics

#plot prices
#par(mfrow=c(2,2))
plot(sp.c, main="Standard & Poor's 500 index", ylab="Price")
#intraday plot
plot(sp_day.c , main="Standard & Poor's 500 index 1 day", ylab="Price")

#plot returns
plot(sp.r*100, main="Standard & Poor's 500 index", ylab="Returns in Percent")
#intraday plot
plot(sp_day.r*100, main="Standard & Poor's 500 index 1 day", ylab="Returns in Percent")

#some basic statistics
basicStats(sp_day.r)

#Compare normal distribution function to SP distribution.
#par(mfrow=c(2,2))
a=density(sp.r)
plot(a,main="Density SP vs Normal")
lines(a$x,dnorm(a$x,mean=mean(sp.r, na.rm=T), sd=sqrt(var(sp.r, na.rm=T))),lty=2, col="red")
chart.QQPlot(sp.r, main = "qq-Plot SP vs Normal", distribution = 'norm', col=c("black","red"))

#intraday version
a=density(sp_day.r)
plot(a,main="Intraday density SP vs Normal")
lines(a$x,dnorm(a$x,mean=mean(sp_day.r, na.rm=T), sd=sqrt(var(sp_day.r, na.rm=T))),lty=2, col="red")
chart.QQPlot(sp_day.r, main = "Intraday qq-Plot SP vs Normal", distribution = 'norm', col=c("black","red"))

#Jarque-Bera test for normality
jarque.bera.test((sp.r))
jarque.bera.test((sp_day.r))

#For the acf and the pacf fct. we need equidistant ticks, this function will use the last realized value for the value of the grid point.
aggmin = aggregatets(sp_day.r,on="minutes",k=1)

#autocorrelation function intraday
#par(mfrow=c(2,2))
acf(aggmin,main="SP intraday log returns")
pacf(aggmin,main="SP intraday log returns")
acf(aggmin^2,main="SP intraday squared returns")
pacf(aggmin^2,main="SP intraday squared returns")

#****************GARCH Analysis-closing prices********************

#We get the closing prices for each day.

sp_sum.c = to.period(sp.c, 'days', drop.time = TRUE)
index(sp_sum.c) = trunc(index(sp_sum.c),"days")

#We get the log returns for each day.

sp_sum.r = diff(log(sp_sum.c[,4]))
sp_sum.r <- sp_sum.r[-1,]

# Descriptive statistics

#plot prices
#par(mfrow=c(2,2))
plot(sp_sum.c, main="Standard & Poor's 500 index closing prices", ylab="Price")

#plot returns
plot(sp_sum.r*100, main="Standard & Poor's 500 index closing prices", ylab="Returns in Percent")

#some basic statistics
basicStats(sp_sum.r)

#Jarque-Bera test for normality
jarque.bera.test((sp_sum.r))

#Compare normal distribution function to SP distribution.
par(mfrow=c(2,2))
a=density(sp_sum.r)
plot(a,main="Density SP vs Normal")
lines(a$x,dnorm(a$x,mean=mean(sp_sum.r, na.rm=T), sd=sqrt(var(sp_sum.r, na.rm=T))),lty=2, col="red")
chart.QQPlot(sp_sum.r, main = "qq-Plot SP vs Normal", distribution = 'norm', col=c("black","red"))

##Function from R-Code of the lecture "Financial Volatility" to calculate the 95% GARCH confidence bounds

#function for gamma
gamma=function(x,h)
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  gamma=sum(x[1:(n-h)]*x[(h+1):n])/n
}

#function for roh
rho=function(x,h)
{
  rho=gamma(x,h)/gamma(x,0)
}

#function for the new acf plot
n1.acf=function(x,main=NULL,method="NP")
{
  x<-as.numeric(x)
  n=length(x)
  nlag=as.integer(min(10*log10(n),n-1))
  acf.val=sapply(c(1:nlag),function(h) rho(x,h))
  x2=x^2
  var= 1+(sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band=sqrt(var/n)
  minval=1.2*min(acf.val,-1.96*band,-1.96/sqrt(n))
  maxval=1.2*max(acf.val,1.96*band,1.96/sqrt(n))
  acf(x, xlab="Lag", ylab="Sample autocorrelations", main=main,
      ylim=c(minval,maxval))
  lines(c(1:nlag),-1.96*band,lty=1,col="red")
  lines(c(1:nlag),1.96*band,lty=1,col="red")
}

# Now we have a look at the box test results
Box.test(sp_sum.r,lag=5,type="Ljung") #Ljung-Box statistic Q(5) Box.test
Box.test(sp_sum.r,lag=10,type="Ljung") #Ljung-Box statistic Q(10) Box.test
Box.test(sp_sum.r,lag=20,type="Ljung") #Ljung-Box statistic Q(20) Box.test
# - > p-values of about 0 for all the tested lags, we therefore reject the null hypothesis of no autocorrelations at the 1% confidence level.

# Do we have ARCH effects?

#autocorrelation function closing prices for S&P
#par(mfrow=c(2,1))
n1.acf(sp_sum.r,main="SP closing price returns")
pacf(sp_sum.r,main="SP closing price returns")
acf(sp_sum.r^2,main="SP closing price squared returns")
pacf(sp_sum.r^2,main="SP closing price squared returns")


##Function from R-Code of the lecture "Financial Volatility" to calculate the Lagrange Multiplier
LM=function(x,h)
{
  n=length(x)
  x2=x^2-mean(x^2)
  dat<-matrix(,n-h,h+1)
  for (i in 1:(h+1))
  {
    dat[,i]=x2[(h+2-i):(n-i+1)]
  }
  a=lm(dat[,1]~dat[,2:(h+1)])
  r2=summary(a)$r.squared
  print(r2 * n)
  print(1-pchisq(r2*n,h))
}

LM(sp_sum.r,5)
LM(sp_sum.r,10)
LM(sp_sum.r,20)
# - > we reject the null hypothesis that alpha_1=alpha_2=...=alpha_q=0 and therefore assume ARCH effects in the series.

# We have a look again at the box test under the assumption of ARCH effect.
##Function from R-Code of the lecture "Financial Volatility" to calculate the corrected ljung box portmanteau test under GARCH

#function for the asymptotic gamma
gamma.asy=function(x,h)
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  x2=x^2
  gamma.asy<-matrix(,h,h)
  for (i in 1:h)
  {
    for (j in i:h)
    {
      gamma.asy[i,j]=gamma.asy[j,i]=sum(
        x[(j-i+1):(n-i)]*x[1:(n-j)]*x2[(j+1):n])/n
    }
  }
  rho.asy=1/gamma(x,0)^2*gamma.asy
  list(gamma.asy=gamma.asy,rho.asy=rho.asy)
}

#function for the corrected box test
corr.Box.test=function(x,h)
{
  n<-length(x)
  a=gamma.asy(x,h)
  acf.val=sapply(c(1:h),function(h) rho(x,h))
  val=n*(acf.val%*%solve(a$rho.asy)%*%acf.val)
  print(val)
  print(1-pchisq(val,h))
}

# Corrected Ljung Box
corr.Box.test(sp_sum.r,5) #Q(5)
corr.Box.test(sp_sum.r,10) #Q(10)
corr.Box.test(sp_sum.r,20) #Q(20)
# - > again p-values of about 0 for all the tested lags, we therefore reject the null hypothesis of no autocorrelations at the 1% confidence level.

# Analysis using a GARCH(1,1) model (as benchmark)
#sp_sum.d11=garchFit(sp_sum.r~arma(1,1) + garch(1,1),trace=F)
sp_sum.d11=garchFit(sp_sum.r~garch(1,1),trace=F)

# Obtain results
summary(sp_sum.d11)

#residual diagnostic of GARCH(1,1)
par(mfrow=c(1,2))
acf(sp_sum.d11@residuals/sp_sum.d11@sigma.t,main="GARCH(1,1) residuals")
acf(sp_sum.d11@residuals^2/sp_sum.d11@sigma.t^2,main="GARCH(1,1) squared residuals")
# - > We reject the null hypothesis of no independence

#Jarque-Bera test for normality
jarque.bera.test((sp_sum.d11@residuals/sp_sum.d11@sigma.t))
# - > we reject normality assumption for the residuals

#Compare normal distribution function to SP residual distribution.
#par(mfrow=c(2,2))
a=density(sp_sum.d11@residuals/sp_sum.d11@sigma.t)
plot(a,main="GARCH(1,1) SP Residuals vs Normal")
lines(a$x,dnorm(a$x,mean=mean(sp_sum.d11@residuals/sp_sum.d11@sigma.t, na.rm=T), sd=sqrt(var(sp_sum.d11@residuals/sp_sum.d11@sigma.t, na.rm=T))),lty=2, col="red")
chart.QQPlot(sp_sum.d11@residuals/sp_sum.d11@sigma.t, main = "qq-Plot GARCH(1,1) SP Residual vs Normal", distribution = 'norm', col=c("black","red"))

AIC = matrix(nrow=4, ncol=4)
BIC = matrix(nrow=4, ncol=4)
LLH = matrix(nrow=4, ncol=4)

# We run different garch model fits and compare the test statistics to find the best model.
for (i in 1:4)
{
  for (j in 0:3)
  {
    sp_sum.dab <- garchFit(formula=substitute(~garch(alpha,beta), list(alpha = i, beta=j)), data=sp_sum.r, trace=F, cond.dist = c("norm"))
    AIC[i,j+1] <- sp_sum.dab@fit$ics[1]
    BIC[i,j+1] <- sp_sum.dab@fit$ics[2]
    LLH[i,j+1] <- sp_sum.dab@fit$params$llh
  }
}

AIC
min(AIC)
# - > GARCH(4,2)
BIC
min(BIC)
# - > GARCH(2,1)
LLH
max(LLH)
# - > GARCH(1,0)

sp_sum.d21 <- garchFit(formula=~garch(2,1), data=sp_sum.r, trace=F, cond.dist = c("norm"))

# Obtain results
summary(sp_sum.d21)
# -> alpha_2 is not significant different from zero, the AIC/BIC for a GARCH (2,1) and a GARCH(1,1) are almost identical, so we go back to a GARCH(1,1) model.

# We run an additional test, to see if we should add an ARMA part to our GARCH part.
for (i in 0:3)
{
  for (j in 0:3)
  {
    print(i)
    print(j)
    sp_sum.dab <- garchFit(formula=substitute(~arma(alpha,beta)+garch(1,1), list(alpha = i, beta=j)), data=sp_sum.r, trace=F, cond.dist = c("norm"))
    AIC[i+1,j+1] <- sp_sum.dab@fit$ics[1]
    BIC[i+1,j+1] <- sp_sum.dab@fit$ics[2]
    LLH[i+1,j+1] <- sp_sum.dab@fit$params$llh
    #summary(sp_sum.d21)
  }
}

AIC
min(AIC)
# - > ARMA(1,1) plus GARCH(1,1)
BIC
min(BIC)
# - > ARMA(0,1) plus GARCH(1,1)
LLH
max(LLH)
# - > ARMA(0,0) plus GARCH(1,1)

# -> We go for an ARMA(1,1)+GARCH(1,1) model.
sp_sum.d1111 <- garchFit(formula=~arma(1,1)+garch(1,1), data=sp_sum.r, trace=F, cond.dist = c("norm"))

# Obtain results
summary(sp_sum.d1111)
# -> alpha_2 is not significant different from zero, the AIC/BIC for a GARCH (2,1) and a GARCH(1,1) are almost identical, so we go back to a GARCH(1,1) model.

#residual diagnostic of ARMA(1,1)+GARCH(1,1)
par(mfrow=c(1,2))
acf(sp_sum.d1111@residuals/sp_sum.d1111@sigma.t,main="GARCH(1,1) residuals")
acf(sp_sum.d1111@residuals^2/sp_sum.d1111@sigma.t^2,main="GARCH(1,1) squared residuals")

#Jarque-Bera test for normality
jarque.bera.test((sp_sum.d1111@residuals/sp_sum.d1111@sigma.t))
# - > we reject normality assumption for the residuals

#Compare normal distribution function to SP residual distribution.
#par(mfrow=c(2,2))
a=density(sp_sum.d1111@residuals/sp_sum.d1111@sigma.t)
plot(a,main="ARMA(1,1)+GARCH(1,1) SP Residuals vs Normal")
lines(a$x,dnorm(a$x,mean=mean(sp_sum.d1111@residuals/sp_sum.d1111@sigma.t, na.rm=T), sd=sqrt(var(sp_sum.d11@residuals/sp_sum.d11@sigma.t, na.rm=T))),lty=2, col="red")
chart.QQPlot(sp_sum.d1111@residuals/sp_sum.d1111@sigma.t, main = "qq-Plot ARMA(1,1)+GARCH(1,1) SP Residual vs Normal", distribution = 'norm', col=c("black","red"))

#****************Asymmetries********************

# We investigate on asymmetries

a=c()
for (i in 1:length(sp_sum.r))
{
  a=c(a,max(sp_sum.r[i],0))
}
a

#Compute sample correlation (code from lecture)
for (h in 1:40)
{
  print(h)
  print(cor(a[(1+h):(length(a))],sp_sum.r[(1):(length(sp_sum.r)-h)])) #Leverage effect
}
# -> For the values lags 1,2,5,10,20,40 we find indication of leverage effects as all calculated sample correlations are negative

#Run a regression (code from lecture)
a=(sp_sum.r<0) #Sign Bias
a=c()
for (i in 1:length(sp_sum.r))
{
  #a=c(a,100*min(sp_sum.r[i],0)) #Negative Size Bias
  #a=c(a,100*max(sp_sum.r[i],0)) #Positive Size Bias
}

h=5
b=lm(100*sp_sum.r[(h+1):length(sp_sum.r)]^2~a[1:(length(sp_sum.r)-h)])
summary(b)
# -> We find for lag 1: sign bias:0.009021* , neg. sign bias:-0.0219***, pos. sign bias:0.006734***
# -> We find for lag 5: sign bias:0.007086. , neg. sign bias:-0.0263***, pos. sign bias:0.018003***
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# -> All significant at the 90% confidence level

# -> We need a model to take asymetries into account


#let's now extend our GARCH(1,1) to incorporate the asymmetries.

#We fit again GARCH(1,1), then EGARCH(1,1), GRJ-GARCH(1,1), T-GARCH(1,1) and P-GARCH(1,1) with p=2, all with an ARMA(1,1)
sp_sum.d1111=ugarchfit(
  data=sp_sum.r,
  spec=ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(1,1), arfima=FALSE), distribution.model="norm")
)

#EGARCH
sp_sum.d1111eg=ugarchfit(
  data=sp_sum.r,
  spec=ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(1,1), arfima=FALSE), distribution.model="norm")
)

#GJR-GARCH
sp_sum.d1111gg=ugarchfit(
  data=sp_sum.r,
  spec=ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(1,1), arfima=FALSE), distribution.model="norm")
)

#TGARCH (= PGARCH with d=1)
sp_sum.d1111tg=ugarchfit(
  data=sp_sum.r,
  spec=ugarchspec(variance.model=list(model="apARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(1,1), arfima=FALSE), fixed.pars=list(delta=1), distribution.model="norm")
)

#PGARCH (d=2)
sp_sum.d1111pg=ugarchfit(
  data=sp_sum.r,
  spec=ugarchspec(variance.model=list(model="apARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(1,1), arfima=FALSE), fixed.pars=list(delta=2), distribution.model="norm")
)

#We compute the AIC, BIC and LLH
AIC <- c(infocriteria(sp_sum.d1111)[1],
                     infocriteria(sp_sum.d1111eg)[1],
                     infocriteria(sp_sum.d1111gg)[1],
                     infocriteria(sp_sum.d1111tg)[1],
                     infocriteria(sp_sum.d1111pg)[1]
)
AIC
min(AIC)
# - > ARMA(1,1) plus PGARCH(1,1)

BIC <- c(infocriteria(sp_sum.d1111)[2],
         infocriteria(sp_sum.d1111eg)[2],
         infocriteria(sp_sum.d1111gg)[2],
         infocriteria(sp_sum.d1111tg)[2],
         infocriteria(sp_sum.d1111pg)[2]
)
BIC
min(BIC)
# - > ARMA(1,1) plus PGARCH(1,1)

LLH <- c(
         sp_sum.d1111@fit$LLH,
         sp_sum.d1111eg@fit$LLH,
         sp_sum.d1111gg@fit$LLH,
         sp_sum.d1111tg@fit$LLH,
         sp_sum.d1111pg@fit$LLH
)
LLH
max(LLH)
# - > ARMA(1,1) plus PGARCH(1,1)

# -> ARMA(1,1) plus PGARCH(1,1) with d=2 is the best model

#****************ARMA(1,1) plus PGARCH(1,1) with d=2 Analysis-closing prices********************

# Obtain results
infocriteria(sp_sum.d1111pg)
sp_sum.d1111pg@fit$matcoef

# -> mu is not significant different from zero, other coefficients are all different from 0 the 98% confidence level.

#residual diagnostic of ARMA(1,1) plus PGARCH(1,1) (d=2)
par(mfrow=c(1,2))
acf(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg),main="PGARCH(1,1) (d=2) residuals")
acf(sp_sum.d1111pg@fit$residuals^2/sigma(sp_sum.d1111pg)^2,main="PGARCH(1,1) (d=2) squared residuals")

# Now we have a look at the box test results
Box.test(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg),lag=5,type="Ljung")
Box.test(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg),lag=10,type="Ljung")
Box.test(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg),lag=20,type="Ljung")
Box.test(sp_sum.d1111pg@fit$residuals^2/sigma(sp_sum.d1111pg)^2,lag=5,type="Ljung")
Box.test(sp_sum.d1111pg@fit$residuals^2/sigma(sp_sum.d1111pg)^2,lag=10,type="Ljung")
Box.test(sp_sum.d1111pg@fit$residuals^2/sigma(sp_sum.d1111pg)^2,lag=20,type="Ljung")
# - > We cannot reject the null hypothesis of no independence

#Jarque-Bera test for normality
jarque.bera.test(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg))
# - > we reject normality assumption for the residuals

#Compare normal distribution function to SP residual distribution.
#par(mfrow=c(2,2))
a=density(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg))
plot(a,main="PARCH(1,1) (d=2) SP Residuals vs Normal")
lines(a$x,dnorm(a$x,mean=mean(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg), na.rm=T), sd=sqrt(var(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg), na.rm=T))),lty=2, col="red")
chart.QQPlot(sp_sum.d1111pg@fit$residuals/sigma(sp_sum.d1111pg), main = "qq-Plot PGARCH(1,1) (d=2) SP Residual vs Normal", distribution = 'norm', col=c("black","red"))

#News impact curve GARCH(1,1) vs the PGARCH(1,1)

ni=newsimpact(z = NULL, sp_sum.d1111)
ni2=newsimpact(z = NULL, sp_sum.d1111tg)
plot(ni$zx, ni$zy, ylab=ni$yexpr, xlab=ni$xexpr, type="l", main = "News Impact Curve")
lines(ni2$zx, ni2$zy, lty=2, col=2)
legend("topleft", leg=c("GARCH(1,1)", "PGARCH(1,1)"), lty=c(1,2), col=1:2, bg="white")

#****************Non-gaussian error distribution****************

#We consider a non-gaussian error distribution

#PGARCH (d=2)
sp_sum.d1111pg=ugarchfit(
  data=sp_sum.r,
  spec=ugarchspec(variance.model=list(model="apARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(1,1), arfima=FALSE), fixed.pars=list(delta=2), distribution.model="norm")
)

#PGARCH (d=2) with student-t distribution
sp_sum.d1111pg.t=ugarchfit(
  data=sp_sum.r,
  spec=ugarchspec(variance.model=list(model="apARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(1,1), arfima=FALSE), fixed.pars=list(delta=2), distribution.model="std")
)

#We compute the AIC, BIC and LLH
AIC <- c(infocriteria(sp_sum.d1111pg)[1],
         infocriteria(sp_sum.d1111pg.t)[1]
)
AIC
min(AIC)
# - > PGARCH(1,1) with t-student distribution

BIC <- c(infocriteria(sp_sum.d1111pg)[2],
         infocriteria(sp_sum.d1111pg.t)[2]
)
BIC
min(BIC)
# - > PGARCH(1,1) with t-student distribution

LLH <- c(
  sp_sum.d1111pg@fit$LLH,
  sp_sum.d1111pg.t@fit$LLH
)
LLH
max(LLH)
# - > PGARCH(1,1) with t-student distribution

# -> ARMA(1,1) plus PGARCH(1,1) with d=2 and student distribution is the best model

#****************ARMA(1,1) plus PGARCH(1,1) with d=2 and student distribution Analysis-closing prices********************

# Obtain results
infocriteria(sp_sum.d1111pg.t)
sp_sum.d1111pg.t@fit$matcoef

# -> mu is not significant different from zero, other coefficients are all different from 0 the 90% confidence level.

#residual diagnostic of ARMA(1,1) plus PGARCH(1,1) (d=2) with student distr.
par(mfrow=c(1,2))
acf(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t),main="PGARCH(1,1) (d=2) residuals")
acf(sp_sum.d1111pg.t@fit$residuals^2/sigma(sp_sum.d1111pg.t)^2,main="PGARCH(1,1) (d=2) squared residuals")

# Now we have a look at the box test results
Box.test(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t),lag=5,type="Ljung")
Box.test(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t),lag=10,type="Ljung")
Box.test(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t),lag=20,type="Ljung")
Box.test(sp_sum.d1111pg.t@fit$residuals^2/sigma(sp_sum.d1111pg.t)^2,lag=5,type="Ljung")
Box.test(sp_sum.d1111pg.t@fit$residuals^2/sigma(sp_sum.d1111pg.t)^2,lag=10,type="Ljung")
Box.test(sp_sum.d1111pg.t@fit$residuals^2/sigma(sp_sum.d1111pg.t)^2,lag=20,type="Ljung")
# - > We cannot reject the null hypothesis of no independence

#Jarque-Bera test for normality
jarque.bera.test(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t))
# - > we reject normality assumption for the residuals

#Compare normal distribution function to SP residual distribution.
#par(mfrow=c(2,2))
a=density(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t))
plot(a,main="PARCH(1,1) (d=2) SP Residuals vs Normal")
lines(a$x,dnorm(a$x,mean=mean(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t), na.rm=T), sd=sqrt(var(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t), na.rm=T))),lty=2, col="red")
chart.QQPlot(sp_sum.d1111pg.t@fit$residuals/sigma(sp_sum.d1111pg.t), main = "qq-Plot PGARCH(1,1) (d=2) SP Residual vs Normal", distribution = 'norm', col=c("black","red"))
