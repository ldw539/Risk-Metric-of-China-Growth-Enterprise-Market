rm(list=ls())

library(TSA)

library(rugarch)

data1 <- read.csv("IDX 20100601 20100624.csv")

data2 <- read.csv("IDX 20100625 20150624.csv")

data3 <- read.csv("IDX 20150625 20200624.csv")

data <- rbind(data1,data2,data3)


inthesample <- data[1:2202,]
outthesample <- data[2203:2447, ]


drtn <- 0.01*inthesample$Idxtrd08


nlogr <- -log(drtn+1)

ts.plot(nlogr)

##---------riskmetic method---------## 
mspec = ugarchspec(mean.model = list(armaOrder=c(0,0),include.mean=
                                       FALSE), variance.model = list(model = "iGARCH",garchOrder=c(1,1)), 
                   fixed.pars = list(omega = 0))

fitmrm = ugarchfit(data = nlogr,spec = mspec)  ##fit igarch model

forecastrm <- ugarchforecast(fitmrm,n.ahead = 1)


sigma <- sigma(forecastrm)

VaRrm95 <- qnorm(0.95)*sigma

print("")
VaRrm95 <- as.numeric(VaRrm95)
VaRrm95

##------------econometric approach-----------------##
acf(nlogr)

Box.test(nlogr,lag=12)

t.test(nlogr)

library(fUnitRoots)
adfTest(nlogr,lags=10,type=c("c"))


eacf(nlogr)   ##eacf find ARMA order

econometricmodel <- arima(nlogr,order = c(1,0,1))
tsdiag(econometricmodel,gof.lag = 10)


residual <- econometricmodel$residuals

Box.test(residual^2,lag = 12,type = "Ljung")   #Ljung-Box test 

archtest <- function(mresiduals,m){
  y <- mresiduals^2
  T <- length(mresiduals)
  atsq = y[(m+1):T]
  x = matrix(0,(T-m),m)
  for (i in 1:m){ 
    x[,i]=y[(m+1-i):(T-i)] 
  }
  md = lm(atsq~x)
  summary(md)
}

archtest(residual,10)   #archtest 


mspec = ugarchspec(mean.model = list(armaOrder=c(1,1),include.mean=
                                       TRUE),variance.model= list(garchOrder=c(1,1)))

fitmea  = ugarchfit(data=nlogr,spec = mspec)

forecastea<- ugarchforecast(fitmea,n.ahead=1)


fitted(forecastea)   ## forecast conditional mean

sigma(forecastea)   ## forecast conditional standard deviation


forecastea

VaRea95  <- fitted(forecastea) + qnorm(0.95)*sigma(forecastea)

print("")

VaRea95 <- as.numeric(VaRea95)

##------------quantile estimation----------------##

#emperical quantile#

prob <- c(0.95)   #set the probability

VaRqe95 <- quantile(nlogr,prob)  # get the emperical quantile


print("")
VaRqe95 <- as.numeric(VaRqe95)


##-----------Extreme value approach--------------------##
library(evir)

n <- 28 #select the length of subperiod

m1 = gev(nlogr,block = n) #fit the  extreme value model

parameter <- m1$par.ests# Obtain the maximum likelihood estimates of shape
# parameter xi ,scale parameter simga ,location parameter mu 

xi <- parameter[1]

sigma <- parameter[2]

mu <- parameter[3]

if(parameter[1] == 0){
  VaRevt95 = mu - sigma*log(-n*log(0.95))
}else{
  VaRevt95 = mu - sigma*(1-(-n*log(0.95))^(-xi))/xi
}

print("")
VaRevt95 <- as.numeric(VaRevt95)
VaRevt95



##------------Peak over threshold-------------##
meplot(nlogr)   ##find the threshold

threshold <- 0.01

gpdmodel <- gpd(nlogr,threshold = 0.01)  #fit the generalized pareto distribution 

gpdxi <- gpdmodel$par.ests[1]

gpdbeta <- gpdmodel$par.ests[2]

VaRpot95 <- threshold -
  gpdbeta*(1-((length(nlogr)*0.05)/gpdmodel$n.exceed)^(-gpdxi))/gpdxi

print("")
VaRpot95 <- as.numeric(VaRpot95)
VaRpot95

##----------backtesting-----------------------##

HitSequence <- function(returns_X, VaR_X) {
  N = length(returns_X)
  Hit_X = numeric(N)
  Hit_X[which(returns_X > VaR_X)] = 1L
  return(Hit_X)
}


Christoffersen <- function(Hit, alpha) {
  N = length(Hit)
  
  n1 <- sum(Hit)
  
  n0 <- N-n1            # n0 and n1 are the number of 0s and 1s in the sample. 
  
  Lp <- (1-alpha)^(n0)*alpha^(n1)   #numberator 
  
  n00 = n01 = n10 = n11 = 0
  
  for (i in 2:N) {
    if (Hit[i] == 0L & Hit[i - 1L] == 0L)
      n00 = n00 + 1
    if (Hit[i] == 0L & Hit[i - 1L] == 1L)
      n01 = n01 + 1
    if (Hit[i] == 1L & Hit[i - 1L] == 0L)
      n10 = n10 + 1
    if (Hit[i] == 1L & Hit[i - 1L] == 1L)
      n11 = n11 + 1
  }
  pi01 = n01/(n00 + n01)
  pi11 = n11/(n10 + n11)
  pi = (n01 + n11)/(n00 + n01 + n10 + n11)
  
  
  
  LP1 = (1 - pi01)^n00 * pi01^n01 * (1 - pi11)^n10 *pi11^n11
  
  LRcc = -2*log(Lp/LP1)
  pvalue = 1 - pchisq(LRcc, df = 2L)
  LRcc = c(LRcc, pvalue)
  names(LRcc) = c("Test", "Pvalue")
  return(LRcc)
}


Backtest <- function(returns_X,VaR_X,alpha){
  Hit <- HitSequence(returns_X,VaR_X)
  LRcc <- Christoffersen(Hit,alpha)
  return(LRcc)
}


nlogrt <- -log(0.01*outthesample$Idxtrd08+1)


var1 <- rep(VaRrm95,times=length(outthesample[,1]))
Backtest(nlogrt,var1,0.05)



forecast <- ugarchforecast(fitmea,n.ahead = length(outthesample[,1]))
            
var2 <- fitted(forecast)+qnorm(0.95)*sigma(forecast)

var2 <- as.vector(var2)

Backtest(nlogrt,var2,0.05)




var3 <-rep(VaRqe95,times=length(outthesample[,1]))

Backtest(nlogrt,var3,0.05)


var4 <- rep(VaRevt95,times=length(outthesample[,1]))
Backtest(nlogrt,var4,0.05)



var5 <- rep(VaRpot95,times=length(outthesample[,1]))

Backtest(nlogrt,var5,0.05)