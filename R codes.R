rm(list=ls())
library(TSA)
library(fUnitRoots)
library(forecast)
library(tseries)
library(lmtest)
library(magrittr)
library(readr)
library(car)

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

# Import the data set
bc <- read_csv("dataset.csv",col_names=TRUE)
bc = bc$Bitcoin
bc %>% str() # Total observation 162

# Convert into a time series object 
bc.ts = ts(bc,start=c(2011,08), frequency=12)
bc.ts %>% class()
bc.ts %>% head()

bc.ts %>% summary()
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#     3988    470981   7930612  19781743  31599541 128015000  

par(mar = c(4, 4, 2.5, 1),  # Adjust margin of plots
    mgp = c(2, 0.5, 0))  # Adjust xlab/ylab position

# Time series plot of raw series
plot(bc.ts,type='o',ylab='Bitcoin Price', xlab='Year', main='Figure 1: Time series plot of Monthly Bitcoin Price Series')
# Upwards quadratic trend exists,  
# No seasonality, change point
# slightly changing variance, but mostly autoregressive behaivor


par(mfrow=c(1,2),
    mar = c(3, 3, 4, 1),  # Adjust margin of plots
    mgp = c(1.5, 0.5, 0)) # Adjust xlab/ylab position
acf(bc.ts, main = "Figure 2: ACF of Bitcoin Price Change") # slowly decaying pattern ->non stationary
pacf(bc.ts, main = "Figure 3: PACF of Bitcoin Price Change") # high first lag ->non stationary
par(mfrow=c(1,1))

# To confirm the stationarity
adf.test(bc.ts) # p-value > 0.05 ->non stationary
pp.test(bc.ts) # p-value > 0.05 ->non stationary
kpss.test(bc.ts) # p-value < 0.05 ->non stationary

# Check normality
qqnorm(bc.ts, ylab="Monthly Bitcoin Price", xlab="Normal Scores", main = "Figure 4: QQ plot of Monthly Bitcoin Price") #Check normality
qqline(bc.ts, col = 2, lwd = 1, lty = 2) # data is not normal distribution
shapiro.test(bc.ts)  # data is not normal distribution

# Box-Cox transformation
BC <- BoxCox.ar(bc.ts, lambda = seq(-1, 1.125, 0.05))
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda 
newBC.ts <- ((bc.ts^lambda)-1)/lambda

# optimal lambda is 0.05
# The optimal lambda value is 0.05 which is close to log transformation

# Time series plot of BoxCox transformed Bitcoin price
plot(newBC.ts,type='o', ylab ="Bitcoin Price", main="Figure 6: Time Series Plot of Box-Cox Transformed Bitcoin Price")

par(mfrow=c(1,2),
    mar = c(3, 3, 4, 1),  # Adjust margin of plots
    mgp = c(1.5, 0.5, 0))
acf(newBC.ts, main = "Figure 7: ACF of Box-Cox Transformed Series") 
pacf(newBC.ts, main = "Figure 8: PACF of Box-Cox Transformed Series")
par(mfrow=c(1,1))
# Upward quadratic trend exists 
adf.test(newBC.ts) # p-value > 0.05 ->non stationary
pp.test(newBC.ts) # p-value > 0.05 ->non stationary
kpss.test(newBC.ts) # p-value < 0.05 ->non stationary

# Check normality again
qqnorm(newBC.ts, main = "Figure 9: QQ plot of Box-Cox Transformed Bitcoin Price", ylab="BC transformed Monthly Bitcoin Price", xlab="Normal Scores")
qqline(newBC.ts,col = 2, lwd = 1, lty = 2)  # slightly improvement of the normality, but both tails are still off
shapiro.test(newBC.ts) # p<0.05 -> data is not normal distribution

# 1st differencing to BoxCox transformed series, because of the improvement of normality
diff.bc.ts = diff(newBC.ts, differences = 1)
plot(diff.bc.ts, type='o',ylab='First Difference of monthly Bitcoin price', 
     main ="Figure 10: First difference of Box-Cox Transformed Bitcoin Price Series") # flat mean level -> stop differencing

par(mfrow=c(1,2),
    mar = c(3, 3, 4, 1),  # Adjust margin of plots
    mgp = c(1.5, 0.5, 0))
acf(diff.bc.ts, main = "Figure 11: ACF of First Differencing") 
pacf(diff.bc.ts, main = "Figure 12: PACF of First Differencing")
par(mfrow=c(1,1))

# To confirm the stationarity of the first differenced series
adf.test(diff.bc.ts) # p < 0.05 -> stationary
pp.test(diff.bc.ts) # p < 0.05 -> stationary
kpss.test(diff.bc.ts) # p > 0.05 -> stationary

# Model specification with the first differenced series (d = 1)

# ACF & PACF
par(mfrow=c(1,2))
acf(diff.bc.ts, main = "Figure 9: ACF of First Differencing") 
pacf(diff.bc.ts, main = "Figure 10: PACF of First Differencing")
par(mfrow=c(1,1))

# From ACF & PACF, p = 2 & q = 2
# possible models: {ARIMA(2,1,2)}

# EACF
eacf(diff.bc.ts)
# From EACF, (0,1), (1,0), (1,1)
# possible models: {ARIMA(0,1,1), ARIMA(1,1,0), ARIMA(1,1,1)}

#Fig 13: BIC table - from ACF, PACF & EACF, max p & q = 2, so set the limit to 5
res = armasubsets(y=diff.bc.ts, nar=5, nma=5,
                  y.name='p',ar.method='ols')
plot(res)
# From BIC table, (1,0) from first row, (1,4) from second row, (5,4) from third row
# possible models: {ARIMA(1,1,0), ARIMA(1,1,4), ARIMA(5,1,4)}

# Hence, Final set of possible model are: 
#{ARIMA(2,1,2), ARIMA(0,1,1), ARIMA(1,1,0), ARIMA(1,1,1), ARIMA(1,1,4), ARIMA(5,1,4)}

# ARIMA(2,1,2) ***
model.212 = Arima(bc.ts, order=c(2,1,2), method='ML')
coeftest(model.212)

model.212CSS = Arima(bc.ts, order=c(2,1,2), method='CSS')
coeftest(model.212CSS)

model.212CSSML = Arima(bc.ts, order=c(2,1,2), method='CSS-ML')
coeftest(model.212CSSML)

# ARIMA(0,1,1) 
model.011 = Arima(bc.ts, order=c(0,1,1), method='ML')
coeftest(model.011)

model.011CSS = Arima(bc.ts, order=c(0,1,1), method='CSS')
coeftest(model.011CSS)

# ARIMA(1,1,0) 
model.110 = Arima(bc.ts,order=c(1,1,0), method='ML')
coeftest(model.110)

model.110CSS = Arima(bc.ts,order=c(1,1,0), method='CSS')
coeftest(model.110CSS)

# ARIMA(1,1,1)
model.111 = Arima(bc.ts,order=c(1,1,1), method='ML')
coeftest(model.111)

model.111CSS = Arima(bc.ts,order=c(1,1,1), method='CSS')
coeftest(model.111CSS)

# ARIMA(1,1,4) 
model.114 = Arima(bc.ts,order=c(1,1,4), method='ML')
coeftest(model.114)

model.114CSS = Arima(bc.ts,order=c(1,1,4), method='CSS')
coeftest(model.114CSS)

# ARIMA(5,1,4) 
model.514 = Arima(bc.ts,order=c(5,1,4), method='ML')
coeftest(model.514)

model.514CSS = Arima(bc.ts,order=c(5,1,4), method='CSS')
coeftest(model.514CSS)

# Goodness-of-fit metrics 
# Check AIC and BIC values of the models to decide the best one within the subset of possible models.
sort.score(AIC(model.212, model.011, model.110, model.111, model.114, model.514), score = "aic")


sort.score(BIC(model.212, model.011, model.110, model.111, model.114, model.514), score = "bic")

# Hence, Final set of possible model are: 
#{ARIMA(0,1,1), ARIMA(1,1,0), ARIMA(1,1,1), ARIMA(1,1,4), ARIMA(5,1,4}

Smodel_212_css <- accuracy(model.212)[1:7]
Smodel_011_css <- accuracy(model.011)[1:7]
Smodel_110_css <- accuracy(model.110)[1:7]
Smodel_111_css <- accuracy(model.111)[1:7]
Smodel_114_css <- accuracy(model.114)[1:7]
Smodel_514_css <- accuracy(model.514)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_212_css, Smodel_011_css, Smodel_110_css, Smodel_111_css, Smodel_114_css, Smodel_514_css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(2,1,2)", "ARIMA(0,1,1)", "ARIMA(1,1,0)", "ARIMA(1,1,1)", 
                          "ARIMA(1,1,4)", "ARIMA(5,1,4)")
round(df.Smodels,  digits = 3)

#                   ME    RMSE     MAE   MPE   MAPE  MASE   ACF1
#ARIMA(2,1,2) 718956.8 5691449 2822234< 2.497 18.725 0.217<  0.003
#ARIMA(0,1,1) 710175.5 5855822 2945637 2.550 17.415< 0.227 -0.008
#ARIMA(1,1,0) 695133.5 5853895 2950675 2.522 17.443 0.227 -0.012
#ARIMA(1,1,1) 693599.2 5853884 2951392 2.518 17.446 0.227 -0.014
#ARIMA(1,1,4) 654373.5 5684633< 2987766 2.448 19.358 0.230 -0.015
#ARIMA(5,1,4) 660717.9 5392454 2829687 2.578 19.371 0.218 -0.014

# The best model is ARIMA(2,1,2). So, the overparametrised models 
# are ARIMA(3,1,2) and ARIMA(2,1,3)

# ARIMA(3,1,2)
model.312 = Arima(bc.ts, order=c(3,1,2), method='ML')
coeftest(model.312)

model.312CSS = Arima(bc.ts, order=c(3,1,2), method='CSS')
coeftest(model.312CSS)

model.312CSSML = Arima(bc.ts, order=c(3,1,2), method='CSS-ML')
coeftest(model.312CSSML)

# ARIMA(2,1,3)
model.213 = Arima(bc.ts, order=c(2,1,3), method='ML')
coeftest(model.213)

model.213CSS = Arima(bc.ts, order=c(2,1,3), method='CSS')
coeftest(model.213CSS)

model.213CSSML = Arima(bc.ts, order=c(2,1,3), method='CSS-ML')
coeftest(model.213CSSML)

# Goodness-of-fit metrics 

sort.score(AIC(model.212, model.312, model.213), score = "aic")
sort.score(BIC(model.212, model.312, model.213), score = "bic")

Smodel_212_css <- accuracy(model.212)[1:7]
Smodel_312_css <- accuracy(model.312)[1:7]
Smodel_213_css <- accuracy(model.213)[1:7]
df.Smodels <- data.frame(
  rbind(Smodel_212_css, Smodel_312_css, Smodel_213_css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(2,1,2)", "ARIMA(3,1,2)", "ARIMA(2,1,3)")
round(df.Smodels,  digits = 3)

# By checking the overparametrised models, we could conclude that the best model is ARIMA(2,1,2).