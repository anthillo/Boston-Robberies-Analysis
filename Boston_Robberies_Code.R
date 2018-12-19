library(ggplot2)
library(forecast)

rob <- read.table("monthly.boston.robberies.txt", header = F)
rob.ts <- ts(rob, frequency = 12)
plot(rob.ts, xlab = "Time (Years) ", ylab = "Robberies in Boston", main = "Monthly Robberies in Boston 1966-1975")

acf(rob.ts, lag.max=50)
title("ACF of Original Data", line = 1)
pacf(rob.ts, lag.max=50)
title("PACF of Original Data", line = 1)

#decomposition plot of original data
decom.plot = decompose(rob.ts)
autoplot(decom.plot, main = "Decomposition Plot", xlab = "Time (years)")

#checking for seasonality component
seasonplot(rob.ts, 12, col=rainbow(12) ,year.labels=TRUE, main="Seasonal Plot")

#Stabilizing the increasing variance over time
# Box-Cox Transformation
require(MASS)
bcTransform <- boxcox(rob.ts~as.numeric(1:length(rob.ts)))

#square root transformation
plot(sqrt(rob.ts), xlab = "Time (Years) ", ylab = "sqrt(Robberies in Boston)", main = "Square-root Transformation of \n Monthly Robberies in Boston 1966-1975")

#cube root transformation
plot((rob.ts)^(1/3), xlab = "Time (Years) ", ylab = "cuberoot(Robberies in Boston)", main = "Cube-root Transformation of \n Monthly Robberies in Boston 1966-1975")

#Comparing the Variances of the initial time series with square-root and cube root transformations
var(rob.ts)
var(sqrt(rob.ts))
var((rob.ts)^(1/3))
# cube root transformation has the smallest variance

#new model with cube root transformation
cub.rob <- ts((rob.ts)^(1/3))

#De-trend the data 
robdiff1 <- diff(cub.rob, lag=1)

#Check the variance to see if it decreased
var(robdiff1)
#the variance did in fact decrease

#Difference the data once more to check the variance to see if it decrease
robdiff1diff1 <- diff(robdiff1, lag=1)
var(robdiff1diff1)
#the variance did not decrease so only revert to differencing data once

#Plot of Transformed and Differenced Data
plot(robdiff1, xlab = "Time(months)", ylab = "cuberoot(# Robberies)", main = expression(nabla[1]~"Transformed Data"))

#ACF of Transformed and Differenced Data
acf(robdiff1, lag.max=50)
title("ACF of Transformed and Differenced Data", line = 1)

#PACF of Transformed and Differenced Data
pacf(robdiff1, lag.max=50)
title("PACF of Transformed and Differenced Data", line = 1)

#AIC VALUES
aiccs <- matrix(NA, nr = 6, nc = 6)
dimnames(aiccs) = list(p=0:5, q=0:5)
for(p in 0:5)
{
  for(q in 0:5)
  {
    aiccs[p+1,q+1] = AIC(arima(robdiff1, order = c(p,1,q), method="ML"))
  }
}
aiccs

(aiccs==min(aiccs))

#BIC VALUES
biccs <- matrix(NA, nr = 6, nc = 6)
dimnames(aiccs) = list(p=0:5, q=0:5)
for(p in 0:5)
{
  for(q in 0:5)
  {
    biccs[p+1,q+1] = BIC(arima(robdiff1, order = c(p,1,q), method="ML"))
  }
}
biccs
(biccs==min(biccs))

#Model Estimation MLE method

#Model 1 - ARIMA(1,1,2)
fit1<- arima(cub.rob, order = c(1,1,2), method = "ML")
fit1

#Model 2 - ARIMA(2,1,2)
fit2<- arima(cub.rob, order = c(2,1,2), method = "ML")
fit2

source("plot.roots.txt")
par(mfrow = c(1,3))
#Roots of ARIMA(1,1,2)
fit1 <- arima(cub.rob, order = c(1,1,2), method="ML")
fit1
plot.roots(NULL, polyroot(c(1, -0.3718)), main = "Roots for AR Part - ARIMA(1,1,2)")
plot.roots(NULL, polyroot(c(1, 0.0913, -0.302)), main = "Roots for MA Part - ARIMA(1,1,2)")

#Roots for ARIMA(2,1,2)
fit2 <- arima(cub.rob, order = c(2, 1, 2), method = "ML")
fit2
plot.roots(NULL, polyroot(c(1, -0.2991, 0.1702)), main = "Roots for AR - ARIMA(2,1,2)")
plot.roots(NULL, polyroot(c(1, 0.0266, -0.4532)), main = "Roots for MA - ARIMA(2,1,2) ")

library(astsa)

model.1 <- arima(robdiff1, order = c(1,1,2))
resids1 = model.1$residuals

model.2 <- arima(robdiff1, order = c(2,1,2))
resids2 = model.2$residuals

## normality checks 

## histogram of residuals of ARMA(1,1,2)
hist(resids1, main="Histogram of Residuals of ARMA(1,1,2)", xlab = "Residuals")

## histogram of residuals of ARMA(2,1,2)
hist(resids2, main="Histogram of Residuals of ARMA(2,1,2)", xlab = "Residuals")

## qq plot of ARMA(1,1,2)
qqnorm(resids1, main = "Nornal Q-Q Plot - ARMA(1,1,2)")
qqline(resids1)

## qq plot of ARMA(2,1,2)
qqnorm(resids2, main = "Nornal Q-Q Plot - ARMA(2,1,2)")
qqline(resids2)

## Shapiro test for ARMA(1,1,2)
shap1 <- shapiro.test(resids1)
shap1$statistic # W Statistic
shap1$p.value # P-value

## Shapiro test for ARMA(2,1,2)
shap2 <- shapiro.test(resids2)
shap2$statistic # W Statistic
shap2$p.value # P-value

## Serial correlation check (Check for Independence)

## Box-Pierce ARIMA(1,1,2)
box_pierce <- Box.test(resids1, lag = 11, type = "Box-Pierce", fitdf = 2)$p.value
box_pierce
## Ljung-Box ARIMA(1,1,2)
ljung_box <- Box.test(resids1, lag = 11, type = "Ljung-Box", fitdf = 2)$p.value
ljung_box

## Box-Pierce ARIMA(2,1,2)
box_pierce <- Box.test(resids2, lag = 11, type = "Box-Pierce", fitdf = 2)$p.value
box_pierce
## Ljung-Box ARIMA(2,1,2)
ljung_box <- Box.test(resids2, lag = 11, type = "Ljung-Box", fitdf = 2)$p.value
ljung_box

## Constant variance Check

## ACF of residuals ARIMA(1,1,2)

acf(resids1, lag.max=500)
title("ACF Plot \n of Residuals - ARIMA(1,1,2)", line = -1)
## PACF of residuals ARIMA(1,1,2)
pacf(resids1, lag.max=500)
title("PACF Plot \n of Residuals - ARIMA(1,1,2)", line = -1)

## ACF of residuals ARIMA(2,1,2)
acf(resids2, lag.max=500)
title("ACF Plot \n of Residuals - ARIMA(2,1,2)", line = -1)
## PACF of residuals ARIMA(2,1,2)
pacf(resids2, lag.max=500)
title("PACF Plot \n of Residuals - ARIMA(2,1,2)", line = -1)

#forecasting 
fit = arima(cub.rob, order = c(1, 1, 2), method = "ML", xreg=1 : length(cub.rob)) 
pred.transf <- predict(fit, n.ahead = 12, newxreg=(length(cub.rob)+1) : length((cub.rob)+12))
# upper bound for the C.I. for transformed data
upper.transf= pred.transf$pred + 2*pred.transf$se 
# lower bound
lower.transf= pred.transf$pred - 2*pred.transf$se 
#plotting cub.rob and forecasting
ts.plot(cub.rob, xlim=c(1,length(cub.rob)+12), ylim = c(0,max(upper.transf)), xlab = "Time(months)", ylab = "cuberoot(Robberies in Boston)", main = "Forecasting Based on Transformed Data") 
lines(upper.transf, col="blue", lty="dashed")
lines(lower.transf, col="blue", lty="dashed")
points((length(cub.rob)+1):(length(cub.rob)+12), pred.transf$pred, col="red")

# RETURN TO ORIGINAL DATA
#get predictions and s.e's of transformed time series
rob.ts2 <- ts(rob)
# back-transform to get predictions of original time series
pred.orig <- pred.transf$pred^3
# bounds of the confidence intervals
upper = upper.transf^3
lower = lower.transf^3

# Plot forecasts using the original data
ts.plot(rob.ts2, xlim=c(1,length(rob.ts2)+12), ylim = c(0,max(upper)), main = "Forecasting Based on Original Data", xlab = "Time (months)", ylab = "Number of Robberies")
lines(upper, col="blue", lty="dashed")
lines(lower, col="blue", lty="dashed")
points((length(rob.ts2)+1):(length(rob.ts2)+12), pred.orig, col="red")

# Zoomed in  plot of  the last 12 values plus forecast:
ts.plot(rob.ts2, xlim=c(length(rob.ts2)-12,length(rob.ts2)+12), ylim = c(0,max(upper)), xlab = "Time (months)", ylab = "Number of Robberies", main = "Zoomed In Plot of Forecasting Original Data")
lines((length(rob.ts2)+1):(length(rob.ts2)+12),upper, lty=2, col="blue")
lines((length(rob.ts2)+1):(length(rob.ts2)+12),lower, lty=2, col="blue")
points((length(rob.ts2)+1):(length(rob.ts2)+12),pred.orig, col="red")
