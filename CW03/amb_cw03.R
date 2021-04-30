## Question 2 ----
# load libraries and data ----
library(forecast)
library(ggplot2)

pine <- scan('pine.dat', skip = 1) 
pine.ts <- ts(pine, start = 1107, frequency = 1)
# subset from 1201 to 1500
pine.ts <- window(pine.ts, start = 1201, end = 1500)
length(pine.ts)

# visualize time plots ----
cbind("Year Rings of Douglas Fir" = pine.ts,
      "Year Rings of Douglas Fir \n (logged)" = log(pine.ts)) %>%  
  autoplot(facets=TRUE) +
  xlab("Year") + ylab("")

# ACF/PACF ----
# lag=1
plot.acf <- ggAcf(diff(log(pine.ts),1)) + ggtitle('differenced')    # q=22
plot.pacf <- ggPacf(diff(log(pine.ts),1)) + ggtitle('differenced')   # p=18
plot.acf   
plot.pacf  
# lag=10
plot.acf2 <- ggAcf(diff(log(pine.ts),10)) + ggtitle('differenced')
plot.pacf2 <- ggPacf(diff(log(pine.ts),10)) + ggtitle('differenced')
plot.acf2 
plot.pacf2  

pine.ts %>% log() %>% diff(10) %>% ggtsdisplay(main="")

# transformation ----
pine.log <- log(pine.ts)
frequency(pine.log)   # =1 (year)

# fit SARIMA on log + 1st diff(10-year) ----
# auto
auto.arima(pine.log)

# SARIMA(2,0,3)(2,1,1)^10
fit <- pine.log %>%
        Arima(order=c(2,0,3),
              seasonal = list(order = c(2,1,1), period = 10)) 
fit   # 336.44
checkresiduals(fit)

par(mfrow = c(2,2))
plot(resid(fit))
qqnorm(resid(fit), main = ''); qqline(resid(fit))
acf(resid(fit))
pacf(resid(fit))

# SARIMA(2,0,2)(2,1,1)^10
fit2 <- pine.log %>%
        Arima(order=c(2,0,2),
              seasonal = list(order = c(2,1,1), period = 10)) 
fit2   # 335.37
checkresiduals(fit2)

# SARIMA(2,0,8)(2,1,1)^10
fit3 <- pine.log %>%
        Arima(order=c(2,0,8),
              seasonal = list(order = c(2,1,1), period = 10)) 
fit3   # 342.68
checkresiduals(fit3)

## change AR's non-seasonal component
# SARIMA(1,0,3)(2,1,1)^10
fit4 <- pine.log %>%
        Arima(order=c(1,0,3),
              seasonal = list(order = c(2,1,1), period = 10)) 
fit4   # 334.08
checkresiduals(fit4)

# SARIMA(3,0,3)(2,1,1)^10
fit5 <- pine.log %>%
        Arima(order=c(3,0,3),
              seasonal = list(order = c(2,1,1), period = 10)) 
fit5   # 328.94
checkresiduals(fit5)

# SARIMA(8,0,3)(2,1,1)^10
fit6 <- pine.log %>%
        Arima(order=c(8,0,3),
              seasonal = list(order = c(2,1,1), period = 10)) 
fit6   # 334.52
checkresiduals(fit6)

# change SEASONAL component
# SARIMA(3,0,3)(2,1,0)^10
fit7 <- pine.log %>%
        Arima(order=c(3,0,3),
              seasonal = list(order = c(2,1,0), period = 10)) 
fit7   # 386.03
checkresiduals(fit7)

# SARIMA(3,0,3)(2,1,2)^10
fit8 <- pine.log %>%
        Arima(order=c(3,0,3),
              seasonal = list(order = c(2,1,2), period = 10)) 
fit8   # 386.03
checkresiduals(fit8)

# SARIMA(3,0,3)(1,1,1)^10
fit9 <- pine.log %>%
        Arima(order=c(3,0,3),
              seasonal = list(order = c(1,1,1), period = 10)) 
fit9   # 332.93
checkresiduals(fit9)

# compare AIC ----
aic <- AIC(fit, fit2, fit3, 
    fit4, fit5, fit6,
    fit7, fit8, fit9)
round(aic[,2],3)

# check residuals of the chosen model ----
fit5   # 328.94
checkresiduals(fit5)

# forecasting ----
frequency(fit5)
fit5 %>% 
        forecast(h = 50, level = c(80,95)) %>%
        autoplot() +
        ylab("Year Rings of Douglas Fir (log)") + xlab("Year")

# get forecast values
forecast <- fit5 %>% forecast(h = 50, level = c(80,95))
summary(forecast)
forecast.log <- forecast$mean
forecast.log
forecast.ts <- exp(forecast.log)
forecast.ts
# get test set
pine.ts.original <- ts(pine, start = 1107, frequency = 1)
test <- window(pine.ts.original, start = 1501, end = 1550)
test
# calculate RMSE
library(ModelMetrics)
rmse(test, forecast.ts)   # = 37.50559

# calculate all the RMSEs ----
func.rmse <- function(fit, test) {
        forecast <- fit %>% forecast(h = 50)
        forecast.log <- forecast$mean
        forecast.ts <- exp(forecast.log)
        rmse <- rmse(test, forecast.ts)
        return(rmse)
}

rmse <- rep(0, 9)
rmse[1] <- func.rmse(fit, test)
rmse[2] <- func.rmse(fit2, test)
rmse[3] <- func.rmse(fit3, test)
rmse[4] <- func.rmse(fit4, test)
rmse[5] <- func.rmse(fit5, test)
rmse[6] <- func.rmse(fit6, test)
rmse[7] <- func.rmse(fit7, test)
rmse[8] <- func.rmse(fit8, test)
rmse[9] <- func.rmse(fit9, test)
rmse
min(rmse)
