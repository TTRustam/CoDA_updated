

new_kt <- ts(fit_one$kt, start = 1965, end = 2018)


plot(diff(new_kt))
acf2(new_kt)
auto.arima(new_kt)
stepwise = FALSE
sarima(new_kt, p = 1, d = 1, q = 1)
sarima(new_kt, p = 0, d = 2, q = 1)


par(mfrow = c(2, 1))
sarima.for(window(new_kt, 1965, 2000), p = 1, d = 1, q = 1, 10)
lines(new_kt)
sarima.for(window(new_kt, 1965, 2000), p = 0, d = 2, q = 1, 10)
lines(new_kt)


# series removed from the excess
# add a honest description on why males do not work

checkresiduals()

accuracy()

autoplot() + autolayer(fitted())

ets()
# Simple eexponential smoothing ay(t) + a(1 - a)y(t-1) + a(1 - a) ^ 2 y(t-2), 0 =< a =< 1 
# a - how much weight is put on the recent observation and how quaicy weight decays away
# bigger a means more weight is put on a recent observation and weights decay very quickly
# small small weight on recent and decays slower
# Component from forecast equation : y(t+h|t) = l(t) // Smoothiing equation : l(t) = ay(t) + (1 - a)l(t-1)
# a and l(0) are two parametres estimated by regression minimizing the sum of squared errors ses function R 
# add trend with holt linear trend b(t) B* is how quincly the slope change, small means
# its hardly changes and almost linear trend. Thus mmodel has 4 parameters.
# phi for damping the larger the more damping it is if  == 1 equivalent for Holt method

par(mfrow = c(2, 1))
sarima.for(window(new_kt, 1980, 2014), p = 1, d = 1, q = 1, 10)
lines(new_kt)
sarima.for(window(new_kt, 1980, 2014), p = 0, d = 2, q = 1, 10)
lines(new_kt)


fcst_one <- forecast(Arima(window(new_kt, 1998, 2018),
                           order         = c(1, 1, 1) ,
                           include.drift = TRUE, 
                           method        = "ML"), h = ih)


fcst_two <- forecast(Arima(window(new_kt, 1998, 2018),
                           order         = c(0, 2, 1),
                           method        = "ML"), h = ih)



fcst_3 <- forecast(Arima(new_kt,
                           order         = c(1, 1, 1) ,
                           include.drift = TRUE, 
                           method        = "ML"), h = ih)


fcst_4 <- forecast(Arima(new_kt,
                           order         = c(0, 2, 1),
                           method        = "ML"), h = ih)

par(mfrow = c(4, 1))
plot(fcst_one)
plot(fcst_two)
plot(fcst_3)
plot(fcst_4)