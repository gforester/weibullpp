# weibullpp
An R package to mirror functionality in Weibul++

* right censoring only
* weibull and exponential distribution only
* fitting by MLE

## Life Data Analysis
Input times and censoring status. 
Status should be a 0 for right censored, otherwise 1.
Optional 3rd parameter for units of time values.
```R
x = lifedata(time, status, units)
```
Fit data to distribution. First input should be result of `lifedata()`.
```R
y = fit_data(x, dist = 'weibull')
```
plot data examples
```R
plot(y, type = 'failure') #probability of failure over time
plot(y, type = 'reliability') #reliability over time
plot(y, type = 'pdf') #pdf of fitted distribution
hist(y, type = 'failure') #histogram of failure points
hist(y, type = 'suspension') #histogram of suspension points
```

Calculation of metrics on fitted data is done through the `calculate(...)` function.
```R
fit = fit_data(x, dist = 'weibull') #input to calculate is the result of fit_data
calculate(fit, 'reliability', 100) #probability of surviving to t=100
calculate(fit, 'failure', 100) #probability of failure before t=100
calculate(fit, 'mean life', NA) #mean time to failure
calculate(fit, 'failure rate',100) #aka hazard function; instantaneous number of failures per unit time at t = 100
```

## Recurrent Event Analysis
At this time, only really set up for 1 unique ID.
General Renewal Process (GRP) is the assumed parameteric model with power law model.

Input: times, ID, and status.
Optional Input: time units.
Status and time units input should be similar to Life Data Analysis.
```R
x = recurrentdata(time, id, status, units)
plot(x)
y = mle_recurrent(x, type = 1)
```
Plotting data will result in a step plot of Number of occurences vs time.
`mle_recurrent` will default to type I but also accepts type II model.

