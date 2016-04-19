# weibullpp
An R package to mirror functionality in Weibul++

* right censoring only
* weibull and exponential distribution only
* fitting by MLE
* prob. of failure and reliability vs time plots

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
plot data, type can be failure or reliability
```R
plot(y, type = 'failure')
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

