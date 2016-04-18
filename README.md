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


