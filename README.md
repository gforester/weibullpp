# weibullpp 
weibullpp is an R package designed for reliability engineering.
It is meant to add to community of packages in R for survival analysis.
[CRAN Task View: Survival Analysis](https://cran.r-project.org/web/views/Survival.html)
Another goal is to offer an open-source alternative to commercial software such as Reliasoft's Weibull++.
To this end, many plots and calculations available in Weibull++ will be made available in the package.

As of version 0.2, the package is restricted to the following
* life data analysis 
* right censoring
* weibull and exponential distribution
* parameteric fitting by MLE

## Life Data Analysis
The `lifedata()` function will take the input times and censoring status to create a `lifedata` object.
There are two required inputs: time-to-event and censoring status, respectively.
Status should be a 0 for right censored, otherwise 1.
An optional 3rd parameter for units of time values is available not used in any way.
Example:
```R
x = lifedata(c(90.7, 114.8, 12.0, 144.35, 199.8), c(0,1,1,1,0), 'hours')
```
To fit data to a distribution, use `lifedata()`. 
The result is a `fitted_life_data` object.
This object constains the original data, the distribution fitted, the fitted values, and goodness-of-fit values.
Note: you may see warnings outputted by the optimization routine when fitting Weibulls.
Example:
```R
y =  fit_data(x) #fit weibull distribution by default
y_exp=fit_data(x, dist = 'exponential') #fit exponential distribution
```
The package provides several plot types of the fitted distribution.
These include: probability on linearized scales, CDF (unreliability), reliability, PDF, and failure rate (hazard function).
The first input to the `plot()` function should be a `fitted_life_data` object.
The second input is a string to denote the plot type.
Linearized scale plots are specific to the fitted distribution and attempt to follow reliability conventions.

The first 3 plots include points representing the calculated median ranks of failure data.
It is referred to as the "standard ranking method" in Weibull++.

Examples:
```R
plot(y) #default to linearized scale plot
plot(y, type = 'failure') #probability of failure over time
plot(y, type = 'reliability') #reliability over time
plot(y, type = 'pdf') #pdf of fitted distribution
```
For styling purposes, an optional `theme` input can accept values of either `'base_r'` or `'weibull++`.
The latter produces a plot that looks similar to those produced by Weibull++.
As of v0.2, graphical parameters can NOT be passed; we will strive to have this functionality in the next release.
```R
plot(y, theme = 'weibull++') #linearized scale plot with Weibull++ styling
```

Histograms can be plotted of the failures or suspension points.
The function can take either `lifedata` objects or `fitted_life_data` objects.
At this time, they do NOT accept `weibull++` theme inputs.
NOTE: this functionality is in its early infancy.
Examples:
```R
hist(x) #input from lifedata object
hist(y, type = 'failure') #histogram of failure points, input is fitted_life_data_object
hist(y, type = 'suspension') #histogram of suspension points
```

Calculation of metrics on fitted data is done through the `calculate(...)` function.
See example for available options.
```R
fit = fit_data(x, dist = 'weibull') #input to calculate is the result of fit_data
calculate(fit, 'reliability', 100) #probability of surviving to t=100
calculate(fit, 'failure', 100) #probability of failure before t=100
calculate(fit, 'mean life', NA) #mean time to failure
calculate(fit, 'failure rate',100) #aka hazard function; instantaneous number of failures per unit time at t = 100
```

## Recurrent Event Analysis
This is a skeleton in what we hope the actual functionality will be.

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

