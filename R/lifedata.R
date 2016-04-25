
#' Create lifedata Object.
#' 
#' @param time numerical vector indicating time-to-event
#' @param status numerical vector indicating status. 0 for right censored. 1 for failure event.
#' @param units optional character to denote units of time
#' @return a \code{lifedata} object
#' @examples
#' lifedata(c(90.7, 114.8, 12.0, 144.35, 199.8), c(0,1,1,1,0), 'hours')
lifedata <- function(time, status, units = NULL){
  #check that time is positive
  if(any(time < 0)) stop('time values must be non-negative')
  #check that status is 1 or 0
  if(any(!((status == 0) | (status == 1)))) stop('status values must be either 0 or 1')
  
  new_lifedata <- list(time = time, status = status)
  class(new_lifedata) <- 'lifedata'
  attr(new_lifedata,'units') <- units
  
  new_lifedata
}
print.lifedata <- function(x){
  pos_or_not <- rep(' ',length(x[[1]]))
  pos_or_not[x[[2]] == 0] <- '+'
  print(paste0(x[[1]],pos_or_not), quote = F)
}

#' @title Calculate Median Ranks
#' 
#' @description 
#' s3 generic
#' 
#' @details
#' refer to median_ranks.lifedata
#' 
#' @param x only method defined for lifedata objects
median_ranks <- function(x, ...) UseMethod('median_ranks')
#' @title Calculate Median Ranks
#' 
#' @description 
#' calculates median ranks of lifedata object.
#' 
#' @details
#' Median ranks are based on. 
#' 
#' By default, ranks are calculated based on solving 
#' 
#' For censored data, the mean order number of the failure points are determined first.
#' Median ranks are only calculated for failure points.
#' 
#' @param x a lifedata object
#' @param method indicate calculation method. see details
#' @return numeric vector corresponding to median ranks for failure data
median_ranks.lifedata <- function(x, method = 'inverse_F'){
  N_total <- length(x$time)         #total number of times
  #ties in rank use default: 'average'
  if(any(x$status == 0)){
    raw_rank <- rank(x$time)          #rank all times, regardless of status
    idx <- x$status == 1              #T/F, T for failure
    temp <- rep(NA, length(raw_rank)) #resulting rank will be NA for censored data
    
    order_failure <- sort(raw_rank[idx])  #rank of just failures
    for(i in seq_along(order_failure)){
      if(i != 1){
        prev <- temp[order_failure[i-1]]  
      }else{
        prev <- 0
      }
      N_beyond <- sum(x$time >= x$time[order_failure[i]])
      temp[order_failure[i]] <- prev + ( (N_total+1-prev) / (1+N_beyond) )
    }
    #for compatibility, pass only non-NA items. in future, redo
    temp <- temp[!is.na(temp)]
  }else{
    temp <- rank(x$time)  
  }

  # METHODS to find median rands
  if(method == 'inverse_F'){
    ms <- 2*(N_total-temp+1)
    ns <- 2*temp
    Fstars <- qf(0.5, ms, ns, lower.tail = F) 
    mrs <- 1/(1+(Fstars*(N_total-temp+1)/temp))
  }
  if(method == 'inverse_binom'){ #should get equivalent to inverse_F w.o. optimization needed.
    mrs <- sapply(X = temp, FUN = function(y) solve_for_p(N_total, y))  
  }
  if(method == 'benard'){
    mrs <- (temp-0.3)/(N_total+0.4)
  }
  return(mrs)
}
solve_for_p <- function(N, i, cdf = 0.5){
  func_to_optim <- function(prob){
    (sum((choose(N,i:N)*prob^(i:N)*(1-prob)^(N-(i:N))))-cdf)^2
  }
  optimize(f = func_to_optim, interval = c(0,1))$minimum
}

#' @title Fit Distribution to lifedata Object
#' 
#' @description 
#' Fits a distribution to lifedata object.
#' Only MLE are done and only weibull and exponential distributions are available
#' 
#' @details
#' For exponential distributions, uses analytic solution.
#' For Weibull distributions, use optim function to find MLE.
#' 
#' @param x a lifedata object
#' @param dist character indicating distribution. only exponential and weibull are available
#' @return 
#' A fitted_life_data object.
#' Object includes original input lifedata object, the distribution fitted, log likelihood, parameter estimates.
#' If all data points are complete, a chi-square goodness-of-fit output is included
#' 
fit_data.lifedata <- function(x, dist = 'weibull'){
  if(!((dist == 'weibull') | (dist == 'exponential'))) stop('dist must be "weibull" or "exponential"')
  
  if(dist == 'exponential'){
    res <- fit_exp(x)
    mle_cdf <- function(z){
      pexp(z, 1/res$scale)
    }
    n_params <- 1
  }else{
    res <- weibull_mle(x)
    mle_cdf <- function(z){
      pweibull(z, res$shape, res$scale)
    }
    n_params <- 2
  }
  
  if(all(x$status==1)){
    chisq_test_res <- chi_square_test(x$time, mle_cdf, n_params)
    gof = list(chisq = chisq_test_res)
  }else{
    gof = "only available when all data complete"
  }
    
  to_return <- list(data = x, dist = dist, fit = res, gof = gof)
  class(to_return) <- 'fitted_life_data'
  
  return(to_return)
}

#' @title Plotting fitted_life_data Objects
#' 
#' @description 
#' Provides various plots of fitted distribution
#' 
#' @details
#' Values for type input
#' \tabular{ll}{
#' probability   \tab linearized scale plot. distribution specific.\cr
#' failure \tab probability of failure over time i.e. CDF\cr
#' reliability \tab reliability over time i.e. 1-CDF \cr
#' pdf \tab fitted probability density function \cr
#' failure rate \tab fitted failure rate i.e. hazard function \cr
#' }
#' First 3 plots will included points corresponding to the 
#' median ranks of the failure data.
#' 
#' By convention, in linearized scales, the exponential distribution plots the reliability while 
#' Weibull plots the unreliability.
#' 
#' When theme = 'weibull++', the plots styling is automatically set to mirror Weibull++.
#' By default, theme = 'base_r' and plots a very simple plot.
#' At this time, no graphical parameters can be passed.
#' 
#' @param x a fitted_life_data object
#' @param type type of plot desired. default to 'probability'. see details
#' @param theme character indicating style. only 'base_r' and 'weibull++' are available. see details.
#' @return 
#' produces desired plot
plot.fitted_life_data <- function(x, type = 'probability', theme = 'base_r', 
                                  line_par = NULL, point_par = NULL, ...){
  if(type == 'probability'){
    if(x$dist == 'weibull') {plot_weibull(x, theme, line_par, point_par, ...)}
    if(x$dist == 'exponential'){plot_exponential(x,theme, line_par, point_par, ...)}
  }
  if(type == 'failure'){
    plot_cdfs(x, lower.tail = T, theme, line_par, point_par, ...)
  }
  if(type == 'reliability'){
    plot_cdfs(x, lower.tail = F, theme, line_par, point_par, ...)
  }
  if(type == 'pdf'){
    plot_pdf(x, theme, ...)
  }
  if(type == 'failure rate'){
    plot_hazard(x, theme, ...)
  }
}
plot_cdfs <- function(x, lower.tail, theme, line_par, point_par, ...){
  xmax <- 1.02 * max(x$data[[1]])
  xs <- seq(from=0, to=xmax, length.out = 5E2)
  if(x$dist == 'exponential'){
    if(lower.tail){
      ys <- pexp(xs, rate = 1/x$fit$scale)  
    }else{
      ys <- pexp(xs, rate = 1/x$fit$scale, lower.tail = F)
    }
  }
  if(x$dist == 'weibull'){
    if(lower.tail){
      ys <- pweibull(xs, shape = x$fit$shape, scale = x$fit$scale)  
    }else{
      ys <- pweibull(xs, shape = x$fit$shape, scale = x$fit$scale, lower.tail = F) 
    }
  }
  
  #get points to plot
  pts_ys <- median_ranks(x$data)
  if(!lower.tail) pts_ys <- 1 - pts_ys
  pts_xs <- x$data$time[x$data$status == 1]
    
  #plot based on theme
  if(theme == 'base_r'){
    do.call(plot,c(list(x = xs, y = ys, type = 'l'), line_par, ...))
    do.call(points, c(list(x = pts_xs, y = pts_ys), point_par, ...))
    # plot(xs, ys, type = 'l', ...)
    # points(pts_xs, pts_ys, ...)
  }
  if(theme == 'weibull++'){
    if(lower.tail){
      cdf_title = 'Unreliability vs. Time'
      cdf_y_lab = 'Unreliability'
    }else{
      cdf_title = 'Reliability vs. Time'
      cdf_y_lab = 'Reliability'
    }
    weibull_theme_plot(xs, ys, cdf_title, 'Time', cdf_y_lab)
    points(pts_xs, pts_ys, pch = 16, col='blue')
  }
}
plot_pdf <- function(x, theme, ...){
  xmax <- 1.02 * max(x$data[[1]])
  xs <- seq(from=0, to=xmax, length.out = 5E2)
  if(x$dist == 'exponential'){
    ys <- dexp(xs, rate = 1/x$fit$scale)
  }
  if(x$dist == 'weibull'){
    ys <- dweibull(xs, shape = x$fit$shape, scale = x$fit$scale)
  }
  
  if(theme == 'base_r'){
    plot(xs, ys, type = 'l', ...)  
  }
  if(theme == 'ggplot'){
    plot(xs, ys, type = 'l') 
  }
  if(theme == 'weibull++'){
    weibull_theme_plot(xs, ys, 'Probability Density Function', 'Time', 'f(t)')
  }
  
}
plot_hazard <- function(x, theme, ...){
  xmax <- 1.02 * max(x$data[[1]])
  xs <- seq(from=0, to=xmax, length.out = 5E2)
  ys <- sapply(X=xs, FUN = function(z){calculate(x,'failure rate',z)})
  
  if(theme == 'base_r'){
    plot(xs, ys, type = 'l', ...)  
  }
  if(theme == 'weibull++'){
    weibull_theme_plot(xs, ys, 'Failure Rate vs. Time', 'Time', 'Failure Rate')
  }
}
plot_weibull <- function(x, theme, line_par, point_par, ...){
  #only need 2 points for straight line
  xmax <- 10^(ceiling(log10(max(x$data[[1]]))))
  xmin <- 10^(floor(log10(min(x$data[[1]]))))
  xs <- c(xmin,xmax)
  ys <- pweibull(xs, shape = x$fit$shape, scale = x$fit$scale)  
  #make sure xmax isn't so large that ys = 1. this will cause inifinities
  if(ys[2] == 1){
    xs[2] <- qweibull(0.999, shape = x$fit$shape, scale = x$fit$scale)
    ys[2] <- 0.999
  }
  
  #y tick locations
  ymin <- (floor(log10(ys[1])))
  ymax <- (ceiling(log10(ys[2])))
  if(ymax == 0) ymax = -1
  y_lbs <- sort(c(10^(ymin:ymax), 1-(10^(ymin:ymax))))
  y_ats <- log10(-log(1-y_lbs))
  y_lbs <- y_lbs[!is.infinite(y_ats)]
  y_ats <- y_ats[!is.infinite(y_ats)]
  y_mini_lbs <- 10^unlist(lapply(X = y_lbs[1:(length(y_lbs)/2)], 
                              FUN = function(y) log10(y)+log10(2:9)))
  at_end <- F
  for(i in ((length(y_lbs)/2)+1):length(y_lbs)){
    if(i == length(y_lbs)){
      y_lbs[i+1] <- as.numeric(paste0(y_lbs[i],'9'))
      at_end <- T
    }
    to_add <- seq(from = y_lbs[i], to = y_lbs[i+1], by = diff(c(y_lbs[i],y_lbs[i+1]))/9)
    to_add <- to_add[-c(1,length(to_add))]
    y_mini_lbs <- c(y_mini_lbs,to_add)
    if(at_end){
      y_lbs <- y_lbs[-length(y_lbs)]
    }
  }
  y_mini_ats <- log10(-log(1-y_mini_lbs))
  
  #linearize data 
  xs <- log10(xs)
  ys <- log10(-log(1-ys))
  x_ats <- xs[1]:xs[2]
  x_lbs <- 10^x_ats
  x_mini_ats <- unlist(lapply(X = x_ats, FUN = function(y) y+log10(2:9)))
  x_mini_lbs <- 10^x_mini_ats
  
  #get points to plot
  pts_ys <- log10(-log(1-median_ranks(x$data)))
  pts_xs <- log10(x$data$time[x$data$status == 1])
  
  if(theme == 'base_r'){
    do.call(plot, c(list(x = xs, y = ys, type = 'l', xaxt = 'n', yaxt = 'n'), line_par, ...))
    # plot(xs, ys, type = 'l', xaxt = 'n', yaxt = 'n', ...)
    axis(side = 1, at = x_ats, labels = x_lbs)
    axis(side = 2, at = y_ats, labels = y_lbs)
    do.call(points, c(list(x = pts_xs, y = pts_ys), point_par, ...))
    # points(pts_xs, pts_ys, ...)
  }
  if(theme == 'weibull++'){
    par(las = 1, xaxs = 'i', yaxs = 'i', tcl = 0, mar = c(1.6, 2.2, 1.6, 2.1), cex.axis = 0.5, 
        mgp = c(0.75,0,0), col.axis = 'blue', col.lab = 'red', col.main = 'red', cex.main = 1)
    plot(xs, ys, type = 'l', xaxt = 'n', yaxt = 'n', ann = F, col='blue')
    axis(side = 1, at = x_ats, labels = x_lbs, lwd = 0, line = -0.4)
    axis(side = 2, at = y_ats, labels = y_lbs, lwd = 0, line = 0.1)
    title(ylab='Unreliability', line = 1.25)
    title(xlab='Time', line = 0.5)
    title(main='Probability - Weibull', line =0.5)
    abline(v = x_mini_ats, col = 'green',lwd= 0.5)
    abline(h = y_mini_ats, col = 'green',lwd= 0.5)
    abline(v = x_ats, col = 'red', lwd = 0.5)
    abline(h = y_ats, col = 'red', lwd = 0.5)
    lines(xs, ys, col='blue')
    points(pts_xs, pts_ys, pch = 16, col = 'blue')
    par(las = 0, xaxs = 'r', yaxs = 'r', tcl = -0.5, mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, 
        mgp = c(3,1,0), col.axis = 'black', col.lab = 'black', col.main='black', cex.main = 1.2)
  }
}
plot_exponential <- function(x, theme, line_par, point_par, ...){
  #get points to plot - needed for y-scaling
  pts_ys <- log(1-median_ranks(x$data))
  pts_xs <- x$data$time[x$data$status == 1]
  
  #find x-range of plot
  x_ats <- pretty(x$data$time)
  if(min(x_ats) > min(x$data$time)){
    x_ats <- c(min(x_ats) - diff(x_ats)[1] , x_ats)
  }
  if(max(x_ats) < max(x$data$time)){
    x_ats <- c(x_ats, max(x_ats) + diff(x_ats)[1])
  }
  x_mini_ats <- pretty(x_ats[1:2])
  x_mini_ats <- seq(from = x_mini_ats[1], to = x_ats[length(x_ats)],by = diff(x_mini_ats)[1])
  x_mini_ats <- x_mini_ats[!(x_mini_ats %in% x_ats)]
  
  #only need 2 points for straight line
  xs <- c( x_ats[1], x_ats[length(x_ats)] )
  rs <- 1-pexp(xs, rate = 1/x$fit$scale)
  ys <- log(rs)
  
  #y_scale
  pts_y_range <- range(exp(pts_ys))
  ymin <- floor(log10(min(c(rs[2],pts_y_range[1]))))
  ymax <- ceiling(log10(max(c(rs[1],pts_y_range[2]))))
  # if(ymax == 0) ymax = -1
  y_lbs <- c(10^(ymin:ymax))
  y_ats <- log(y_lbs)
  y_mini_lbs <- 10^unlist(lapply(X = y_lbs[-length(y_lbs)], FUN = function(y) log10(y)+log10(2:9)))
  y_mini_ats <- log(y_mini_lbs)
  
  
  if(theme == 'base_r'){
    do.call(plot, c(list(x = xs, y = ys, type = 'l', xaxt = 'n', yaxt = 'n', ylim = c(min(y_ats), max(y_ats))),
                    line_par, ...))
    # plot(xs, ys, type = 'l', xaxt = 'n', yaxt = 'n', ylim = c(min(y_ats),max(y_ats)), ...)
    axis(side = 1, at = x_ats)
    axis(side = 2, at = y_ats, labels = y_lbs)
    axis(side = 2, at = y_mini_ats, labels = F)
    do.call(points, c(list(x = pts_xs, y = pts_ys), point_par, ...))
    # points(pts_xs, pts_ys, ...)
  }
  if(theme == 'weibull++'){
    par(las = 1, xaxs = 'i', yaxs = 'i', tcl = 0, mar = c(1.6, 2.2, 1.6, 2.1), cex.axis = 0.5, 
        mgp = c(0.75,0,0), col.axis = 'blue', col.lab = 'red', col.main = 'red', cex.main = 1)
    plot(xs, ys, type = 'l', xaxt = 'n', yaxt = 'n', ann = F, col='blue',
         ylim = c(min(y_ats),max(y_ats)))
    axis(side = 1, at = x_ats, lwd = 0, line = -0.4)
    axis(side = 2, at = y_ats, labels = y_lbs, lwd = 0, line = 0.1)
    title(ylab='Reliability', line = 1.25)
    title(xlab='Time', line = 0.5)
    title(main='Probability - Exponential', line =0.5)
    abline(v = x_mini_ats, col = 'green',lwd= 0.5)
    abline(h = y_mini_ats, col = 'green',lwd= 0.5)
    abline(v = x_ats, col = 'red', lwd = 0.5)
    abline(h = y_ats, col = 'red', lwd = 0.5)
    lines(xs, ys, col='blue')
    points(pts_xs, pts_ys, pch = 16, col = 'blue', xpd = T)
    par(las = 0, xaxs = 'r', yaxs = 'r', tcl = -0.5, mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, 
        mgp = c(3,1,0), col.axis = 'black', col.lab = 'black', col.main='black', cex.main = 1.2)
  }
}

weibull_theme_plot <- function(x,y,lab1, lab2, lab3){
  par(las = 1, xaxs = 'i', yaxs = 'i', tcl = 0, mar = c(1.6, 2.2, 1.6, 2.1), cex.axis = 0.5, 
      mgp = c(0.75,0,0), col.axis = 'blue', col.lab = 'red', col.main = 'red', cex.main = 1)
  xmarks <- pretty(x)
  xminis <- seq(from = xmarks[1], to = xmarks[length(xmarks)], by = diff(pretty(xmarks[1:2]))[1])
  ymarks <- pretty(y)
  yminis <- seq(from = ymarks[1], to = ymarks[length(ymarks)], by = diff(pretty(ymarks[1:2]))[1])
  plot(x, y, type = 'l', main = '', xlab = "", ylab = "",
       xlim = c(min(xmarks),max(xmarks)), ylim = c(min(ymarks),max(ymarks)),yaxt = 'n', xaxt = 'n',
       col = 'blue') 
  axis(side = 2, at = ymarks, tick=F, line = 0.1)
  axis(side = 1, at = xmarks, tick=F, line = -0.4)
  title(ylab=lab3, line = 1.25)
  title(xlab=lab2, line = 0.5)
  title(main=lab1, line =0.5)
  abline(v = xminis, col = 'green',lwd= 0.5)
  abline(h = yminis, col = 'green',lwd= 0.5)
  abline(v = xmarks, col = 'red', lwd = 0.5)
  abline(h = ymarks, col = 'red', lwd = 0.5)
  lines(x, y, col='blue')
  par(las = 0, xaxs = 'r', yaxs = 'r', tcl = -0.5, mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, 
      mgp = c(3,1,0), col.axis = 'black', col.lab = 'black', col.main='black', cex.main = 1.2)
}

hist.fitted_life_data <- function(x, type = 'failure', theme = 'base_r'){
  times <- x$data$time
  status <- x$data$status
  hist_lifedata(times, status, type, theme) 
}
hist.lifedata <- function(x, type = 'failure', theme = 'base_r'){
  times <- x$time
  status <- x$status
  hist_lifedata(times, status, type, theme)
}
hist_lifedata <- function(x, status, type, theme){
  if(!(type %in% c('failure','suspension'))){
    stop('type must be either failure or suspension')
  }else{
    if(type == 'failure'){
      times <- x[status == 1]
    }else{
      times <- x[status == 0]
    }
  }
  
  if(!(theme %in% c('base_r'))){
    warning('theme unknown. defaulting to base_r')
    theme <- 'base_r'
  }
  if(theme == 'base_r'){
    hist(times)  
  }
}

#' @title Calculate Metrics
#' 
#' @description 
#' generic function for fitted_life_data only at this time
#' 
#' @details
#' see calculate.fitted_life_data
#' 
#' @param x only defined for a fitted_life_data object
#' @return 
#' returns numeric value
calculate <- function(x, ...) UseMethod('calculate')
#' @title Calculate Metrics From Fitted Distribution
#' 
#' @description 
#' calculate various values based on input and fitted distribution
#' 
#' @details
#' Values for value input
#' \tabular{ll}{
#' reliability  \tab reliability value at input. Prob(T>input)\cr
#' failure \tab probability of failure by input. Prob(T<input)\cr
#' mean life \tab mean of fitted distribution\cr
#' failure rate \tab failure/hazard at input. Prob(T=input+delta|T>input) \cr
#' }
#' 
#' @param x a fitted_life_data object
#' @param value desired metric. see details
#' @param input numeric input for certain values. ignored if not needed.
#' @return 
#' desired numeric value
calculate.fitted_life_data <- function(x, value, input  = NA, cond_input = NA){
  if(!(value %in% c('reliability', 'failure', 'mean life', 'failure rate', 'reliable life', 'bx life',
                    'cond reliab', 'cond fail'))){
    stop(paste('value',value,'not recognized'))
    to_return <- NULL
  }else{
    #probability of reliability and failure -------------------
    if(value %in% c('reliability', 'failure')){
      if(value == 'reliability'){
        use_lower_tail = F
      }else{
        use_lower_tail = T
      } 
      
      if(x$dist == 'exponential'){
        to_return <- pexp(input, 1/x$fit$scale, lower.tail = use_lower_tail)
      }
      if(x$dist == 'weibull'){
        to_return <- pweibull(input, x$fit$shape, x$fit$scale, lower.tail = use_lower_tail)
      }
    }
    #mean life ---------------
    if(value == 'mean life'){
      if(x$dist == 'exponential'){
        to_return <- x$fit$scale
      }
      if(x$dist == 'weibull'){
        to_return <- x$fit$scale * gamma(1+(1/x$fit$shape))
      }  
    }
    #failure rate ------------
    if(value == 'failure rate'){
      if(x$dist == 'exponential'){
        to_return <- 1/x$fit$scale 
      }
      if(x$dist == 'weibull'){
        to_return <- x$fit$shape * (input^(x$fit$shape-1)) * ((1/x$fit$scale)^x$fit$shape)
      }       
    }
    #warranty life --------
    if(value == 'reliable life'){
      if(x$dist == 'exponential'){
        to_return <- qexp(p = input, rate = 1/x$fit$scale, lower.tail = F)
      }
      if(x$dist == 'weibull'){
        to_return <- qweibull(p = input, shape = x$fit$shape, scale = x$fit$scale, lower.tail = F)
      }
    }
    #warranty life --------
    if(value == 'bx life'){
      if(x$dist == 'exponential'){
        to_return <- qexp(p = input, rate = 1/x$fit$scale)
      }
      if(x$dist == 'weibull'){
        to_return <- qweibull(p = input, shape = x$fit$shape, scale = x$fit$scale)
      }
    }
    #conditional reliability --------
    if(value == 'cond reliab'){
      if(x$dist == 'exponential'){
        to_return <- pexp(q = input, rate = 1/x$fit$scale, lower.tail = F)
      }
      if(x$dist == 'weibull'){
        to_return <- pweibull(q = input+cond_input, shape = x$fit$shape, scale = x$fit$scale, lower.tail = F) / 
          pweibull(q = cond_input, shape = x$fit$shape, scale = x$fit$scale, lower.tail = F)
      }
    }
    #conditional failure --------
    if(value == 'cond fail'){
      if(x$dist == 'exponential'){
        to_return <- pexp(q = input, rate = 1/x$fit$scale)
      }
      if(x$dist == 'weibull'){
        num <- pweibull(q = input+cond_input, shape = x$fit$shape, scale = x$fit$scale)  - 
          pweibull(q = cond_input, shape = x$fit$shape, scale = x$fit$scale)
        to_return <- num / pweibull(q = cond_input, shape = x$fit$shape, scale = x$fit$scale, lower.tail = F)
      }
    }
  }
  
  return(to_return)
}







