



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

fit_data.lifedata <- function(x, dist = 'weibull'){
  if(!((dist == 'weibull') | (dist == 'exponential'))) stop('dist must be "weibull" or "exponential"')
  
  if(dist == 'exponential'){
    res <- fit_exp(x)
  }else{
    res <- weibull_mle(x)
  }
  
  to_return <- list(data = x, dist = dist, fit = res)
  class(to_return) <- 'fitted_life_data'
  
  return(to_return)
}



plot.fitted_life_data <- function(x, type = 'failure'){
  if(type == 'failure'){
    plot_cdfs(x, lower.tail = T)
  }
  if(type == 'reliability'){
    plot_cdfs(x, lower.tail = F)
  }
}
plot_cdfs <- function(x, lower.tail){
  #get Kaplan-Meier
  times <- Surv(x$data[[1]],x$data[[2]])
  fit <- survfit(times ~ 1)
  #finish later with Modified Kaplan Meier
  
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
  plot(xs, ys, type = 'l')
}
