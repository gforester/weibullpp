



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



plot.fitted_life_data <- function(x, type = 'failure', theme = 'base_r'){
  if(type == 'failure'){
    plot_cdfs(x, lower.tail = T, theme)
  }
  if(type == 'reliability'){
    plot_cdfs(x, lower.tail = F, theme)
  }
  if(type == 'pdf'){
    plot_pdf(x, theme)
  }
}
plot_cdfs <- function(x, lower.tail, theme){
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
plot_pdf <- function(x, theme){
  xmax <- 1.02 * max(x$data[[1]])
  xs <- seq(from=0, to=xmax, length.out = 5E2)
  if(x$dist == 'exponential'){
    ys <- dexp(xs, shape = x$fit$shape, scale = x$fit$scale)
  }
  if(x$dist == 'weibull'){
    ys <- dweibull(xs, shape = x$fit$shape, scale = x$fit$scale)
  }
  
  if(theme == 'base_r'){
    plot(xs, ys, type = 'l')  
  }
  if(theme == 'ggplot'){
    plot(xs, ys, type = 'l') 
  }
  if(theme == 'weibull++'){
    par(las = 1, xaxs = 'i', yaxs = 'i', tcl = 0, mar = c(2.1, 2.1, 2.1, 4.1), cex.axis = 0.5, mgp = c(1,0,0))
    xmarks <- pretty(xs)
    xminis <- seq(from = xmarks[1], to = xmarks[length(xmarks)], by = diff(pretty(xmarks[1:2]))[1])
    ymarks <- pretty(ys)
    yminis <- seq(from = ymarks[1], to = ymarks[length(ymarks)], by = diff(pretty(ymarks[1:2]))[1])
    plot(xs, ys, type = 'l', main = 'Probability Density Function', xlab = "Time", ylab = "f(t)",
         xlim = c(min(xmarks),max(xmarks)), ylim = c(min(ymarks),max(ymarks)),
         col = 'blue') 
    abline(v = xminis, col = 'green',lwd= 0.5)
    abline(h = yminis, col = 'green',lwd= 0.5)
    abline(v = xmarks, col = 'red', lwd = 0.5)
    abline(h = ymarks, col = 'red', lwd = 0.5)
    lines(xs, ys, col='blue')
    par(las = 0, xaxs = 'r', yaxs = 'r', tcl = -0.5, mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, mgp = c(3,1,0))
  }
  
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

calculate <- function(x, ...) UseMethod('calculate')
calculate.fitted_life_data <- function(x, value, input){
  if(!(value %in% c('reliability', 'failure', 'mean life', 'failure rate'))){
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

  }
  
  return(to_return)
}







