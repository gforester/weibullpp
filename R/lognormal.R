

# for complete data only ------------------------------
lognormal_ld_like_complete <- function(theta, lifedata_object){
  log_t <- log(lifedata_object[[1]])
  Y <- sum(log_t)
  X <- sum(log_t^2)
  N <- length(log_t)
  full_like <- -(0.5*N*log(2*pi)) - (N*log(theta[2])) - Y - (0.5*X/(theta[2]^2)) + 
    (theta[1]*Y/(theta[2]^2)) - ((N*theta[1]^2)/(2*theta[2]^2))
  return(full_like)
}

# for data with censored points -----------------------
#brute-force R method for log-likelihood - use for censored data
lognormal_ld_like <- function(theta, lifedata_object){
  data_full <- lifedata_object[[1]][lifedata_object[[2]] == 1]
  data_cens <- lifedata_object[[1]][lifedata_object[[2]] == 0]
  full_like <- cens_like <- 0
  
  if(length(data_full)>0){
    full_like <- sum(sapply(X = data_full, FUN = dlnorm, meanlog = theta[1], sdlog = theta[2], log = T))  
  }
  if(length(data_cens)>0){
    cens_like <- sum(sapply(X = data_cens, FUN = plnorm, meanlog = theta[1], sdlog = theta[2], lower.tail=F, log.p = T))
  }
  
  return(full_like + cens_like)
}

# MLE --------------
#use Nelder-Mead default in R
#assume the_data a lifedata object
lognormal_mle <- function(the_data){
  log_t <- log(the_data[[1]])
  if(sum(the_data[[2]]) == length(the_data[[2]])){
    #complete - solution analytical
    Y <- sum(log_t)
    X <- sum(log_t^2)
    N <- length(log_t)
    theta <- c(Y/N, sqrt((X/N)-((Y/N)^2)))
    mle = list(logmean = theta[1], sdlog = theta[2],log_like = lognormal_ld_like_complete(theta,the_data),
               hessian = optimHess(theta, lognormal_ld_like_complete,lifedata_object = the_data))
  }else{
    func_for_optim <- function(th) -1*lognormal_ld_like(th, the_data)
    guesses <- c(mean(log_t), sd(log_t))
    res <- optim(guesses, func_for_optim, hessian = T, control = list(reltol = 1E-10, maxit = 1E3))
    mle = list(logmean = res$par[1], sdlog = res$par[2], log_like = -1*res$value, hessian = -1*res$hessian)
  }

  return(mle)
  
}

lognormal_fisher_ci <- function(mle_res){
  neg_l <- -1*mle_res$hessian
  fisher <- solve(neg_l)
  return(list(logmean_delta = sqrt(fisher[1,1]), sdlog_delta = sqrt(fisher[2,2])))
}

# calculations -----------------

lognormal_calc <- function(fit, type, input = NA, cond_input = NA, alpha = 0.05){
  to_return <- 'type not available for distribution'
  
  #mapping of parameter estimates and convariance matrix
  logmean <- fit$fit$logmean
  sdlog <- fit$fit$sdlog
  cov_mat <- solve(-1*fit$fit$hessian)
  
  #critical value
  z_star <- abs(qnorm(alpha/2))
  
  # reliabililty and probability of failure
  if(type %in% c('reliability','failure')){
    lower_tail <- switch(type,'reliability' = F, 'failure' = T)
    to_return <- plnorm(input, logmean, sdlog, lower.tail = lower_tail)
    z <- (log(input)-logmean)/sdlog
    var_z <- ( cov_mat[1,1] + ((z^2)*cov_mat[2,2]) + (2*z*cov_mat[1,2]) )/(sdlog^2)
    limits <- pnorm(q = z+c(-1,1)*z_star*sqrt(var_z),lower.tail = lower_tail)
    upper <- max(limits)
    lower <- min(limits)
  }
  # mean time to failure
  if(type == 'mean life'){
    to_return <- exp(logmean+0.5*sdlog^2)
    delta <- to_return * sqrt(cov_mat[1,1] + ((sdlog^2) * cov_mat[2,2]) )
    limits <- to_return + c(1,-1)*z_star*delta
    upper <- max(limits)
    lower <- min(limits)
  }
  # # failure rate aka hazard function
  if(type == 'failure rate'){
    f <- dlnorm(input,logmean,sdlog)
    R <- plnorm(input, logmean, sdlog, lower.tail = F)
    to_return <- f/R
    z <- (log(input)-logmean)/sdlog
    var_z <- ( cov_mat[1,1] + ((z^2)*cov_mat[2,2]) + (2*z*cov_mat[1,2]) )/(sdlog^2)
    var_R <- var_z * exp(-z^2) / (2*pi)
    var_f <- var_z * z^2 * f^2
    var_h <- (to_return^2)*((var_f/f^2)+(var_R/R^2))
    # NOTE: this is a great simplication - covariance term removed
    upper <- to_return+z_star*sqrt(var_h)
    lower <- to_return-z_star*sqrt(var_h)
  }
  # reliable life aka warranty life and BX% life
  if(type %in% c('reliable life','bx life')){
    lower_tail <- switch(type,'reliable life' = F, 'bx life' = T)
    to_return <- qlnorm(input, logmean, sdlog, lower.tail = lower_tail)
    z <- qnorm(input, lower.tail = lower_tail)
    var_t <- cov_mat[1,1] + (z^2)*cov_mat[2,2] + 2*z*cov_mat[1,2]
    limits <- exp(log(to_return)+(c(-1,1)*z_star*sqrt(var_t)))
    upper <- max(limits)
    lower <- min(limits)
  }
  # # cond distribution
  if(type %in% c('cond reliab','cond fail')){
    lower_tail <- switch(type,'cond reliab' = F, 'cond fail' = T)
    num <- plnorm(input+cond_input, logmean, sdlog, lower.tail = lower_tail)
    den <- plnorm(cond_input, logmean, sdlog,lower.tail = F)
    to_return <-  num / den 
    z <- (log(cond_input+input)-logmean)/sdlog
    var_z <- ( cov_mat[1,1] + ((z^2)*cov_mat[2,2]) + (2*z*cov_mat[1,2]) )/(sdlog^2)
    var_num <- var_z * exp(-z^2) / (2*pi)
    z <- (log(cond_input)-logmean)/sdlog
    var_z <- ( cov_mat[1,1] + ((z^2)*cov_mat[2,2]) + (2*z*cov_mat[1,2]) )/(sdlog^2)
    var_den <- var_z * exp(-z^2) / (2*pi)
    var_cond <-  (to_return^2)*((var_num/num^2)+(var_den/den^2))
  #   limits <- pexp(q = input, rate = 1/c(scale_upper, scale_lower), lower.tail = lower_tail)
    # NOTE: this is a great simplication - covariance term removed
    upper <- to_return+z_star*sqrt(var_cond)
    lower <- to_return-z_star*sqrt(var_cond)
  }
  
  return(list(value = to_return, upper_limit = upper, lower_limit = lower))
  
}

# linearized plot --------------
# x is log base e of times
# y is inverse standard normal CDF of F(t) estimates
# axes will be marked as in log base 10
plot_lognormal <- function(x, theme, alpha, line_par, point_par, ...){
  #get points to plot - needed for y-scaling
  pts_ys <- qnorm((median_ranks(x$data)))
  pts_xs <- log(x$data$time[x$data$status == 1])

  # find x-range of plot
  # lowest of powers of 10 down.
  # weibull++ is set to 0, this is weird?
  min_power <- floor(log10(min(x$data$time)))
  max_power <- ceiling(log10(max(x$data$time)))
  x_lbs <- 10^(min_power:max_power)
  x_ats <- log( x_lbs )
  x_mini_lbs <- sapply(X = x_lbs, FUN = function(y) y*2:9)
  x_mini_ats <- log(x_mini_lbs)

  #calculate a-b line of fitted distribution
  slope <- 1 / x$fit$sdlog
  intercept <- -x$fit$logmean / x$fit$sdlog 
  
  #y_scale
  # cat(median_ranks(x$data),'\n')
  # cat(pts_ys, '\n')
  pts_y_cdf_range <- range(pnorm(pts_ys))
  ymin <- floor(log10(pts_y_cdf_range[1]))
  ymax <- ceiling(log10(pts_y_cdf_range[2]))
  if(ymax == 0) ymax = -1
  y_lbs <- c(10^(ymin:ymax))
  y_lbs <- c(y_lbs, 1 - y_lbs) #max should be at 1-lowest tick mark
  y_ats <- qnorm(y_lbs)
#   y_mini_lbs <- 10^unlist(lapply(X = y_lbs[-length(y_lbs)], FUN = function(y) log10(y)+log10(2:9)))
#   y_mini_ats <- log(y_mini_lbs)
  
  
  if(theme == 'base_r'){
    to_add = list(xlab = 'Time', ylab = 'Probability')
    to_add$main = 'Unreliability'
    to_add$ylim = range(y_ats)
    to_add$xlim = range(x_ats)
    if(any(log(x$data$time) == to_add$xlim[2])){
      to_add$xlim[2] = to_add$xlim[2] + 0.04*diff(range(x_ats))
    }
    more_pars = list(...)
    if(any(names(list(...)) == 'xlab')){
      to_add$xlab = NULL
    }
    if(any(names(list(...)) == 'ylab')){
      to_add$ylab = NULL
    }
    if(any(names(list(...)) == 'main')){
      to_add$main = NULL
    }
    if(any(names(list(...)) == 'xlim')){
      to_add$xlim = NULL
      more_pars$xlim = log(more_pars$xlim)
    }
    if(any(names(list(...)) == 'ylim')){
      to_add$ylim = NULL
    }
    if(length(to_add) == 0){
      to_add = NULL
    }
    do.call(plot, c(list(x = pts_xs, y = pts_ys, type = 'l', xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i'),
                    line_par, to_add, more_pars))
    do.call(abline,c(list(a = intercept, b = slope), line_par))
    # bound_xs <- seq(from=par()$usr[1], to=par()$usr[2], length.out = 5E2)
    # calc_res <- lapply(X = bound_xs, FUN = function(z){calculate(x,'reliability',z,alpha=alpha)})
    # uppers <- sapply(X = calc_res, FUN = function(z){log(z$upper_limit)})
    # lowers <- sapply(X = calc_res, FUN = function(z){log(z$lower_limit)})
    # lines(x=bound_xs, y=uppers, lty = 2)
    # lines(x=bound_xs, y=lowers, lty = 2)
    axis(side = 1, at = x_ats, labels = x_lbs)
    axis(side = 2, at = y_ats, labels = y_lbs)
    # axis(side = 2, at = y_mini_ats, labels = F)
    do.call(points, c(list(x = pts_xs, y = pts_ys), point_par, more_pars))
  }
#   if(theme == 'ggplot'){
#     #point styling
#     if(any(names(point_par) == 'col')){
#       point_colour <- point_par$col
#     }else{
#       point_colour <- 'black'
#     }
#     if(any(names(point_par) == 'fill')){
#       point_fill <- point_par$col
#     }else{
#       point_fill <- 'black'
#     }
#     if(any(names(point_par) == 'cex')){
#       point_size <- point_par$cex
#     }else{
#       point_size <- 1.5
#     }
#     #line styling
#     if(any(names(line_par) == 'lwd')){
#       line_size <- line_par$lwd
#     }else{
#       line_size <- 0.5
#     }
#     if(any(names(line_par) == 'col')){
#       line_colour <- line_par$col
#     }else{
#       line_colour <- 'black'
#     }
#     
#     to_plot <- ggplot()+
#       geom_abline(slope=slope, intercept=intercept, size = line_size, colour = line_colour)+
#       geom_point(data=data.frame(x=pts_xs,y=pts_ys),mapping=aes(x=x,y=y),
#                  fill = point_fill, colour = point_colour, size = point_size)+
#       scale_y_continuous(breaks=y_ats, minor_breaks = y_mini_ats, labels = y_lbs)+
#       scale_x_continuous(breaks=x_ats)
#     to_add <- character(0)
#     par_list <- list(...)
#     if(any(names(par_list) == 'xlab')){
#       to_add <- c(to_add,paste0('xlab("',par_list$xlab,'")'))
#     }
#     if(any(names(par_list) == 'ylab')){
#       to_add <- c(to_add,paste0('ylab("',par_list$ylab,'")'))
#     }
#     if(length(par_list)>0){
#       eval(parse(text = paste0('to_plot <- to_plot + ',paste0(to_add,collapse=" + "))))
#     }
#     bound_xs <- seq(from=par()$usr[1], to=par()$usr[2], length.out = 5E2)
#     calc_res <- lapply(X = bound_xs, FUN = function(z){calculate(x,'reliability',z,alpha=alpha)})
#     uppers <- sapply(X = calc_res, FUN = function(z){log(z$upper_limit)})
#     lowers <- sapply(X = calc_res, FUN = function(z){log(z$lower_limit)})
#     to_plot <- to_plot+geom_line(data=data.frame(x=bound_xs,y=uppers),mapping=aes(x=x,y=y),linetype=2)+
#       geom_line(data=data.frame(x=bound_xs,y=lowers),mapping=aes(x=x,y=y),linetype=2)
#     print(to_plot)
#   }
  if(theme == 'weibull++'){
    par(las = 1, xaxs = 'i', yaxs = 'i', tcl = 0, mar = c(1.6, 2.2, 1.6, 2.1), cex.axis = 0.5, 
        mgp = c(0.75,0,0), col.axis = 'blue', col.lab = 'red', col.main = 'red', cex.main = 1)
    to_add_xlim = range(x_ats)
    if(any(log(x$data$time) == to_add_xlim[2])){
      to_add_xlim[2] = to_add_xlim[2] + 0.04*diff(range(x_ats))
    }
    plot(pts_xs, pts_ys, xaxt = 'n', yaxt = 'n', ann = F, col='blue',pch = 16,xlim = to_add_xlim,
         ylim = c(min(y_ats),max(y_ats)))
    axis(side = 1, at = x_ats, labels = x_lbs, lwd = 0, line = -0.4)
    axis(side = 2, at = y_ats, labels = y_lbs, lwd = 0, line = 0.1)
    title(ylab='Unreliability', line = 1.25)
    title(xlab='Time', line = 0.5)
    title(main='Probability - Lognormal', line =0.5)
    abline(v = x_mini_ats, col = 'green',lwd= 0.5)
#     abline(h = y_mini_ats, col = 'green',lwd= 0.5)
    abline(v = x_ats, col = 'red', lwd = 0.5)
    abline(h = y_ats, col = 'red', lwd = 0.5)
    abline(intercept, slope, col='blue')
#     bound_xs <- seq(from=par()$usr[1], to=par()$usr[2], length.out = 5E2)
#     calc_res <- lapply(X = bound_xs, FUN = function(z){calculate(x,'reliability',z,alpha=alpha)})
#     uppers <- sapply(X = calc_res, FUN = function(z){log(z$upper_limit)})
#     lowers <- sapply(X = calc_res, FUN = function(z){log(z$lower_limit)})
#     lines(x=bound_xs, y=uppers, col='red')
#     lines(x=bound_xs, y=lowers, col='red')
    # points(pts_xs, pts_ys, pch = 16, col = 'blue', xpd = T)
    par(las = 0, xaxs = 'r', yaxs = 'r', tcl = -0.5, mar = c(5.1, 4.1, 4.1, 2.1), cex.axis = 1, 
        mgp = c(3,1,0), col.axis = 'black', col.lab = 'black', col.main='black', cex.main = 1.2)
  }
}


