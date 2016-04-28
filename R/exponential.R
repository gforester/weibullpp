
exp_log_like <- function(x, scale){
  full_like <- 0
  cens_like <- 0
  if(any(x[[2]]==1)){
    full_like <- sum(dexp(x[[1]][x[[2]]==1],rate = 1/scale, log = T))  
  }
  if(any(x[[2]]==0)){
    cens_like <- sum(pexp(x[[1]][x[[2]]==0],rate = 1/scale, lower.tail=F, log.p = T))  
  }
  return(full_like + cens_like)
}


#right censored only
fit_exp <- function(x){
  mle_scale <- sum(x[[1]]/sum(x[[2]]))
  # full_like <- sum(dexp(x[[1]][x[[2]]==1],rate = 1/mle_scale, log = T))
  # cens_like <- sum(pexp(x[[1]][x[[2]]==0],rate = 1/mle_scale, lower.tail=F, log.p = T))
  l <- exp_log_like(x, mle_scale)
  #in future, replace above likelihood calculation with analytic expression
  return(list(scale = mle_scale,
              # log_like = full_like + cens_like,
              log_like = l,
              se = exp_se(mle_scale, sum(x[[2]]))
              ))
}

#Fisher se
exp_se <- function(mle,N){
  return(list(fisher = sqrt((mle^2)/N)))
}

#rank regression
exp_rr <- function(x){
  #get median ranks estimate of probability of failure
  xs <- x$time[x$status == 1]
  ps <- median_ranks(x)
  
  #linearize
  ys <- log(1-ps)
  
  #est slopes
  rho <- cor(xs,ys)
  sdx <- sd(xs)
  sdy <- sd(ys)
  b_y <- rho*sdy/sdx
  b_x <- rho*sdx/sdy
  
  #est scale
  scale_rry <- -1/b_y
  scale_rrx <- -b_x
  
  return(list(scale_rry = scale_rry, scale_rrx = scale_rrx, 
              log_like_rry = exp_log_like(x, scale_rry),
              log_like_rrx = exp_log_like(x, scale_rrx),
              rho = rho))
}

#calculations
exp_calc <- function(fit, type, input = NA, cond_input = NA, alpha = 0.05){
  to_return <- 'type not available for distribution'
  
  #scale breadth
  z_star <- abs(qnorm(alpha/2))
  scale_upper <- fit$fit$scale + z_star*fit$std_err
  scale_lower <- fit$fit$scale - z_star*fit$std_err
  
  # reliabililty and probability of failure
  if(type %in% c('reliability','failure')){
    lower_tail <- switch(type,'reliability' = F, 'failure' = T)
    to_return <- pexp(input, 1/fit$fit$scale, lower.tail = lower_tail)
    limits <-pexp(input, 1/c(scale_upper,scale_lower), lower.tail = lower_tail)
    upper <- max(limits)
    lower <- min(limits)
  }
  # mean time to failure
  if(type == 'mean life'){
    to_return <- fit$fit$scale
    upper <- scale_upper
    lower <- scale_lower
  }
  # failure rate aka hazard function
  if(type == 'failure rate'){
    to_return <- 1/fit$fit$scale 
    upper <- 1/scale_lower
    lower <- 1/scale_upper
  }
  # reliable life aka warranty life and BX% life
  if(type %in% c('reliable life','bx life')){
    lower_tail <- switch(type,'reliable life' = F, 'bx life' = T)
    to_return <- qexp(p = input, rate = 1/fit$fit$scale, lower.tail = lower_tail)
    limits <- qexp(p = input, rate = 1/c(scale_upper, scale_lower), lower.tail = lower_tail)
    upper <- max(limits)
    lower <- min(limits)
  }
  # cond distribution 
  if(type %in% c('cond reliab','cond fail')){
    lower_tail <- switch(type,'cond reliab' = F, 'cond fail' = T)  
    to_return <- pexp(q = input, rate = 1/fit$fit$scale, lower.tail = lower_tail)
    limits <- pexp(q = input, rate = 1/c(scale_upper, scale_lower), lower.tail = lower_tail)
    upper <- max(limits)
    lower <- min(limits)
  }
  
  return(list(value = to_return, upper_limit = upper, lower_limit = lower))
  
}


