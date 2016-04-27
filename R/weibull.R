
#brute-force R method for log-likelihood
weibull_ld_like <- function(theta, lifedata_object){
  data_full <- lifedata_object[[1]][lifedata_object[[2]] == 1]
  data_cens <- lifedata_object[[1]][lifedata_object[[2]] == 0]
  full_like <- cens_like <- 0
  
  if(length(data_full)>0){
    full_like <- sum(sapply(X = data_full, FUN = dweibull, shape = theta[1], scale = theta[2], log = T))  
  }
  if(length(data_cens)>0){
    cens_like <- sum(sapply(X = data_cens, FUN = pweibull, shape = theta[1], scale = theta[2], lower.tail=F, log.p = T))
  }
  
  return(full_like + cens_like)
}

#use Nelder-Mead default in R
#assume the_data a lifedata object
weibull_mle <- function(the_data){
  
  func_for_optim <- function(th) -1*weibull_ld_like(th, the_data)
  guesses <- c(1,mean(the_data[[1]]))
  res <- optim(guesses, func_for_optim, hessian = T, control = list(reltol = 1E-10, maxit = 1E3))
  
  mle = list(shape = res$par[1], scale = res$par[2], log_like = -1*res$value, hessian = -1*res$hessian)
  return(mle)
  
}

weibull_fisher_ci <- function(mle_res){
  neg_l <- -1*mle_res$hessian
  fisher <- solve(neg_l)
  return(list(shape_delta = sqrt(fisher[1,1]), scale_delta = sqrt(fisher[2,2])))
}

weibull_hessian <- function(x, shape, scale){
  return(-1*optimHess(c(shape, scale), weibull_ld_like, lifedata_object = x))
}

#rank regression
weibull_rr <- function(x){
  
  #get median ranks estimate of probability of failure
  ts <- x$time[x$status == 1]
  ps <- median_ranks(x)
  
  #linearize
  ys <- log10(-log(1-ps))
  xs <- log10(ts)
  
  #est slope-intercept (b-a)
  rho <- cor(xs,ys)
  sdx <- sd(xs)
  sdy <- sd(ys)
  b_y <- rho*sdy/sdx
  b_x <- rho*sdx/sdy
  a_y <- mean(ys)-b_y*mean(xs)
  a_x <- mean(xs)-b_x*mean(ys)
  
  #est shape and scale
  shape_rry <- b_y
  shape_rrx <- 1/b_x
  scale_rry <- 10^(-a_y/b_y) 
  scale_rrx <- 10^(a_x/(b_x*shape_rrx))
    
  return(list(scale_rry = scale_rry, scale_rrx = scale_rrx, 
              shape_rry = shape_rry, shape_rrx = shape_rrx,
              log_like_rry = exp_log_like(x, scale_rry),
              log_like_rrx = exp_log_like(x, scale_rrx),
              rho = rho))
  
}
