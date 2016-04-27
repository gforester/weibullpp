
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


