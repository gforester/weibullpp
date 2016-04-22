
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
  res <- optim(guesses, func_for_optim, control = list(reltol = 1E-10, maxit = 1E3))
  
  mle = list(shape = res$par[1], scale = res$par[2], log_like = -1*res$value)
  return(mle)
  
}
