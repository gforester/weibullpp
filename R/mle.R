



#use Nelder-Mead default in R
#assume the_data a Surv object
find_mle <- function(the_data){
  
  func_for_optim <- function(th) -1*weibull_2p_like_1(th, the_data)
  guesses <- c(1,mean(the_data[,1]))
  res <- optim(guesses, func_for_optim, control = list(reltol = 1E-10, maxit = 1E3))
  
  mle = list(shape = res$par[1], scale = res$par[2], log_like = -1*res$value)
  return(mle)
  
}

#unbiased estimator -> see weibull

