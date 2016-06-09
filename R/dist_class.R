#dist_class.R

#dist class structure: 
#  th - character vector containing names of params in order of input into functions 
#  d,p,r,g - functions: should be the typical d,p,q,r calculation structure except parameter inputs as one vector
#  g - function to determine guess for mle based on lifedata object input
#  lins - list of 2: function to linearize y given ranks and function to linearize x given times
dist_class <- function(th,d,p,q,r,g,lins=NA){
  if(is.na(lins)){
    lins = list('use'='no')
  }
  return(structure(list(thetas = th, d=d,p=p,q=q,r=r,guess = g, lins = lins), class="dist_class"))
}

#x is lifedata object
#y is dist_class object
#theta is parameters
log_like_generic <- function(x,y,theta){
  l_cens <- 0
  l_full <- 0
  tmp <- split(x[[1]],x[[2]])
  full <- tmp$`1`
  cens <- tmp$`0`
  if(!is.null(full)){
    l_full <- Reduce("+",y$d(full,theta, log = T),0)
  }
  if(!is.null(cens)){
    l_cens <- Reduce("+",y$p(cens,theta, lower.tail = F, log.p = T),0)
  }
  return(l_cens + l_full)
}

#mle generic
#x is lifedata object
#y is dist_class object
mle_generic <- function(x,y){
  to_optim <- function(theta){
    return(-1*log_like_generic(x,y,theta))
  }
  guesses <- y$guess(x)
  res <- optim(guesses, to_optim, hessian = T, control = list(reltol = 1E-10, maxit = 1E3))
  if(res$convergence != 0 ){
    code_in_english <- switch(as.character(res$convergence),
                              '1' = 'max iteration reached',
                              '10' = 'degeneracy of Nelder-Mead',
                              '51' = 'warning in L-BFGS-B',
                              '52' = 'error in L-BFGS-B',
                              paste('uknown flag', res$convergence))
    mle = list('error' = 'optimization error', 'code' = code_in_english, 'message' = res$message)
  }else{
    mle_params <- res$par
    names(mle_params) <- y$thetas
    mle = c(mle_params,list(log_like = -1*res$value, hessian = -1*res$hessian))
  }
  return(mle)
}

rr_generic <- function(x,y, method = 'median_ranks'){
  #get x,y for linearized plot
  ranks <- estimate_ranks(x,method)
  lin_ys <- y$lins[[2]](ranks)
  lin_xs <- y$lins[[1]](x$time[x$status == 1])
  
  #possible that ranks could result in non-finite numbers
  to_keep <- is.finite(lin_ys)
  lin_ys <- lin_ys[to_keep]
  lin_xs <- lin_xs[to_keep]
  
  #correlation, sds, and means
  rho <- cor(xs,ys)
  sdx <- sd(xs)
  sdy <- sd(ys)
  b_y <- rho*sdy/sdx
  b_x <- rho*sdx/sdy
  a_y <- mean(ys)-b_y*mean(xs)
  a_x <- mean(xs)-b_x*mean(ys)
  
  #connect to parameters
}

#y is a dist_class object
test_fit_data <- function(x, dist = y, method = 'mle', rank_method = 'median_ranks'){
  
}
