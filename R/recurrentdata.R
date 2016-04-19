
recurrentdata <- function(x, ids, status, units = NULL){
  if(!class(x) %in% c('numeric', 'integer')) stop('input must be numeric')
  new_recurrentdata <- structure(list(ID = ids, times = x, status = status), class='recurrentdata')
  attr(new_recurrentdata,'units') <- units
  return(new_recurrentdata)
}

plot.recurrentdata <- function(x){
  plot(c(x$times[1],x$times), 0:length(x$times),type = 'S', ylab = 'Num Occurences', xlab = 'Time')
}

calculate_v_age <- function(x,q, model_type = 1){
  if(model_type == 1){
    v_ages <- x$times[x$status == 1]*q
  }else{
    if(model_type == 2){
      xs <- diff(c(0,x$times[x$status == 1]))
      v_ages <- rep(0, sum(x$status))
      v_ages[1] <- q*xs[1]
      for(i in 2:length(v_ages)){
        v_ages[i] <- q*(xs[i]+v_ages[i-1])
      }
    }else{
      stop('model_type must either be 1 or 2')
    }
  }
  return(v_ages)
}

#note typo in Reliasoft documentation
grp_like <- function(theta, the_data){
  xs <- the_data$xs
  vs <- the_data$vs
  vn <- vs[length(vs)]
  n <- the_data$n
  deltaT <- the_data$deltaT
  
  beta <- theta$beta
  lambda <- theta$lambda
  
  stepped_vs <- c(0,vs[-length(vs)])
  l <- n*(log(lambda)+log(beta)) - 
    lambda*((deltaT+vn)^beta - vn^beta) - 
    lambda*sum((xs+stepped_vs)^beta-stepped_vs^beta) + 
    (beta-1)*sum(log(xs+stepped_vs))
  return(l)
}

mle_recurrent <- function(x, model_type = 1){
  if(length(unique(x$ID))>1) stop("Can't handle more than 1 ID at this time")
  the_data <- list(deltaT = max(x$times) - max(x$times[x$status == 1]),
                   n = sum(x$status),
                   xs = diff(c(0,x$times)))
  func_for_optim <- function(th) {
    q <- th[1]
    the_big_list <- c(the_data, list(vs = calculate_v_age(x,q,model_type)))
    to_return = -1*grp_like(list(beta = th[2], lambda = th[3]), the_big_list)
    if(is.na(to_return)) to_return <- rnorm(1,0,1)+100
    # print(paste(paste(th,collapse=" "),to_return))
    return(to_return)
  }
  guesses <- c(0.05,1.1,1/(10*mean(diff(x$times))))
#   res <- optim(guesses, func_for_optim,method = "L-BFGS-B",
#                lower = c(0,0,0),upper = c(1,1E4, 1E4),
#                control = list(factr = 1E-10, maxit = 1E3, trace=1, REPORT=1))
#   res
  nlminb(guesses, func_for_optim, lower = c(0,0,0), upper = c(1,2.5,1))
  # mle = list(q = res$par[1], beta = res$par[2], lambda = res$par[3], log_like = -1*res$value)
  # return(mle)
  
#   constrOptim(guesses, func_for_optim, grad = NULL, method = 'L-BFGS-B',
#               ui = matrix(c(1,0,0,-1,0,0,0,1,0,0,0,1),4,3,byrow = T),
#               ci = c(0,-1,0,0), control = list(factr = 1E-10, maxit = 1E3))
  
  # func_for_optim <- function(x = NULL, index = NULL, fmsfundata = NULL){
    # return(list(f = ))
  # }
}

