

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




