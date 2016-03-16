
#brute-force R method for log-likelihood
#theta is a vector with shape in first entry and shape in second entry
#the data should be a Surv object in surv_object
weibull_2p_like_1 <- function(theta, surv_object){
  full_like <- sum(sapply(X = surv_object[surv_object[,2] == 1,1], FUN = dweibull, shape = theta[1], scale = theta[2], log = T))
  cens_like <- sum(sapply(X = surv_object[surv_object[,2] == 0,1], FUN = pweibull, shape = theta[1], scale = theta[2], lower.tail=F, log.p = T))
  return(full_like + cens_like)
}
