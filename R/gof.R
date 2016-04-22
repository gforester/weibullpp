

# modified_ks_test <- function(fail_times, fitted_cdf){
#   #create step function
#   ts <- sort(fail_times)
#   steps <- c(0,cumsum(as.numeric(table(ts)/length(ts))))
#   ts <- unique(ts)
#   Sn <- stepfun(ts, steps)
#   to_optim <- function(x){abs(Sn(x)-fitted_cdf(x))}
#   res <- optimize(f = to_optim, interval = c(0.99*min(ts),1.01*max(ts)), maximum = T)
#   D <- res$objective
#   .Call(stats:::C_pKolmogorov2x,STAT=D,n=length(fail_times))
# }


chi_square_test <- function(fail_times, fitted_cdf, n_fitted){
  N <- length(fail_times)
  breaks <- pretty(fail_times, 1+3.322*log10(N))
  chisq_stat <- 0
  for(i in 1:(length(breaks)-1)){
    actual <- sum(fail_times >= breaks[i] & fail_times < breaks[i+1])
    expected <- N*( fitted_cdf(breaks[i+1]) - fitted_cdf(breaks[i]) )
    chisq_stat <- chisq_stat + (((expected - actual)^2)/expected)
  }
  dof <- length(breaks)-1-n_fitted-1
  p_value <- pchisq(chisq_stat,dof,lower.tail=F)
  chisq_crit <- qchisq(0.05, dof, lower.tail=F)
  return(list(chisq_stat = chisq_stat, chisq_crit = chisq_crit, dof = dof,p_value = p_value))
}

