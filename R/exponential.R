
#right censored only
fit_exp <- function(x){
  mle_scale <- sum(x[[1]]/sum(x[[2]]))
  full_like <- sum(dexp(x[[1]][x[[2]]==1],rate = 1/mle_scale, log = T))
  cens_like <- sum(pexp(x[[1]][x[[2]]==0],rate = 1/mle_scale, lower.tail=F, log.p = T))
  return(list(scale = mle_scale,
              log_like = full_like + cens_like
              ))
}
# 

