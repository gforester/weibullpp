
# test 1 - complete data
run_test_gof_1 <- function(){
  set_gof_1 <- lifedata(c(43,68,74,77,80,91,99,103,103,166), rep(1,10))
  fit_gof_1_exp <- fit_data(set_gof_1, dist = 'exponential')
  fit_gof_1_wei <- fit_data(set_gof_1)
  cat('Exponential Fit \n')
  cat('expect scale', round(1/0.0111,4), 'got', round(fit_gof_1_exp$fit$scale,4), "\n")
  cat('expect log-like', -55.04, 'got', round(fit_gof_1_exp$fit$log_like,2), "\n")
  cat('Weibull Fit \n')
  cat('expect shape', 3.03, 'got', round(fit_gof_1_wei$fit$shape,2), "\n")
  cat('expect scale', 100.99, 'got', round(fit_gof_1_wei$fit$scale,2), "\n")
  cat('expect log-like', -48.42, 'got', round(fit_gof_1_wei$fit$log_like,2), "\n")
}


