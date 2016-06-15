# 
# default_dist_list <- list(
#   'weibull' = structure(list(thetas = c('shape','scale'), 
#                              d = function(x,theta, log = F) {dweibull(x,theta[1],theta[2], log)},
#                              p = function(x,theta, lower.tail= T, log.p = F) {pweibull(x, theta[1], theta[2], lower.tail, log.p)},
#                              q = function(x,theta, lower.tail= T, log.p = F) {qweibull(x, theta[1], theta[2], lower.tail, log.p)},
#                              r=r,guess = g, rrs = rrs, lins = lins), 
#                         class="dist_class"),
#   'exponential' = structure(list(thetas = c('scale'), 
#                                  d=d,p=p,q=q,r=r,guess = g, rrs = rrs, lins = lins), 
#                             class="dist_class"),
#   'lognormal' = structure(list(thetas = c('logmean','sdlog'), 
#                                d=d,p=p,q=q,r=r,guess = g, rrs = rrs, lins = lins), 
#                           class="dist_class")
#   
# )

#check that dist is in proper format!!
