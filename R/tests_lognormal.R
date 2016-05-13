
# from http://reliawiki.org/index.php/The_Lognormal_Distribution#RRY_Example
# lognorm_test1 <- lifedata(c(5,10,15,20,25,30,35,40,50,60,70,80,90,100),rep(1,14))
# to_test <- fit_data(lognorm_test1, 'lognormal')
# cat('Expect logmean',3.516,'got',round(to_test$fit$logmean,3),'\n')
# cat('Expect sdlog',0.849,'got',round(to_test$fit$sdlog,3),'\n')
# cat('Expect se on logmean',round(sqrt(0.0515),3),'got',round(to_test$std_error$logmean_delta,3),'\n')
# cat('Expect se on sdlog',round(sqrt(0.0258),3),'got',round(to_test$std_error$sdlog_delta,3),'\n')

# # from http://reliawiki.org/index.php/The_Lognormal_Distribution#Suspension_Data_Example
# #this was base 10 log so don't use 
# lognorm_test2 <- function(){
  # x <- lifedata(c(22.5,37.5,46,48.5,51.5,53,54.5,57.5,66.5,68,69.5,76.5,77,78.5,80,81.5,82,83,84,
  #   91.5,93.5,102.5,107,108.5,112.5,113.5,116,117,118.5,119,120,122.5,123,127.5,131,132.5,134,rep(135, 59)),
  #   c(rep(1,37),rep(0,59))
  # )
  # to_test <- fit_data(x, 'lognormal')
#   cat('Expect logmean',2.2223,'got',round(to_test$fit$logmean,4),'\n')
#   cat('Expect sdlog',0.3064,'got',round(to_test$fit$sdlog,4),'\n')
#   cat('Expect se on logmean',round(sqrt(0.0019),4),'got',round(to_test$std_error$logmean_delta,4),'\n')
#   cat('Expect se on sdlog',round(sqrt(0.0015),4),'got',round(to_test$std_error$sdlog_delta,4),'\n')
# }
# 
