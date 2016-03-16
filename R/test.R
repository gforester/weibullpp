

test_weibullpp = 5

#Reference Case
#The data set in Table C.5 on page 633 of the book Statistical Methods for Reliability Data by Dr. Meeker and Dr. Escobar, John Wiley & Sons, 1998 is used.
# Result
# The results are given on page 193 in Example 8.16.
# The MLE solution for  is 2.035.
# The B10 life  is 3903. Its 95% 2-sided likelihood ratio bounds are [2093, 22, 144].
test1 <- data.frame(
  n_in_state = c(288, 148, 1, 124, 1, 111, 1, 106, 99, 110, 114, 119, 127, 1, 1, 123, 93, 47, 41, 27, 1, 11, 6, 1, 2),
  full_or_cens = c(0,0,1,0,1,0,1,0,0,0,0,0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
  time = c(50, 150, 230, 250, 334, 350, 423, 450, 550, 650, 750, 850, 950, 990, 1009, 1050, 1150, 1250, 1350, 1450, 1510, 1550, 1650, 1850, 2050)
)
test1_surv <- Surv(rep(test1$time, test1$n_in_state), rep(test1$full_or_cens, test1$n_in_state))
#weibull++ answers
#shape = 2.035247
#scale = 11793.21105
#like = -76.436896
#survival survreg answer #res = survreg(test1_surv ~ 1)
#shape = 2.035319 #1/exp(res$icoef[2])
#scale = 11792.17817
#like = -76.436896
#weibulpp answers
#shape = 2.035229
#scale = 11793.46551
#like = -76.436896


