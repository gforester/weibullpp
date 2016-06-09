library(weibullpp)
context("dist_class object test")


test_that("object creation working", {
  set.seed(10)
  test_class = dist_class(c('mean','sd'),
                          function(x,theta, ...){dnorm(x,theta[1],theta[2],...)}, 
                          function(x,theta, ...){pnorm(x,theta[1],theta[2],...)}, 
                          function(x,theta, ...){qnorm(x,theta[1],theta[2],...)},
                          function(x,theta, ...){rnorm(x,theta[1],theta[2],...)},
                          function(x){return(c(mean(x[[1]]),sd(x[[1]])))}
                          )
  expect_equal(test_class$d(1,c(0,1)), dnorm(1))
  expect_equal(test_class$p(1,c(0,1)), pnorm(1))
  expect_equal(test_class$q(0.5,c(0,1)), 0)
})

test_that("likelihood calculation", {
  x <- lifedata(c(22.5,37.5,46,48.5,51.5,53,54.5,57.5,66.5,68,69.5,76.5,77,78.5,80,81.5,82,83,84,
    91.5,93.5,102.5,107,108.5,112.5,113.5,116,117,118.5,119,120,122.5,123,127.5,131,132.5,134,rep(135, 59)),
    c(rep(1,37),rep(0,59))
  )
  test_class = dist_class(c('meanlog','sdlog'),
                          function(x,theta, ...){dlnorm(x,theta[1],theta[2],...)}, 
                          function(x,theta, ...){plnorm(x,theta[1],theta[2],...)}, 
                          function(x,theta, ...){qlnorm(x,theta[1],theta[2],...)},
                          function(x,theta, ...){rlnorm(x,theta[1],theta[2],...)},
                          function(x){return(c(mean(log(x[[1]])), sd(log(x[[1]]))))}
                          )
  expect_equal(round(log_like_generic(x,test_class,c(5.12,0.706)),2),-237.09)
})

test_that("maximum likelihood fitting", {
  x <- lifedata(c(22.5,37.5,46,48.5,51.5,53,54.5,57.5,66.5,68,69.5,76.5,77,78.5,80,81.5,82,83,84,
                  91.5,93.5,102.5,107,108.5,112.5,113.5,116,117,118.5,119,120,122.5,123,127.5,131,132.5,134,rep(135, 59)),
                c(rep(1,37),rep(0,59))
  )
  test_class = dist_class(c('meanlog','sdlog'),
                          function(x,theta, ...){dlnorm(x,theta[1],theta[2],...)}, 
                          function(x,theta, ...){plnorm(x,theta[1],theta[2],...)}, 
                          function(x,theta, ...){qlnorm(x,theta[1],theta[2],...)},
                          function(x,theta, ...){rlnorm(x,theta[1],theta[2],...)},
                          function(x){return(c(mean(log(x[[1]])), sd(log(x[[1]]))))}
  )
  test_res <- mle_generic(x, test_class)
  expect_equal(round(test_res$meanlog,3),5.117)
  expect_equal(round(test_res$sdlog,4),0.7055)
  expect_equal(round(test_res$log_like,2),-237.09)
})