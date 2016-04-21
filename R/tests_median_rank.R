
# DATA AND RESULTS SOURCE - Reliasoft's Life Data Analysis Reference, 
#     Reliasoft, Reliability Engineering eTextbook Library, 
#     url:http://reliawiki.org/index.php/Parameter_Estimation

# TEST 1 - All complete data
# set1 <- lifedata(c(16,34,53,75,93,120), rep(1,6))
# res1 <- c(0.1091, 0.2644, 0.4214, 0.5786, 0.7356, 0.8910)
# tes1a <- median_ranks(set1, method = 'inverse_binom')
# tes1b <- median_ranks(set1)
# tes1c <- median_ranks(set1, method = 'benard')
# cat(paste0(paste('expect', res1, 'got', round(tes1a,4)), collapse="\n"))
# cat(paste0(paste('expect', res1, 'got', round(tes1b,4)), collapse="\n"))
# cat(paste0(paste('expect', res1, 'got', round(tes1c,4)), collapse="\n"))

# TEST 2 - Complete and Right Censored
# set2 <- lifedata(c(5100,9500, 15000, 22000,40000), c(1,0,1,0,1))
# res2 <- c(0.13, 0.36, 0.71)
# tes2 <- median_ranks(set2)
# cat(paste0(paste('expect', res2, 'got', round(tes2,2)), collapse="\n"))
