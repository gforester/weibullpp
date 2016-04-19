


cum_times <- c(50, 94, 196, 268, 290, 329, 332, 347, 544, 732, 811, 899, 945,
               950, 955, 991, 1013, 1152, 1362, 1459, 1489, 1512, 1525, 1539)

#Cox, F. R., and Lewis, P.A. W. (1966), The Statistical Analysis of Series of Events, London: Methuen.



reliasoft_1 <- sort(c(
  0.7,	121.9,	285,
  3.7,	125.5,	304,
  13.2,	133.4,	315.4,
  15,	151,	317.1,
  17.6,	163,	320.6,
  25.3,	164.7,	324.5,
  47.5,	174.5,	324.9,
  54,	177.4,	342,
  54.5,	191.6,	350.2,
  56.4,	192.7,	355.2,
  63.6,	213,	364.6,
  72.2,	244.8,	364.9,
  99.2,	249,	366.3,
  99.6,	250.8,	373,
  100.3,	260.1,	379.4,
  102.5,	263.5,	389,
  112,	273.1,	394.9,
  112.2,	274.7,	395.2,
  120.9,	282.8	
))

#http://www.reliasoft.com/Weibull/examples/rc9/index.htm
# reliasoft solution - q = 0.931569 beta = 0.894429 lambda =  0.247257, log like = -165.247296
# likelihood value of package matches but finds another optimum
# q = 0.005887538 beta =  0.981035391 lambda = 0.148333906, log like = -165.4113