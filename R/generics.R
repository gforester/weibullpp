#' @title Fit Distribution
#' 
#' @description 
#' fits distribution to data
#' 
#' @details
#' refer to fit_data.lifedata
#' 
#' @param x only defined for a lifedata object
fit_data <- function(x,...){
  UseMethod('fit_data')
}

pieplot <- function(x, ...) UseMethod('pieplot')
