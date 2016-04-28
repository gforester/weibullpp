#' @title Fit Distribution
#' 
#' @description 
#' fits distribution to data
#' 
#' @details
#' refer to \code{\link{fit_data.lifedata}}
#' 
#' @param x only defined for a lifedata object
#' 
#' @seealso \code{\link{fit_data.lifedata}}
fit_data <- function(x,...){
  UseMethod('fit_data')
}

pieplot <- function(x, ...) UseMethod('pieplot')
pieplot.default <- function(x, ...) pie(x, ...)
