#' @docType package
#' @name vasicekreg-package
#' @aliases vasicekreg-package
#'
#' @title Overview of the vasicekreg package
#'
#' @description The \pkg{vasicekreg} package implements the probability density function, quantile function, cumulative distribution function and random number generation function for Vasicek distribution parameterized, either, as a function of its mean or its \eqn{\tau}-th quantile, \eqn{0 <\tau<1}. In addition, two gamlss frameworks for regression analysis are available. Some function are written in \proglang{C++} using \pkg{Rcpp}. 
#' 
#' @details 
#' 
#' \code{\link[vasicekreg]{bodyfat}}: Body fat data set.
#' 
#' \code{\link[vasicekreg]{VASIM}}: For mean modeling (con/in)ditional on covariate(s).
#' 
#' \code{\link[vasicekreg]{VASIQ}}: For quantile modeling (con/in)ditional on covariate(s).
#' 
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Bruna Alves \email{pg402900@uem.br}
#' 
#'
#' @useDynLib vasicekreg
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
    library.dynam.unload("vasicekreg", libpath)
}
