#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist checklink
#' @importFrom stats dnorm qnorm runif
#' @name VASIM
#' @aliases VASIM dVASIM pVASIM qVASIM rVASIM
#'
#' @title The Vasicek distribution - mean parameterization
#'
#' @description The function \code{VASIM()} define the Vasicek distribution for a \code{gamlss.family} object to be used in GAMLSS fitting. \code{VASIM()} has mean equal to the parameter mu and sigma as shape parameter. The functions \code{dVASIM}, \code{pVASIM}, \code{qVASIM} and \code{rVASIM} define the density, distribution function, quantile function and random generation for Vasicek distribution.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#' @author Bruna Alves \email{pg402900@uem.br}
#' 
#' @references
#'  
#' Hastie, T. J. and Tibshirani, R. J. (1990). \emph{Generalized Additive Models}. Chapman and Hall, London.
#' 
#' Mazucheli, J., Alves, B. and Korkmaz, M. Ç. (2021). The Vasicek quantile regression model. \emph{(under review)}.
#' 
#' Rigby, R. A. and  Stasinopoulos, D. M. (2005). Generalized additive models for location, scale and shape (with discussion). \emph{Applied. Statistics}, \bold{54}(3), 507--554.
#' 
#' Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z. and De Bastiani, F. (2019). \emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}. Chapman and Hall/CRC.
#' 
#' Stasinopoulos, D. M. and Rigby, R. A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. \emph{Journal of Statistical Software}, \bold{23}(7), 1--45.
#' 
#' Stasinopoulos, D. M., Rigby, R. A., Heller, G., Voudouris, V. and De Bastiani F. (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC.  
#' 
#' Vasicek, O. A. (1987). Probability of loss on loan portfolio. \emph{KMV Corporation}.
#'
#' Vasicek, O. A. (2002). The distribution of loan portfolio value. \emph{Risk}, \bold{15}(12), 1--10.
#' 
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.
#' @param mu.link the mu link function with default logit.
#' @param sigma.link the sigma link function with default logit.
#' @param mu vector of location parameter values, the mean.
#' @param sigma vector of shape parameter values.
#'
#' @return \code{VASIM()} return a gamlss.family object which can be used to fit a Vasicek distribution in the gamlss() function.
#' 
#' @note Note that for \code{VASIQ()}, mu is the \eqn{\tau}-th quantile and sigma a shape parameter. The \code{\link[gamlss]{gamlss}} function is used for parameters estimation.
#' 
#' @seealso \code{\link[vasicekreg]{VASIQ}}.
#'
#' @details
#' Probability density function 
#' \deqn{f(y\mid \mu ,\sigma )=\sqrt{\frac{1-\sigma }{\sigma }}\exp \left\{ \frac{1}{2}\left[ \Phi ^{-1}\left( y\right) ^{2}-\left( \frac{\Phi ^{-1}\left(  y\right)    \sqrt{1-\sigma }-\Phi ^{-1}\left( \mu \right) }{\sqrt{\sigma }}\right) ^{2}\right] \right\}}
#' 
#' Cumulative distribution function
#' \deqn{F(y\mid \mu ,\sigma )=\Phi \left( \frac{\Phi ^{-1}\left( y\right) \sqrt{1-\sigma }-\Phi ^{-1}\left( \mu \right) }{\sqrt{\sigma }}\right)}
#' 
#' Quantile function
#' \deqn{Q(\tau \mid \mu ,\sigma )=F^{-1}(\tau \mid \mu ,\sigma )=\Phi \left(\frac{\Phi ^{-1}\left(\mu\right) +\Phi ^{-1}\left( \tau \right) \sqrt{\sigma }}{\sqrt{1-\sigma }}\right) }
#' where \eqn{\mu = E(Y)}, \eqn{\sigma} is the shape parameter and  \eqn{0<\tau<1} is the \eqn{\tau}-th quantile.

#' @examples
#'
#' set.seed(123)
#' x <- rVASIM(n = 1000, mu = 0.50, sigma = 0.69)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], length.out = 1000)
#' 
#' hist(x, prob = TRUE, main = 'Vasicek')
#' lines(S, dVASIM(x = S, mu = 0.50, sigma = 0.69), col = 2)
#' 
#' plot(ecdf(x))
#' lines(S, pVASIM(q = S, mu = 0.50, sigma = 0.69), col = 2)
#' 
#' plot(quantile(x, probs = S), type = "l")
#' lines(qVASIM(p = S, mu = 0.50, sigma = 0.69), col = 2)
#' 
#' library(gamlss)
#' set.seed(123)
#' data <- data.frame(y =  rVASIM(n = 100, mu = 0.5, sigma = 0.69))
#' 
#' fit <- gamlss(y ~ 1, data = data, mu.link = 'logit', sigma.link = 'logit', family = VASIM)
#' 1 /(1 + exp(-fit$mu.coefficients))
#' 1 /(1 + exp(-fit$sigma.coefficients))
#' 
#' set.seed(123)
#' n <- 100
#' x <- rbinom(n, size = 1, prob = 0.5)
#' eta <- 0.5 + 1 * x;
#' mu <- 1 / (1 + exp(-eta));
#' sigma <- 0.1;
#' y <- rVASIM(n, mu, sigma)
#' data <- data.frame(y, x)
#' 
#' fit <- gamlss(y ~ x, data = data, family = VASIM, mu.link = 'logit', sigma.link = 'logit');
#' 
##################################################
#' @rdname VASIM
#' @export
#
dVASIM <- function (x, mu, sigma, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, sigma > 0, sigma < 1);
  cpp_dvasicekmean (x, mu, sigma, log[1L]);
}
##################################################
#' @rdname VASIM
#' @export
#'
pVASIM <- function (q, mu, sigma, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, sigma > 0, sigma < 1)
  cpp_pvasicekmean (q, mu, sigma, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname VASIM
#' @export
#'
qVASIM <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, sigma > 0, sigma < 1)
  cpp_qvasicekmean (p, mu, sigma, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname VASIM
#' @export
#'
rVASIM <- function(n, mu, sigma)
{
  cpp_qvasicekmean (runif(n), mu, sigma, TRUE, FALSE)
}

##################################################
#' @rdname VASIM
#' @export
#' 
VASIM <- function (mu.link = "logit", sigma.link = "logit")
{
  mstats <- checklink("mu.link", "VASIM", substitute(mu.link), c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "VASIM", substitute(sigma.link), c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  structure(list(family = c("VASIM", "Vasicekm"),
                 parameters = list(mu = TRUE, sigma = TRUE), nopar = 2, type = "Continuous", mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)), mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv, mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta,
                 dldm = function(y, mu, sigma){
                   t2 <- sqrt(0.1e1 - sigma);
                   t4 <- qnorm(mu);
                   t8 <- 0.1e1 / dnorm(t4);
                   qnormx <- qnorm(y);
                   return(0.10e1 * (qnormx * t2 - t4) / sigma * t8);
                 },
                 d2ldm2 = function(y, mu, sigma){
                   t9 <- qnorm(mu);
                   t1 <- 0.1e1 / dnorm(t9);
                   t2 <- t1 * t1;
                   t3 <- 0.1e1 / sigma;
                   t7 <- sqrt(0.1e1 - sigma);
                   t12 <-  t9 * (t1 ^ 0.2e1);
                   qnormx <- qnorm(y);
                   return(-0.10e1 * t2 * t3 + 0.10e1 * (qnormx * t7 - t9) * t3 * t12);
                 },
                 dldd = function(y, mu, sigma){
                   qnormx <- qnorm(y);
                   t1 <- 0.1e1 - sigma;
                   t4 <- 0.1e1 / sigma;
                   t6 <- sqrt(t1);
                   t8 <- qnorm(mu);
                   t9 <- qnormx * t6 - t8;
                   t15 <- t9 * t9;
                   t16 <- sigma * sigma;
                   return(-0.5e0 / t1 - 0.5e0 * t4 + 0.5e0 * t9 * t4 * qnormx / t6 + 0.5e0 * t15 / t16);
                 },
                 d2ldd2 = function(y, mu, sigma){
                   qnormx <- qnorm(y);
                   t1 <- 0.1e1 - sigma;
                   t2 <- t1 * t1;
                   t5 <- sigma * sigma;
                   t6 <- 0.1e1 / t5;
                   t8 <- qnormx * qnormx;
                   t11 <- 0.1e1 / sigma;
                   t14 <- sqrt(t1);
                   t16 <- qnorm(mu);
                   t17 <- qnormx * t14 - t16;
                   t29 <- t17 * t17;
                   return(-0.5e0 / t2 + 0.5e0 * t6 - 0.2500000000e0 * t8 / t1 * t11 - 0.10e1 * t17 * t6 * qnormx / t14 + 
                            0.2500000000e0 * t17 * t11 * qnormx / t14 / t1 - 0.10e1 * t29 / t5 / sigma);
                 },
                 d2ldmdd = function(y, mu, sigma){
                   t1 <- qnorm(mu);
                   t2 <- sqrt(0.1e1 - sigma);
                   t6 <- 0.1e1 / dnorm(t1)
                   t13 <- sigma * sigma;
                   qnormx <- qnorm(y)
                   return(-0.5000000000e0 * qnormx / t2 / sigma * t6 - 0.10e1 * (qnormx * t2 - t1) / t13 * t6);
                 },
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * dVASIM(y, mu, sigma, log = TRUE),
                 rqres = expression(rqres(pfun = "pVASIM", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                 mu.initial = expression({mu <- (y + mean(y))/2}),
                 sigma.initial = expression({sigma <- rep(0.5, length(y))}),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
                 y.valid = function(y) all(y > 0 & y < 1),
                 mean = function(mu, sigma) mu, variance = function(mu,sigma) sigma^2 * mu * (1 - mu)),class = c("gamlss.family","family"))
}


