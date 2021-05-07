#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist checklink
#' @importFrom stats dnorm qnorm runif
#' @name VASIQ
#' @aliases VASIQ dVASIQ pVASIQ qVASIQ rVASIQ
#'
#' @title The Vasicek distribution - quantile parameterization
#'
#' @description The function \code{VASIQ()} define the Vasicek distribution for a \code{gamlss.family} object to be used in GAMLSS fitting. \code{VASIQ()} has mean equal to the parameter mu and sigma as shape parameter. The functions \code{dVASIQ}, \code{pVASIQ}, \code{qVASIQ} and \code{rVASIQ} define the density, distribution function, quantile function and random generation for Vasicek distribution.
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
#' @param nu.link the \eqn{\tau}-th fixed quantile level in \code{VASIQ()} family.
#' @param mu vector of location parameter values.
#' @param sigma vector of shape parameter values.
#' @param nu the \eqn{\tau}-th fixed quantile in \[d-p-q-r\]-VASIQ.
#'
#' @return \code{VASIQ()} return a gamlss.family object which can be used to fit a Vasicek distribution in the gamlss() function.
#' 
#' @note Note that for \code{VASIQ()}, mu is the \eqn{\tau}-th quantile and sigma a shape parameter. The \code{\link[gamlss]{gamlss}} function is used for parameters estimation.
#' 
#' @seealso \code{\link[vasicekreg]{VASIM}}.
#'  
#' @details
#' Probability density function 
#' \deqn{f\left( {y\mid \mu ,\sigma ,\tau } \right) = \sqrt {{\textstyle{{1 - \sigma } \over \sigma }}} \exp \left\{ {\frac{1}{2}\left[ {\Phi ^{ - 1} \left( y \right)^2  - \left( {\frac{{\sqrt {1 - \sigma } \left[ {\Phi ^{ - 1} \left( y \right) - \Phi ^{ - 1} \left( \mu  \right)} \right] - \sqrt \sigma  \,\Phi ^{ - 1} \left( \tau  \right)}}{{\sqrt \sigma  }}} \right)^2 } \right]} \right\}}
#' 
#' Cumulative distribution function
#' \deqn{F\left({y\mid \mu ,\sigma ,\tau } \right) = \Phi \left[ {\frac{{\sqrt {1 - \sigma } \left[ {\Phi ^{ - 1} \left( y \right) - \Phi ^{ - 1} \left( \mu  \right)} \right] - \sqrt \sigma  \,\Phi ^{ - 1} \left( \tau  \right)}}{{\sqrt \sigma  }}} \right]}
#' 
#' where \eqn{\mu} is the \eqn{\tau}-th quantile and  \eqn{\sigma} is the shape parameter.

#' @examples
#'
#' set.seed(123)
#' x <- rVASIQ(n = 1000, mu = 0.50, sigma = 0.69, nu = 0.50)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], length.out = 1000)
#' 
#' hist(x, prob = TRUE, main = 'Vasicek')
#' lines(S, dVASIQ(x = S, mu = 0.50, sigma = 0.69, nu = 0.50), col = 2)
#' 
#' plot(ecdf(x))
#' lines(S, pVASIQ(q = S, mu = 0.50, sigma = 0.69, nu = 0.50), col = 2)
#' 
#' plot(quantile(x, probs = S), type = "l")
#' lines(qVASIQ(p = S, mu = 0.50, sigma = 0.69, nu = 0.50), col = 2)
#' 
#' library(gamlss)
#' set.seed(123)
#' data <- data.frame(y =  rVASIQ(n = 100, mu = 0.5, sigma = 0.69, nu = 0.5))
#' 
#' fit <- gamlss(y ~ 1, data = data, mu.link = 'logit', sigma.link = 'logit', 
#' family = VASIQ, nu.start = 0.1)
#' 1 /(1 + exp(-fit$mu.coefficients))
#' 1 /(1 + exp(-fit$sigma.coefficients))
#' 
#' set.seed(123)
#' n <- 100
#' x <- rbinom(n, size = 1, prob = 0.5)
#' eta <- 0.5 + 1 * x;
#' mu <- 1 / (1 + exp(-eta));
#' sigma <- 0.5;
#' y <- rVASIQ(n, mu, sigma, nu = 0.5)
#' data <- data.frame(y, x)
#' 
#' VASIQ(mu.link = "logit", sigma.link = "identity")
#' fit <- gamlss(y ~ x, data = data, family = VASIQ, mu.link = 'probit', 
#' nu.start = 0.5, sigma.link = "log");
##################################################
#' @rdname VASIQ
#' @export
#
dVASIQ <- function (x, mu, sigma, nu = 0.5, log = FALSE)
{
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, sigma > 0, sigma < 1, nu > 0, nu < 1);
  cpp_dvasicekquant(x, mu, sigma, nu, log[1L]);
}
##################################################
#' @rdname VASIQ
#' @export
#'
pVASIQ <- function (q, mu, sigma, nu = 0.50, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(q > 0, q < 1, mu > 0, mu < 1, sigma > 0, sigma < 1, nu > 0, nu < 1);
  cpp_pvasicekquant(q, mu, sigma, nu, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname VASIQ
#' @export
#'
qVASIQ <- function(p, mu, sigma, nu = 0.50, lower.tail = TRUE, log.p = FALSE)
{
  stopifnot(p > 0, p < 1, mu > 0, mu < 1, sigma > 0, sigma < 1, nu > 0, nu < 1);
  cpp_qvasicekquant(p, mu, sigma, nu, lower.tail[1L], log.p[1L])
}
##################################################
#' @rdname VASIQ
#' @export
#'
rVASIQ <- function(n, mu, sigma, nu = 0.50)
{
  cpp_qvasicekquant(runif(n), mu, sigma, nu, TRUE, FALSE)
}

##################################################
#' @rdname VASIQ
#' @export
#' 
VASIQ <- function (mu.link = "logit", sigma.link = "logit", nu.link = "identity")
{
  mstats <- checklink("mu.link", "VasicekQ", substitute(mu.link), c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "VasicekQ", substitute(sigma.link), c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  vstats <- checklink("nu.link", "VasicekQ", substitute(nu.link), c("identity"))
  structure(
    list(family     = c("VASIQ", "VasicekQ"),
         parameters = list(mu = TRUE, sigma = TRUE, nu = FALSE),
         nopar      = 3,
         type       = "Continuous",
         mu.link    = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         nu.link    = as.character(substitute(nu.link)),
         mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
              mu.dr = mstats$mu.eta,
           sigma.dr = dstats$mu.eta,
              nu.dr = vstats$mu.eta,
               dldm = function(y, mu, sigma, nu){
                   t2 <- sqrt(0.1e1 - sigma);
                   t4 <- qnorm(mu);
                   t6 <- qnorm(nu);
                   t7 <- sqrt(sigma);
                   t12 <- 0.1e1 / dnorm(t4);
                   qnormx <- qnorm(y);
                   return(0.10e1 * (qnormx * t2 - t4 * t2 + t6 * t7) / sigma * t12 * t2);
               },
             d2ldm2 = function(y, mu, sigma, nu){
                   t10 <- qnorm(mu);
                   t1 <- 0.1e1 / dnorm(t10)
                   t2 <- t1 * t1;
                   t3 <- 0.1e1 - sigma;
                   t5 <- 0.1e1 / sigma;
                   t8 <- sqrt(t3);
                   t12 <- qnorm(nu);
                   t13 <- sqrt(sigma);
                   t17 <- -t10;
                   qnormx <- qnorm(y);
                   return(-0.10e1 * t2 * t3 * t5 + 0.10e1 * (qnormx * t8 - t10 * t8 + t12 * t13) * t5 * t17 * t8);
               },
               dldd = function(y, mu, sigma, nu){
                   qnormx <- qnorm(y);
                   t1 <- 0.1e1 - sigma;
                   t4 <- 0.1e1 / sigma;
                   t6 <- sqrt(t1);
                   t8 <- qnorm(mu);
                   t10 <- qnorm(nu);
                   t11 <- sqrt(sigma);
                   t13 <- qnormx * t6 - t8 * t6 + t10 * t11;
                   t15 <- 0.1e1 / t6;
                   t23 <- t13 * t13;
                   t24 <- sigma * sigma;
                   return(-0.1e1 / t1 / 0.2e1 - t4 / 0.2e1 - 0.5000000000e0 * t13 * t4 *
                            (-qnormx * t15 + t8 * t15 + t10 / t11) + 0.5e0 * t23 / t24);
               },
               d2ldmdd = function(y, mu, sigma, nu){
                   t2 <- sqrt(0.1e1 - sigma);
                   t3 <- 0.1e1 / t2;
                   t5 <- qnorm(mu);
                   t7 <- qnorm(nu);
                   t8 <- sqrt(sigma);
                   t12 <- 0.1e1 / sigma;
                   t14 <- 0.1e1 / dnorm(t5);
                   t15 <- t14 * t2;
                   qnormx <- qnorm(y);
                   t21 <- qnormx * t2 - t5 * t2 + t7 * t8;
                   t22 <- sigma * sigma;
                   return(0.5000000000e0 * (-qnormx * t3 + t5 * t3 + t7 / t8) * t12 * t15 - 0.10e1 *
                            t21 / t22 * t15 - 0.5000000000e0 * t21 * t12 * t14 * t3);
               },
               d2ldd2 = function(y, mu, sigma, nu){
                   qnormx <- qnorm(y);
                   t1 <- 0.1e1 - sigma;
                   t2 <- t1 * t1;
                   t5 <- sigma * sigma;
                   t6 <- 0.1e1 / t5;
                   t8 <- sqrt(t1);
                   t9 <- 0.1e1 / t8;
                   t11 <- qnorm(mu);
                   t13 <- qnorm(nu);
                   t14 <- sqrt(sigma);
                   t17 <- -qnormx * t9 + t11 * t9 + t13 / t14;
                   t19 <- 0.1e1 / sigma;
                   t25 <- qnormx * t8 - t11 * t8 + t13 * t14;
                   t31 <- 0.1e1 / t8 / t1;
                   t40 <- t25 * t25;
                   return(-0.1e1 / t2 / 0.2e1 + t6 / 0.2e1 - 0.2500000000e0 * t17 * t17 * t19 + 0.1000000000e1 *
                            t25 * t6 * t17 - 0.2500000000e0 * t25 * t19 * (-qnormx * t31 + t11 * t31 - t13 / t14 / sigma) -
                            0.10e1 * t40 / t5 / sigma);
               },
              dldv = function(y)  rep(0,length(y)),
            d2ldv2 = function(y)  rep(0,length(y)),
           d2ldmdv = function(y)  rep(0,length(y)),
           d2ldddv = function(y)  rep(0,length(y)),
        G.dev.incr = function(y, mu, sigma, nu, w, ...) -2 * dVASIQ(y, mu, sigma, nu, log = TRUE),
             rqres = expression(rqres(pfun = "pVASIQ", type = "Continuous", y = y, mu = mu, sigma = sigma, nu = nu)),
        mu.initial = expression({mu <- (y + mean(y))/2}),
     sigma.initial = expression({sigma <- rep(0.5, length(y))}),
        nu.initial = expression({nu <- rep(0.5,length(y))}),
          mu.valid = function(mu) all(mu > 0 & mu < 1),
       sigma.valid = function(sigma) all(sigma > 0 & sigma < 1),
          nu.valid = function(nu) all(nu > 0 & nu < 1),
           y.valid = function(y) all(y > 0 & y < 1),
              mean = function(mu) mu,
          variance = function(mu,sigma) sigma^2 * mu * (1 - mu)),
    class = c("gamlss.family","family"))
}
