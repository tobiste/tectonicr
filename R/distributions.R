# Distribution ####

## Circular Uniform ------------------------------------------------------------

#' The Circular Uniform Distribution
#'
#' Density, probability distribution function, quantiles, and random generation
#' for the circular uniform distribution.
#'
#' @inheritParams stats::runif
#' @param p numeric. Vector of probabilities with values in \eqn{[0,1]}{[0,1]}.
#' @param theta numeric. Angular value in degrees
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @returns `dcunif` gives the density,
#' `pcunif` gives the probability of circular uniform distribution function,
#' `rcunif` generates random deviates (in degrees), and
#' `qcunif` provides quantiles (in degrees).
#'
#' @seealso [wnorm], [wcauchy], and [vonmises]
#'
#' @name cunif
#'
#' @importFrom stats runif
#'
#' @examples
#' set.seed(1)
#' x <- rcunif(5)
#'
#' dcunif(x)
#' dcunif(x, axial = TRUE)
#'
#' pcunif(x)
#' qcunif(c(.25, .5, .75))
NULL

#' @rdname cunif
#' @export
rcunif <- function(n, axial = FALSE){
  f <- if(isTRUE(axial)) 1 else 2
  stats::runif(n = n, min = 0, max = f * 180)
}

#' @rdname cunif
#' @export
dcunif <- function(theta, axial = FALSE, log = FALSE){
  f <- if(isTRUE(axial)) 1 else 2
  d <- rep(1/(f * 180), length(theta))
  if(isTRUE(log)) log(d) else d
}

#' @rdname cunif
#' @export
pcunif <- function(theta, axial = FALSE, lower.tail = TRUE, log.p = FALSE) {
  f <- if (isTRUE(axial)) 1 else 2
  period <- f * 180
  p <- (theta %% period) / period
  if (!isTRUE(lower.tail)) p <- 1 - p
  if (isTRUE(log.p)) p <- log(p)
  p
}

#' @rdname cunif
#' @export
qcunif <- function(p, axial = FALSE, lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(log.p)) p <- exp(p)
  if (!isTRUE(lower.tail)) p <- 1 - p
  bad <- p < 0 | p > 1
  if (any(bad, na.rm = TRUE)) warning("NaNs produced")
  p[bad] <- NaN
  f <- if (isTRUE(axial)) 1 else 2
  p * f * 180
}

## von Mises -------------------------------------------------------------------

#' The von Mises Distribution
#'
#' Density, probability distribution function, quantiles, and random generation
#' for the circular normal distribution with mean \eqn{\mu}{mu} and kappa \eqn{\kappa}{k}.
#'
#' @inheritParams cunif
#' @param mean numeric. The mean vector in degrees.
#' @param p numeric. Vector of probabilities with values in \eqn{[0,1]}{[0,1]}.
#' @param kappa numeric. Concentration parameter in the range \eqn{[0, Inf]}
#' @param from if `NULL` is set to \eqn{\mu-\pi}{mu-pi}. This is the value from
#' which the `pvm` and `qvm` are evaluated. in degrees.
#' @param log logical. If `TRUE`, probabilities p are given as \eqn{\log(p)}{log(p)}.
#' @param tol numeric. The precision in evaluating the distribution function or the quantile.
#' @param ... parameters passed to [stats::integrate()].
#'
#' @returns `dvm` gives the density,
#' `pvm` gives the probability of the von Mises distribution function,
#' `rvm` generates random deviates (in degrees), and
#' `qvm` provides quantiles (in degrees).
#'
#' @seealso [cunif], [wnorm], and [wcauchy]
#' @name vonmises
#'
#' @importFrom circular circular rvonmises pvonmises qvonmises daxialvonmises
#'
#' @examples
#' set.seed(1)
#' x <- rvm(5, mean = 90, kappa = 2)
#'
#' dvm(x, mean = 90, kappa = 2)
#' dvm(x, mean = 90, kappa = 2, axial = TRUE)
#'
#' pvm(x, mean = 90, kappa = 2)
#' qvm(c(.25, .5, .75), mean = 90, kappa = 2)
NULL

#' @rdname vonmises
#' @export
rvm <- function(n, mean, kappa) {
  mu <- circular::circular(mean, units = "degrees", modulo = "2pi")
  circular::rvonmises(n, mu, kappa) |> as.numeric()
}

#' @rdname vonmises
#' @export
dvm <- function(theta, mean, kappa, axial = FALSE, log = FALSE) {
  f <- if(axial) 2 else 1
  two_pi <- 2 * pi
  x <- deg2rad(f * theta) %% two_pi
  mu <- deg2rad(f * mean) %% two_pi

  delta <- x - mu
  delta_mod <- delta %% two_pi

  n <- length(x)

  if (kappa == 0) {
    vm <- rep(1 / two_pi, n)
  } else if (kappa < 1e+05) {
    bessel_val <- besselI(kappa, nu = 0, expon.scaled = TRUE)
    vm <- (1 / (two_pi * bessel_val)) * (exp(cos(delta) - 1))^kappa
  } else {
    vm <- ifelse(delta_mod == 0, Inf, 0)
  }

  if (isTRUE(axial)) vm <- f * vm
  if (isTRUE(log)) vm <- log(vm)

  return(vm)
}

#' @rdname vonmises
#' @export
pvm <- function(theta, mean, kappa, from = NULL, tol = 1e-20) {
  theta <- circular::circular(theta, units = "degrees", modulo = "2pi")
  mu <- circular::circular(mean, units = "degrees", modulo = "2pi")

  if (!is.null(from)) {
    from <- circular::circular(from, units = "degrees", modulo = "2pi")
  }

  circular::pvonmises(theta, mu, kappa, from = from, tol = tol)
}

#' @rdname vonmises
#' @export
qvm <- function(p, mean, kappa, from = NULL, tol = .Machine$double.eps^(0.6), ...) {
  mu <- circular::circular(mean, units = "degrees", modulo = "2pi")

  if (!is.null(from)) {
    from <- circular::circular(from, units = "degrees", modulo = "2pi")
  }

  circular::qvonmises(p, mu, kappa, from, tol = tol, ...) |>
    as.numeric()
}


#' @keywords internal
A1inv <- function(x) {
  stopifnot(is.numeric(x), length(x) == 1, !is.na(x))
  if (0 <= x && x < 0.53) {
    2 * x + x^3 + (5 * x^5) / 6
  } else if (x < 0.85) {
    -0.4 + 1.39 * x + 0.43 / (1 - x)
  } else {
    1 / (x^3 - 4 * x^2 + 3 * x)
  }
}

#' Concentration parameter of von Mises distribution
#'
#' Estimates the concentration parameter of a von Mises distribution, given a
#' set of angular measurements.
#'
#' @param x numeric. Angles in degrees
#' @param w numeric. Weightings
#' @param bias logical parameter determining whether a bias correction is used
#' in the computation of the MLE. Default for bias is `FALSE` for no bias
#' correction.
#'
#' @details
#' `est.kappa.MLE()` is the maximum likelihood estimate for MLE for \eqn{\kappa}.
#'
#' `est.kappa()` uses an approximation based on the empirical equation:
#' \deqn{\kappa =
#'     \frac{\bar{R}(p-\bar{R}^2)}{1-\bar{R}^2}
#' }
#' where \eqn{\bar{R}} is the mean resultant length and \eqn{p} is the
#' dimensionality of the data (2 for circular data).
#'
#' @returns numeric. Concentration of a von Mises distribution
#' @name estimate-kappa
#'
#' @seealso [vonmises()]
#'
#' @examples
#' set.seed(123)
#' x <- rvm(100, 90, 10)
#' w <- weighting(runif(100, 0, 10))
#'
#' est.kappa(x, w)
#'
#' est.kappa.MLE(x, w)
NULL

#' @rdname estimate-kappa
#' @export
est.kappa.MLE <- function(x, w = NULL, bias = FALSE) {
  # Default weights
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)
  }

  #f <- if (axial) 2 else 1
  # x <- (x * f) %% 360

  # Remove NA pairs
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]

  mean.dir <- circular_mean(x, w = w, axial = FALSE, na.rm = FALSE)
  mean_cos <- mean(cosd(x - mean.dir))
  kappa <- abs(A1inv(mean_cos))

  if (bias) {
    n <- sum(w)
    if (kappa < 2) {
      kappa <- max(kappa - 2 / (n * kappa), 0)
    }
    if (kappa >= 2) {
      kappa <- ((n - 1)^3 * kappa) / (n^3 + n)
    }
  }
  kappa
}

#' @rdname estimate-kappa
#' @param p integer. Number of parameters in the data space: 2 for circle (the default), 3 for a sphere.
#' @export
est.kappa <- function(x, w = NULL, p = 2) {
  Rbar <- mean_resultant_length(x, w, na.rm = TRUE)
  (Rbar * (p - Rbar^2)) / (1 - Rbar^2)
}

## Wrapped Normal --------------------------------------------------------------

#' The Wrapped Normal Distribution
#'
#' Density, probability distribution function, quantiles, and random generation
#' for the circular normal distribution with mean \eqn{\mu}{mu} and standard deviation \eqn{\sigma}{sd}.
#'
#' @inheritParams vonmises
#' @param sd 	numeric. standard deviation of the (unwrapped) normal distribution in degrees.
#' @param ... optional parameters passed to underlying circular functions
#' [circular::pwrappednormal()] and [circular::qwrappednormal()]
#'
#' @returns `dwnorm` gives the density,
#' `pwnorm` gives the probability of the wrapped normal distribution function,
#' `rwnorm` generates random deviates (in degrees), and
#' `qwnorm` provides quantiles (in degrees).
#'
#' @importFrom stats rnorm runif
#' @importFrom circular circular pwrappednormal qwrappednormal
#' @name wnorm
#'
#' @seealso [cunif], [wnorm], [wcauchy], and [vonmises]
#'
#' @examples
#' set.seed(1)
#' x <- rwnorm(5, mean = 90, sd = 5)
#'
#' dwnorm(x, mean = 90, sd = 5, axial = FALSE)
#' dwnorm(x, mean = 90, sd = 5, axial = TRUE)
#'
#' pwnorm(x, mean = 90, sd = 5)
#' qwnorm(c(.25, .5, .75), mean = 90, sd = 5)
NULL

#' @rdname wnorm
#' @export
rwnorm <- function(n, mean = 0, sd = 1){
  two_pi <- 2 * pi
  mean_rad <- deg2rad(mean) %% two_pi
  sd_rad <- deg2rad(sd) %% two_pi

  rho <- exp(-sd_rad^2/2)
  stopifnot(rho > 0 | rho < 1)

  if (rho == 0)
    result <- stats::runif(n, 0, two_pi)
  else if (rho == 1)
    result <- rep(mean_rad, n)
  else {
    result <- stats::rnorm(n, mean_rad, sd_rad) %% two_pi
  }

  return(rad2deg(result))
}

# #' @rdname wnorm
# #' @export
# dwnorm <- function(theta, mean = 0, sd = 1, axial = FALSE, ..., log = FALSE){
#   f <- if(isTRUE(axial)) 2 else 1
#   two_pi <- 2 * pi
#   mu <- deg2rad(f * mean) %% two_pi
#   sd_rad <- deg2rad(f * sd) %% two_pi
#   x <- deg2rad(f * theta)  %% two_pi
#
#   d <- circular::dwrappednormal(circular::circular(x), mu = circular::circular(mu), sd = sd_rad, ...)
#
#   if(isTRUE(log)) d <- log(d)
#   return(f * d)
# }

#' @rdname wnorm
#' @export
dwnorm <- function(theta, mean = 0, sd = 1, axial = FALSE, log = FALSE){
  f <- if(isTRUE(axial)) 2 else 1
  two_pi <- 2 * pi
  mu <- deg2rad(f * mean) %% two_pi
  sd_rad <- deg2rad(f * sd) %% two_pi
  x <- deg2rad(f * theta)  %% two_pi
  delta <- x - mu

  k <- -8:8
  dens <- Reduce(`+`, lapply(k, function(kk) exp(-(delta + 2 * pi * kk)^2 / (2 * sd_rad^2))))
  d <- pi / 180 * dens / (sd_rad * sqrt(2 * pi))

  if(isTRUE(log)) d <- log(d)
  return(f * d)
}

#' @rdname wnorm
#' @export
pwnorm <- function(theta, mean = 0, sd = 1, axial = FALSE, from = NULL, ...){
  f <- if(isTRUE(axial)) 2 else 1
  two_pi <- 2 * pi
  mu <- deg2rad(f * mean) %% two_pi
  sd_rad <- deg2rad(f * sd) %% two_pi
  q <- deg2rad(f * theta)  %% two_pi

  circular::pwrappednormal(circular::circular(q), circular::circular(mu), sd = sd_rad, from = from, ...)
}

#' @rdname wnorm
#' @export
qwnorm <- function(p, mean = 0, sd = 1, axial = FALSE, from = NULL, tol = .Machine$double.eps^(0.6), ...){
  f <- if(isTRUE(axial)) 2 else 1
  two_pi <- 2 * pi
  mu <- deg2rad(f * mean) %% two_pi
  sd_rad <- deg2rad(f * sd) %% two_pi

  circular::qwrappednormal(p, circular::circular(mu), sd = sd_rad, from = from, tol = tol, ...) |>
     as.numeric()
}

## Wrapped Cauchy --------------------------------------------------------------

#' The Wrapped Cauchy Distribution
#'
#' Density, probability distribution function, quantiles, and random generation
#' for the circular wrapped Cauchy distribution with mean \eqn{\mu} and rho \eqn{\rho}
#'
#' @inheritParams vonmises
#' @inheritParams cunif
#' @param rho numeric. Concentration parameter in the range (0, 1)
#' @param log.p logical. If `TRUE`, probabilities p are given as log(p).
#'
#' @returns `dwcauchy` gives the density,
#' `pwcauchy` gives the probability of the wrapped Cauchy distribution function,
#' `rwcauchy` generates random deviates (in degrees), and
#' `qrwcauchy` provides quantiles (in degrees).
#'
#' @name wcauchy
#' @seealso [cunif], [wnorm], and [vonmises]
#'
#' @importFrom stats rcauchy
#'
#' @examples
#' set.seed(1)
#' x <- rwcauchy(5, mean = 90, rho = exp(-1))
#'
#' dwcauchy(x, mean = 90, rho = exp(-1))
#' dwcauchy(x, mean = 90, rho = exp(-1), axial = TRUE)
#'
#' pwcauchy(x, mean = 90, rho = exp(-1))
#' qwcauchy(c(.25, .5, .75), mean = 90, rho = exp(-1))
NULL

#' @rdname wcauchy
#' @export
rwcauchy <- function(n, mean, rho){
  two_pi <- 2 * pi
  mu <- deg2rad(mean) %% two_pi

  if (rho == 0)
    result <- runif(n, 0, two_pi)
  else if (rho == 1)
    result <- rep(mu, n)
  else {
    scale <- -log(rho)
    result <- stats::rcauchy(n, mu, scale)%%(two_pi)
  }
  return(rad2deg(result))
}

#' @rdname wcauchy
#' @export
dwcauchy <- function(theta, mean, rho, axial = FALSE, log = FALSE){
  stopifnot(rho >= 0 & rho <= 1)
  f <- if (axial) 2 else 1
  x <- deg2rad((f * theta) %% 360)
  mu <- deg2rad((f * mean) %% 360)
  d <- (1 - rho^2)/((2 * pi) * (1 + rho^2 - 2 * rho * cos(x - mu))) * f * pi / 180
  if(isTRUE(log)) log(d) else d
}

#' @rdname wcauchy
#' @export
pwcauchy <- function(theta, mean, rho, axial = FALSE, from = NULL, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(rho >= 0 & rho <= 1)
  fac <- if (isTRUE(axial)) 2 else 1
  period <- 360 / fac
  mu <- deg2rad((fac * mean) %% 360)
  if (is.null(from)) from <- mean - period / 2
  k <- (1 + rho) / (1 - rho)

  F0 <- function(th_deg) {
    x <- deg2rad((fac * th_deg) %% 360)
    d <- ((x - mu + pi) %% (2 * pi)) - pi   # wrap into [-pi, pi)
    if (isTRUE(all.equal(rho, 1))) return(ifelse(d >= 0, 1, 0))
    0.5 + atan(k * tan(d / 2)) / pi
  }

  p <- (F0(theta) - F0(from)) %% 1
  if (!isTRUE(lower.tail)) p <- 1 - p
  if (isTRUE(log.p)) p <- log(p)
  p
}

#' @rdname wcauchy
#' @export
qwcauchy <- function(p, mean, rho, axial = FALSE, from = NULL, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(rho >= 0 & rho <= 1)
  if (isTRUE(log.p)) p <- exp(p)
  if (!isTRUE(lower.tail)) p <- 1 - p
  bad <- p < 0 | p > 1
  if (any(bad, na.rm = TRUE)) warning("NaNs produced")
  p[bad] <- NaN

  fac <- if (isTRUE(axial)) 2 else 1
  period <- 360 / fac
  mu <- deg2rad((fac * mean) %% 360)
  if (is.null(from)) from <- mean - period / 2
  k <- (1 + rho) / (1 - rho)

  x0 <- deg2rad((fac * from) %% 360)
  d0 <- ((x0 - mu + pi) %% (2 * pi)) - pi
  F0_from <- if (isTRUE(all.equal(rho, 1))) as.numeric(d0 >= 0) else 0.5 + atan(k * tan(d0 / 2)) / pi

  q <- (F0_from + p) %% 1
  d <- if (isTRUE(all.equal(rho, 1))) ifelse(q == 0, 0, NA) else 2 * atan(tan(pi * (q - 0.5)) / k)
  x <- mu + d
  (rad2deg(x) %% 360) / fac
}

# Wrapped Levy------------------------------------------------------------------
.levy_wrap_series_u <- function(u, cc, P) {
  p <- 1:P
  sp <- sqrt(cc * p)
  ang <- outer(u, p, function(uu, pp) sqrt(cc * pp) - pp * uu)
  s <- rowSums(matrix(exp(-sp), nrow = length(u), ncol = P, byrow = TRUE) * cos(ang))
  1 + 2 * s
}

.levy_choose_P <- function(cc, tol = 1e-9, P0 = 100, P.max = 20000) {
  P <- P0
  prev <- .levy_wrap_series_u(0, cc, P)
  repeat {
    P <- P * 2
    if (P > P.max) return(P.max)
    cur <- .levy_wrap_series_u(0, cc, P)
    if (abs(cur - prev) < tol) return(P)
    prev <- cur
  }
}

## the u* in [0, 2*pi) that maximises the wrapped kernel -- computed once per
## bandwidth (cc = f*bw), reused for every theta/data point in that call.
.levy_wrap_offset <- function(cc, tol = 1e-9) {
  P <- .levy_choose_P(cc, tol = tol)
  fn <- function(u) .levy_wrap_series_u(u, cc, P)
  grid <- seq(0, 2 * pi, length.out = 721)[-721]
  vals <- fn(grid)
  i0 <- which.max(vals)
  lo <- grid[max(1, i0 - 2)]
  hi <- grid[min(length(grid), i0 + 2)] + (2 * pi / 720) * 2
  stats::optimize(function(u) -fn(u), interval = c(lo, hi), tol = 1e-10)$minimum
}


.dwlevy_series <- function(theta, mean, c, axial, P, offset) {
  f  <- if (isTRUE(axial)) 2 else 1
  cc <- f * c
  dth <- deg2rad(f * theta) - deg2rad(f * mean) + offset
  p <- 1:P
  sp <- sqrt(cc * p)
  ang <- outer(dth, p, function(d, pp) sqrt(cc * pp) - pp * d)
  s <- rowSums(matrix(exp(-sp), nrow = length(dth), ncol = P, byrow = TRUE) * cos(ang))
  pmax(f * pi / 180 / (2 * pi) * (1 + 2 * s), 0)
}


#' Wrapped Lévy Distribution
#'
#' Probability density function of the wrapped Lévy distribution
#'
#' @inheritParams vonmises
#' @param c numeric. Scale factor of the Lévy distribution. Small values indicate concentrated data near the mean.
#' @param P0,P.max numeric. Arguments Control how the truncated series
#' (the Fourier/characteristic-function sum) is grown and when it gives up.
#' In practice you'll almost never touch either one, as for typical `c` (\eqn{c \le 0.1})
#' convergence happens within the first one or two doublings, well under `P0`.
#'
#' @noRd
#'
#' @seealso [cunif], [wnorm], [wcauchy], and [vonmises]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x <- rwcauchy(5, mean = 90, rho = exp(-1))
#'
#' dwlevy(x, mean = 80, c = 10)
#' }
dwlevy <- function(theta, mean, c, axial = FALSE, tol = 1e-8, P0 = 100, P.max = 5000) {
  stopifnot(all(c > 0))
  f  <- if (isTRUE(axial)) 2 else 1
  cc <- f * c
  offset <- .levy_wrap_offset(cc, tol = tol)
  P <- P0
  prev <- .dwlevy_series(theta, mean, c, axial, P, offset)
  repeat {
    P <- P * 2
    if (P > P.max) {
      warning(sprintf("dwlevy: not converged to tol=%.1e by P.max=%d for c=%.4g", tol, P.max, c))
      return(prev)
    }
    cur <- .dwlevy_series(theta, mean, c, axial, P, offset)
    if (max(abs(cur - prev)) < tol) return(cur)
    prev <- cur
  }
}


# Matrix-vectorized version for .calc_circular_density():
.wrappedlevy_kernel_matrix <- function(delta, bw, f, tol = 1e-8, block = 25, P.max = 5e4) {
  cc <- f * bw
  offset <- .levy_wrap_offset(cc, tol = tol)   # one scalar optimisation per call
  dth <- -delta + offset
  acc <- matrix(0, nrow(dth), ncol(dth))
  p_lo <- 1
  repeat {
    p <- p_lo:(p_lo + block - 1)
    sp <- sqrt(cc * p)
    block_sum <- matrix(0, nrow(dth), ncol(dth))
    for (idx in seq_along(p)) {
      block_sum <- block_sum + exp(-sp[idx]) * cos(sp[idx] - p[idx] * dth)
    }
    acc <- acc + block_sum
    p_lo <- p_lo + block
    if (max(abs(block_sum)) < tol) break
    if (p_lo > P.max) {
      warning(sprintf("wrappedlevy kernel: not converged to tol=%.1e by P.max=%d (bw=%.4g)", tol, P.max, bw))
      break
    }
  }
  pmax(f * pi / 180 / (2 * pi) * (1 + 2 * acc), 0)
}

# Kernel density ---------------------------------------------------------------

.calc_circular_density <- function(x, z, bw, axial,  kernel = c("vonmises", "wrappedcauchy", "wrappednormal", "wrappedlevy"), weights = NULL) {
  kernel <- match.arg(kernel)
  nx <- length(x)
  # if (kernel == "vonmises") {
  # y <- sapply(z, FUN = dvm, mean = x, kappa = kappa, axial = axial, log = FALSE)
  # } else {
  #    #rho <- exp(-bw^2/2)
  #    rho <- kappa
  #    y <- sapply(z, dwcauchy, mean = x, rho = rho, axial = axial, log = FALSE)
  #  }
  # apply(y, 2, sum) / nx
  if (is.null(weights)) weights <- rep.int(1 / nx, nx)

  f <- if (isTRUE(axial)) 2 else 1

  X <- deg2rad(f * x) %% (2 * pi)
  Z <- deg2rad(f * z) %% (2 * pi)
  delta <- outer(X, Z, "-") # nx x n matrix, built once

  y <- switch(kernel,
    "vonmises" = {
      kappa <- bw
      bessel_val <- besselI(kappa, nu = 0, expon.scaled = TRUE) # once, not per z
      f * exp(kappa * (cos(delta) - 1)) / (2 * pi * bessel_val)
    },
   "wrappedcauchy" = {
     rho <- bw
     f * (1 - rho^2) / (2 * pi * (1 + rho^2 - 2 * rho * cos(delta)))
   },
   "wrappednormal" = {
     sd_rad <- deg2rad(f * bw)
     k <- -8:8
     dens <- Reduce(`+`, lapply(k, function(kk) exp(-(delta + 2 * pi * kk)^2 / (2 * sd_rad^2))))
     f * pi / 180 * dens / (sd_rad * sqrt(2 * pi))
   },
   "wrappedlevy" = .wrappedlevy_kernel_matrix(delta, bw, f)
  )

  #colSums(y) / nx
  colSums(y * weights)  # weights already sum to 1 (or intentionally don't, if subdensity)
}


#' Circular Kernel Density Estimation
#'
#' Kernel density estimates for circular data from a given kernel (von Mises,
#' wrapped Cauchy, and wrapped Normal distribution and bandwidth
#'
#' @param x numeric. A vector of angles (in degrees) from which the estimate is to be computed.
#' @param z numeric. Angles where the density is estimated. If `NULL` equally
#' spaced angles are used according to the parameters `from`, `to` and `n`.
#' @param bw,kappa,rho,sd,c numeric. Smoothing bandwidth expressed as the concentration
#' parameter \eqn{\kappa}{k} for the von Mises distribution, \eqn{\rho}{rho} for the
#' wrapped Cauchy distribution, or \eqn{\sigma}{sd} for the wrapped normal distribution.
#' Small and large values for the von Mises and wrapped normal/Cauchy distribution,
#' respectively, gives smooth density lines. If not specified, parameter will be estimated using
#' [est.kappa()]  for the von Mises distribution, or set to \eqn{p \exp(-1)}{p exp(-1)} and `1`
#' for the wrapped Cauchy and wrapped Normal distribution (where \eqn{p = 2}{p=2}
#' when `axial=TRUE` and 1 otherwise), respectively.
#' @param kernel character. The smoothing kernel to be used; one of `"vonmises"`
#' (the default), `"wrappedcauchy"`, `"wrappednormal`, for the von Mises,
#' the wrapped Cauchy, and the wrapped Normal distribution.
#' @param weights numeric. A vector of observation weights, of the same length as `x`,
#' to give individual observations weight in the density estimate. Should sum to 1;
#' a warning is issued if it doesn't (unless `subdensity = TRUE`). Defaults to
#' equal weight `1/length(x)` per observation, matching [stats::density()].
#' @param subdensity logical. If `TRUE`, suppress the "sum(weights) != 1" warning,
#' for when a deliberately partial (sub-)density is wanted, e.g. one group's
#' contribution to a shared total. See [stats::density()].
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param n integer. Number of equally spaced angles at which the density is
#'  to be estimated.
#' @inheritParams circular_plot
#' @inheritParams stats::density
#'
#' @seealso [stats::density()], [dvm()], [dwcauchy()], [dwnorm()], and [plot_density()]
#' @return Object of class `"density"`
#' @export
#'
#' @examples
#' w <-  weighting(san_andreas$unc)
#'
#' # von Mises kernel density
#' circular_density(san_andreas$azi, kappa = 100)
#' circular_density(san_andreas$azi, weights = w, kappa = 100)
#'
#' # wrapped Cauchy kernel density
#' circular_density(san_andreas$azi, rho = 0.9, kernel = "wrappedcauchy")
#'
#' # wrapped Normal kernel density
#' circular_density(san_andreas$azi, sd = 5, kernel = "wrappednormal")
circular_density <- function(x, z = NULL, bw = NULL, weights = NULL, na.rm = TRUE,
                             from = 0, to = 360, n = 512L,
                             axial = TRUE,
                             kappa = NULL, rho = NULL, sd = NULL, c = NULL,
                             kernel = c("vonmises", "wrappedcauchy", "wrappednormal", "wrappedlevy"),
                             adjust = 1, subdensity = FALSE) {
  kernel <- match.arg(kernel)
  f <- if (isTRUE(axial)) 2 else 1
  bw <- if (is.null(bw)) {
    if (kernel == "vonmises") {
      ifelse(is.null(kappa), 10, kappa)
    } else if (kernel == "wrappedcauchy"){
      ifelse(is.null(rho), exp(-1) * f, rho)
    } else if(kernel == "wrappednormal"){
      ifelse(is.null(sd), 1, sd)
    } else {
      ifelse(is.null(c), 5, c)
      #NULL
    }
  } else {
    bw
  }

  if (!is.numeric(x)) stop("argument 'x' must be numeric")
  N <- length(x)
  has.wts <- !is.null(weights)
  if (has.wts && length(weights) != N) {
    stop("'x' and 'weights' have unequal length")
  }

  x.na <- is.na(x)

  if (any(x.na)) {
    if (isTRUE(na.rm)) {
      N <- length(x <- x[!x.na])
      if (has.wts) {
        trueD <- isTRUE(all.equal(1, sum(weights)))
        weights <- weights[!x.na]
        if (trueD) weights <- weights / sum(weights)
      }
    } else {
      stop("'x' contains missing values")
    }
  }
  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
    if (has.wts) weights <- weights[x.finite]
  }
  nx <- length(x)

  if (!has.wts) {
    weights <- rep.int(1 / nx, nx)
  } else {
    if (!all(is.finite(weights))) stop("'weights' must all be finite")
    if (any(weights < 0)) stop("'weights' must not be negative")
    if (!subdensity && !isTRUE(all.equal(1, sum(weights)))) {
      warning("sum(weights) != 1  -- will not get true density")
    }
  }

  if (is.null(z)) {
    z <- seq(from = from, to = to, length = n)
  } else {
    if (!is.numeric(z)) {
      stop("argument 'z' must be numeric")
    }
    #namez <- deparse(substitute(z))
    z.na <- is.na(z)
    if (any(z.na)) {
      if (isTRUE(na.rm)) {
        z <- z[!z.na]
      } else {
        stop("z contains missing values")
      }
    }
    z.finite <- is.finite(z)
    if (any(!z.finite)) {
      z <- z[z.finite]
    }
  }

  d <- .calc_circular_density(x, z, bw = bw, axial = axial, kernel = kernel, weights = weights)

  structure(
    list(
      x = z,
      y = d * adjust,
      bw = bw,
      n = N,
      call = match.call(),
      data.name = deparse1(substitute(x)),
      has.na = any(x)
    ),
    class = "density"
  )
}
