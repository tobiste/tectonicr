# Distribution ####

## Circular Uniform ------------------------------------------------------------

#' The Circular Uniform Distribution
#'
#' Density, probability distribution function, quantiles, and random generation
#' for the circular uniform distribution.
#'
#' @param n integer. Number of observations in degrees
#' @param p numeric. Vector of probabilities with values in \eqn{[0,1]}{[0,1]}.
#' @param theta numeric. Angular value in degrees
#' @param log,log.p logical. If `TRUE`, probabilities p are given as log(p).
#' @param lower.tail logical. If `TRUE` (default), probabilities are \eqn{P(\Theta \le \theta)},
#' otherwise \eqn{P(\Theta > \theta)}.
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @returns `dcunif` gives the density,
#' `pcunif` gives the probability of circular uniform distribution function,
#' `rcunif` generates random deviates (in degrees), and
#' `qcunif` provides quantiles (in degrees).
#'
#' @seealso [wcauchy] and [vonmises]
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

# pvm.mu0 <- function(theta, kappa, acc) {
#   flag <- TRUE
#   p <- 1
#   sum <- 0
#   while (flag) {
#     term <- (besselI(x = kappa, nu = p, expon.scaled = FALSE) *
#       sin(p * theta)) / p
#     sum <- sum + term
#     p <- p + 1
#     if (abs(term) < acc) {
#       flag <- FALSE
#     }
#   }
#   theta / (2 * pi) + sum / (pi * besselI(
#     x = kappa, nu = 0,
#     expon.scaled = FALSE
#   ))
# }



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
#' @seealso [wcauchy] and [cunif]
#'
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
  if (axial) {
    x <- circular::circular(theta, units = "degrees", modulo = "pi")
    mu <- circular::circular(mean, units = "degrees", modulo = "pi")
    d <- circular::daxialvonmises(x, mu, kappa)
    if (log) d <- log(d)
    return(d)
  } else {
    two_pi <- 2 * pi
    x <- deg2rad(theta) %% two_pi
    mu <- deg2rad(mean) %% two_pi

    delta <- x - mu
    delta_mod <- delta %% two_pi

    n <- length(x)

    if (log) {
      if (kappa == 0) {
        vm <- rep(-log(two_pi), n)
      } else if (kappa < 1e+05) {
        log_bessel <- log(besselI(kappa, nu = 0, expon.scaled = TRUE))
        vm <- -(log(two_pi) + log_bessel + kappa) + kappa * cos(delta)
      } else {
        vm <- ifelse(delta_mod == 0, Inf, -Inf)
      }
    } else {
      if (kappa == 0) {
        vm <- rep(1 / two_pi, n)
      } else if (kappa < 1e+05) {
        bessel_val <- besselI(kappa, nu = 0, expon.scaled = TRUE)
        vm <- (1 / (two_pi * bessel_val)) * (exp(cos(delta) - 1))^kappa
      } else {
        vm <- ifelse(delta_mod == 0, Inf, 0)
      }
    }
    return(vm)
  }
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

  circular::qvonmises(p, mu, kappa, from, tol = tol, ...) |> as.numeric()
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
#' @seealso [vonmises] and [cunif]
#'
#' @importFrom circular circular rwrappedcauchy
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
  mu <- circular::circular(mean, units = "degrees", modulo = "2pi")
  circular::rwrappedcauchy(n, mu, rho) |> as.numeric()
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

# Kernel density ---------------------------------------------------------------

.calc_circular_density <- function(x, z, bw, axial, kernel = c("vonmises", "wrappedcauchy")) {
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

  f <- if (isTRUE(axial)) 2 else 1

  X <- deg2rad(f * x) %% (2 * pi)
  Z <- deg2rad(f * z) %% (2 * pi)
  delta <- outer(X, Z, "-") # nx x n matrix, built once

  if (kernel == "vonmises") {
    kappa <- bw
    # kappa <- ifelse(is.null(bw), 10, bw)
    bessel_val <- besselI(kappa, nu = 0, expon.scaled = TRUE) # once, not per z
    y <- f * exp(kappa * (cos(delta) - 1)) / (2 * pi * bessel_val)
  } else {
    rho <- bw
    # rho <- exp(-bw^2/2)
    # rho <- ifelse(is.null(bw), exp(-1) * f, bw)
    y <- f * (1 - rho^2) / (2 * pi * (1 + rho^2 - 2 * rho * cos(delta)))
  }
  colSums(y) / nx
}


#' Circular Kernel Density
#'
#' Multiples of a von Mises density or wrapped Cauchy distribution for circular data
#'
#' @param x numeric. A vector of angles (in degrees) from which the estimate is to be computed.
#' @param z numeric. Angles where the density is estimated. If `NULL` equally
#' spaced angles are used according to the parameters `from`, `to` and `n`.
#' @param bw,kappa,rho numeric. Smoothing bandwidth expressed as the concentration
#' parameter \eqn{\kappa}{k} for the von Mises distribution or \eqn{\rho}{rho} for the
#' wrapped Cauchy distribution.
#' Small and large values for the von Mises and wrapped Cauchy distribution,
#' respectively, gives smooth density lines. If not specified, parameter will be estimated using
#' [est.kappa()]  for the von Mises distribution, or set to \eqn{p \exp(-1)}{p exp(-1)}
#' for the wrapped Cauchy distribution (where \eqn{p = 2}{p=2} when `axial=TRUE` and 1 otherwise).
#' @param kernel character. The smoothing kernel to be used; one of `"vonmises"`
#' (the default) or `"wrappedcauchy"` for the von Mises or the Wrapped Cauchy
#' distribution.
#' @param axial Logical. Whether data are uniaxial (`axial=FALSE`)
#' or biaxial (`TRUE`, the default).
#' @param n integer. Number of equally spaced angles at which the density is
#'  to be estimated.
#' @inheritParams circular_plot
#' @inheritParams stats::density
#'
#' @seealso [stats::density()], [dvm()] and [dwcauchy()]
#' @return Object of class `"density"`
#' @export
#'
#' @examples
#' # von Mises kernel density
#' circular_density(san_andreas$azi, kappa = 100)
#'
#' # wrapped Cauchy kernel density
#' circular_density(san_andreas$azi, rho = 0.9, kernel = "wrappedcauchy")
circular_density <- function(x, z = NULL, bw = NULL, na.rm = TRUE, from = 0,
                             to = 360, n = 512, axial = TRUE,
                             kappa = NULL, rho = NULL, kernel = c("vonmises", "wrappedcauchy"),
                             adjust = 1) {
  kernel <- match.arg(kernel)
  f <- if (isTRUE(axial)) 2 else 1

  bw <- if (is.null(bw)) {
    if (kernel == "vonmises") {
      ifelse(is.null(kappa), 10, kappa)
    } else {
      ifelse(is.null(rho), exp(-1) * f, rho)
    }
  } else {
    bw
  }


  if (is.null(z)) {
    z <- seq(from = from, to = to, length = n)
  } else {
    if (!is.numeric(z)) {
      stop("argument 'z' must be numeric")
    }
    namez <- deparse(substitute(z))
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

  d <- .calc_circular_density(x, z, bw = bw, axial = axial, kernel = kernel)

  structure(
    list(
      x = z,
      y = d * adjust,
      bw = bw,
      n = length(na.omit(z)),
      call = match.call(),
      data.name = deparse1(substitute(x)),
      has.na = any(is.na(x))
    ),
    class = "density"
  )
}
