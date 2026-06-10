# Distribution ####
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
#' for the circular normal distribution with mean and kappa.
#'
#' @param n integer. Number of observations in degrees
#' @param p numeric. Vector of probabilities with values in \eqn{[0,1]}{[0,1]}.
#' @param mean numeric. Mean angle in degrees
#' @param kappa numeric. Concentration parameter in the range (0, Inf]
#' @param theta numeric. Angular value in degrees
#' @param from if `NULL` is set to \eqn{\mu-\pi}{mu-pi}. This is the value from
#' which the pvm and qvm are evaluated. in degrees.
#' @param tol numeric. The precision in evaluating the distribution function or the quantile.
#' @param log logical. If `TRUE`, probabilities p are given as log(p).
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param ... parameters passed to [stats::integrate()].
#'
#' @returns `dvm` gives the density,
#' `pvm` gives the probability of the von Mises distribution function,
#' `rvm` generates random deviates (in degrees), and
#' `qvm` provides quantiles (in degrees).
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
dvm <- function(theta, mean, kappa, log = FALSE, axial = FALSE) {
  if (axial) {
    x <- circular::circular(theta, units = "degrees", modulo = "pi")
    mu <- circular::circular(mean, units = "degrees", modulo = "pi")
    d <- circular::daxialvonmises(x, mu, kappa)
    if (log) d <- log(d)
    return(d)
  } else {
    # x <- circular::circular(theta, units = "degrees", modulo = "2pi")
    # mu <- circular::circular(mean, units = "degrees", modulo = "2pi")
    # circular::dvonmises(x, mu = mu, kappa = kappa, log = log)

    two_pi <- 2 * pi
    x <- deg2rad(theta) %% two_pi
    mu <- deg2rad(mean) %% two_pi

    delta <- x - mu
    delta_mod <- delta %% two_pi

    n <- length(x)

    # stopifnot(length(mu==1))
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

  circular::pvonmises(theta, mu, kappa, from = NULL, tol = tol)
}

#' @rdname vonmises
#' @export
qvm <- function(p, mean = 0, kappa, from = NULL, tol = .Machine$double.eps^(0.6), ...) {
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
#' @param x numeric. angles in degrees
#' @param w numeric. weightings
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
