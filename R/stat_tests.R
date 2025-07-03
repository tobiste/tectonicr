# Tests ####
#' @keywords internal
nchisq_eq <- function(obs, prd, unc) {
  if (is.na(unc) || unc == 0) unc <- 1 # uncertainty cannot be 0
  w <- deviation_norm(obs, prd)
  x <- (w / unc)^2
  y <- (90 / unc)^2
  return(c(x, y))
}

#' Normalized Chi-Squared Test for Circular Data
#'
#' A quantitative comparison between the predicted and observed directions of
#' \eqn{\sigma_{Hmax}}{SHmax} is obtained by the calculation of the average
#' azimuth and by a normalized \eqn{\chi^2}{chi-squared} test.
#'
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics. *Journal of Geophysical Research: Solid Earth*, **103**,
#'   5037-5059, doi: 10.1029/97JB03390.
#' @param prd Numeric vector containing the modeled azimuths of
#' \eqn{\sigma_{Hmax}}{SHmax}, i.e.
#' the return object from \code{model_shmax()}
#' @param obs Numeric vector containing the observed azimuth of
#' \eqn{\sigma_{Hmax}}{SHmax},
#' same length as \code{prd}
#' @param unc Uncertainty of observed \eqn{\sigma_{Hmax}}{SHmax}, either a
#' numeric vector or a number
#'
#' @returns Numeric vector
#'
#' @details
#' The normalized \eqn{\chi^2}{chi-squared} test is
#' \deqn{ {Norm} \chi^2_i =
#'  = \frac{
#'    \sum^M_{i = 1} \left( \frac{\alpha_i - \alpha_{{predict}}}{\sigma_i}
#'    \right) ^2}
#'    {\sum^M_{i = 1} \left( \frac{90}{\sigma_i} \right) ^2 }}{
#'    (sum( ((obs-prd)/unc)^2 ) / sum( (90/unc)^2 )
#'    }
#' The value of the chi-squared test statistic is a number between 0 and 1
#' indicating the quality of the predicted \eqn{\sigma_{Hmax}}{SHmax}
#' directions. Low values
#' (\eqn{\le 0.15}) indicate good agreement,
#' high values (\eqn{> 0.7}) indicate a systematic misfit between predicted and
#' observed \eqn{\sigma_{Hmax}}{SHmax} directions.
#'
#' @export
#'
#' @examples
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to
#' # Pacific plate
#' data(san_andreas)
#' point <- data.frame(lat = 45, lon = 20)
#' prd <- model_shmax(point, PoR)
#' norm_chisq(obs = c(50, 40, 42), prd = prd$sc, unc = c(10, NA, 5))
#'
#' data(san_andreas)
#' prd2 <- PoR_shmax(san_andreas, PoR, type = "right")
#' norm_chisq(obs = prd2$azi.PoR, 135, unc = san_andreas$unc)
norm_chisq <- function(obs, prd, unc) {
  N <- length(obs)
  if (length(prd) == 1) {
    prd <- rep(prd, N)
  }

  if (length(unc) == 1) {
    unc <- rep(unc, N)
  }
  stopifnot(length(prd) == N, length(unc) == N)

  keep <- !is.na(obs) & !is.na(prd) & !is.na(unc)
  obs <- obs[keep]
  prd <- prd[keep]
  unc <- unc[keep]

  xy <- mapply(
    FUN = nchisq_eq,
    obs = obs, prd = prd, unc = unc
  )

  sum(xy[1, ], na.rm = FALSE) / sum(xy[2, ], na.rm = FALSE)
}


#' Rayleigh Test of Circular Uniformity
#'
#' Performs a Rayleigh test for uniformity of circular/directional data by
#' assessing the significance of the mean resultant length.
#'
#' @param x numeric vector. Values in degrees
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param mu (optional) The specified or known mean direction (in degrees) in
#' alternative hypothesis
#' @param quiet logical. Prints the test's decision.
#'
#' @details \describe{
#' \item{\eqn{H_0}{H0}:}{angles are randomly distributed around the circle.}
#' \item{\eqn{H_1}{H1}:}{angles are from non-uniformly distribution with unknown mean
#' direction and mean resultant length (when `mu` is `NULL`. Alternatively (when
#' `mu` is specified),
#' angles are non-uniformly distributed around a specified direction.}
#' }
#' If `statistic > p.value`, the null hypothesis is rejected,
#' i.e. the length of the mean resultant differs significantly from zero, and
#' the angles are not randomly distributed.
#'
#' @note Although the Rayleigh test is consistent against (non-uniform)
#' von Mises alternatives, it is not consistent against alternatives with
#' `p = 0` (in particular, distributions with antipodal symmetry, i.e. axial
#' data). Tests of non-uniformity which are consistent against all alternatives
#' include Kuiper's test ([kuiper_test()]) and Watson's \eqn{U^2} test
#' ([watson_test()]).
#'
#' @returns a list with the components:
#' \describe{
#'  \item{`R` or `C`}{mean resultant length or the dispersion (if `mu` is
#'  specified). Small values of `R` (large values of `C`) will reject
#'  uniformity. Negative values of `C` indicate that vectors point in opposite
#'  directions (also lead to rejection).}
#'  \item{`statistic`}{test statistic}
#'  \item{`p.value`}{significance level of the test statistic}
#' }
#'
#' @references
#' Fisher, N. I. (1993) Statistical Analysis of Circular Data, Cambridge
#' University Press.
#'
#' @seealso [mean_resultant_length()], [circular_mean()], [norm_chisq()],
#' [kuiper_test()], [watson_test()], [weighted_rayleigh()]
#'
#' @export
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' rayleigh_test(pidgeon_homing, axial = FALSE) # Do not reject null hypothesis.
#' # R = 0.22; stat = 0.497, p = 0.62
#'
#' # Example data from Davis (1986), pp. 316
#' finland_striae <- c(
#'   23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
#'   113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
#'   132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
#'   165, 171, 172, 179, 181, 186, 190, 212
#' )
#' rayleigh_test(finland_striae, axial = FALSE) # reject null hypothesis
#' rayleigh_test(finland_striae, mu = 105, axial = FALSE) # reject null hypothesis
#'
#' # Example data from Mardia and Jupp (1999), pp. 99
#' atomic_weight <- c(
#'   rep(0, 12), rep(3.6, 1), rep(36, 6), rep(72, 1),
#'   rep(108, 2), rep(169.2, 1), rep(324, 1)
#' )
#' rayleigh_test(atomic_weight, 0, axial = FALSE) # reject null hypothesis
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' rayleigh_test(san_andreas$azi) # reject null hypothesis
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' rayleigh_test(sa.por$azi.PoR, mu = 135) # reject null hypothesis
rayleigh_test <- function(x, mu = NULL, axial = TRUE, quiet = FALSE) {
  f <- if (axial) 2 else 1

  if (is.null(mu)) {
    x <- x[!is.na(x)]
    xf <- (x * f) %% 360
    n <- length(x)

    R <- mean_resultant_length(xf, na.rm = FALSE)
    S <- 2 * n * R^2
    Z <- S/2
    # S_mod <- (1 - 1 / (2 * n)) * S + (n * R^4) / 2
    # if(n <= 10){
    #  p.value <- p_value3(R, n)
    # } else  {
    p.value <- rayleigh_p_value1(Z, n)
    # }
    # p.value2 <- rayleigh_p_value1(S_mod /2 , n)

    result <- list(
      R = R,
      statistic = Z ,
      # statistic_mod = S_mod,
      p.value = p.value#,
      # p.value_mod = p.value2
    )
    if (!quiet) {
      if (result$statistic >= p.value) {
        message("Reject Null Hypothesis\n")
      } else {
        message("Do Not Reject Null Hypothesis\n")
      }
    }
  } else {

    #remove NA's
    keep <- !is.na(x) & !is.na(mu)
    x <- x[keep]
    # mu <- mu[keep]

    x <- x * f
    mu <- mu * f
    xmu <- x - mu
    n <- length(x)

    C <- (sum(cosd(xmu))) / n
    s <- sqrt(2 * n) * C
    p.value <- rayleigh_p_value2(s, n)

    result <- list(
      C = C,
      statistic = s,
      p.value = p.value
    )
    if (!quiet) {
      if (s >= p.value) {
        message("Reject Null Hypothesis\n")
      } else {
        message("Do Not Reject Null Hypothesis\n")
      }
    }
  }

  return(result)
}

#' @keywords internal
rayleigh_p_value1 <- function(K, n, wilkie = FALSE) {
  if (!wilkie) {
    # Pearson. 1906; Greenwood and Durand, 1955
    P <- exp(-K)
    if (n < 50) {
      temp <- 1 +
        (2 * K - K^2) / (4 * n) -
        (24 * K - 132 * K^2 + 76 * K^3 - 9 * K^4) / (288 * n^2)
    } else {
      temp <- 1
    }
    P * temp
    min(max(P * temp, 0), 1)
  } else {
    # Wilkie 1983
    Rn <- K * n
    temp <- sqrt(1 + 4 * n + 4 * (n^2 - Rn^2)) - (1 + 2 * n)
    round(exp(temp), 3)
  }
}

#' @keywords internal
rayleigh_p_value2 <- function(K, n) {
  # Greenwood and Durand, 1957
  pK <- stats::pnorm(K) # distribution function of standard normal distribution
  fK <- stats::dnorm(K) # density function of standard normal distribution
  P <- 1 - pK + fK * (
    (3 * K - K^3) / (16 * n) +
      (15 * K + 305 * K^3 - 125 * K^5 + 9 * K^7) / (4608 * n^2)
  )
  min(max(P, 0), 1)
}

#' Weighted Goodness-of-fit Test for Circular Data
#'
#' Weighted version of the Rayleigh test (or V0-test) for uniformity against a
#' distribution with a priori expected von Mises concentration.
#' @param x numeric vector. Values in degrees
#' @param w numeric vector weights of length `length(x)`. If `NULL`, the
#' non-weighted Rayleigh test is performed.
#' @param mu The *a priori* expected direction (in degrees) for the alternative
#' hypothesis.
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param quiet logical. Prints the test's decision.
#'
#' @details
#' The Null hypothesis is uniformity (randomness). The alternative is a
#' distribution with a (specified) mean direction (`mu`).
#' If `statistic >= p.value`, the null hypothesis of randomness is rejected and
#' angles derive from a distribution with a (or the specified) mean direction.
#'
#' @returns a list with the components:
#' \describe{
#'  \item{`R` or `C`}{mean resultant length or the dispersion (if `mu` is
#'  specified). Small values of `R` (large values of `C`) will reject
#'  uniformity. Negative values of `C` indicate that vectors point in opposite
#'  directions (also lead to rejection).}
#'  \item{`statistic`}{Test statistic}
#'  \item{`p.value`}{significance level of the test statistic}
#' }
#'
#' @seealso [rayleigh_test()]
#'
#' @export
#'
#' @examples
#' # Load data
#' data("cpm_models")
#' data(san_andreas)
#' PoR <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "na", "pa")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' data("iceland")
#' PoR.ice <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "eu", "na")
#' ice.por <- PoR_shmax(iceland, PoR.ice, "out")
#' data("tibet")
#' PoR.tib <- equivalent_rotation(cpm_models[["NNR-MORVEL56"]], "eu", "in")
#' tibet.por <- PoR_shmax(tibet, PoR.tib, "in")
#'
#' # GOF test:
#' weighted_rayleigh(tibet.por$azi.PoR, mu = 90, w = 1 / tibet$unc)
#' weighted_rayleigh(ice.por$azi.PoR, mu = 0, w = 1 / iceland$unc)
#' weighted_rayleigh(sa.por$azi.PoR, mu = 135, w = 1 / san_andreas$unc)
weighted_rayleigh <- function(x, mu = NULL, w = NULL, axial = TRUE, quiet = FALSE) {
  if (is.null(w)) {
    rayleigh_test(x, mu = mu, axial = axial)
  } else {
    #remove NA's
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]

    Z <- sum(w)
    n <- length(w)

    if (is.null(mu)) mu <- circular_mean(x, w, axial, na.rm = FALSE)

    d <- x - mu
    f <- if (axial) 2 else 1

    m <- mean_SC(f * d, w = w, na.rm = FALSE)
    C <- as.numeric(m["C"])
    s <- sqrt(2 * n) * C
    p.value <- rayleigh_p_value2(s, n)

    result <- list(
      C = C,
      statistic = s,
      p.value = p.value
    )
    if (!quiet) {
      if (s >= p.value) {
        message("Reject Null Hypothesis\n")
      } else {
        message("Do Not Reject Null Hypothesis\n")
      }
    }
    return(result)
  }
}

#' Kuiper Test of Circular Uniformity
#'
#' Kuiper's test statistic is a rotation-invariant Kolmogorov-type test statistic.
#' The critical values of a modified Kuiper's test statistic are used according
#' to the tabulation given in Stephens (1970).
#'
#' @param x numeric vector containing the circular data which are expressed in degrees
#' @param alpha Significance level of the test. Valid levels are `0.01`, `0.05`, and `0.1`.
#' This argument may be omitted (`NULL`, the default), in which case, a range for the p-value will be returned.
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or circular, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @returns list containing the test statistic `statistic` and the significance
#' level `p.value`.
#' @param quiet logical. Prints the test's decision.
#'
#' @details
#'
#' If `statistic > p.value`, the null hypothesis is rejected.
#' If not, randomness (uniform distribution) cannot be excluded.
#'
#' @export
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' kuiper_test(pidgeon_homing, alpha = .05)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' kuiper_test(sa.por$azi.PoR, alpha = .05)
kuiper_test <- function(x, alpha = 0, axial = TRUE, quiet = FALSE) {
  allowed_alphas <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15)
  if (!(alpha %in% allowed_alphas)) {
    stop("'alpha' must be one of: 0, 0.01, 0.025, 0.05, 0.1, 0.15")
  }

  thresholds <- c(1.537, 1.62, 1.747, 1.862, 2.001)

  kuiper.crits <- cbind(
    rev(allowed_alphas[-1]),
    thresholds
  )
  f <- if (axial) 2 else 1

  x <- x[!is.na(x)] # remove NA's
  x <- (x * f) %% 360
  u <- sort(deg2rad(x) %% (2 * pi)) / (2 * pi)
  n <- length(x)
  i <- 1:n
  D.P <- max(i / n - u)
  D.M <- max(u - (i - 1) / n)
  sqrt_n <- sqrt(n)
  V <- D.P + D.M
  V <- V * (sqrt_n + 0.155 + 0.24 / sqrt_n)

  if (alpha == 0) {
    labels <- c(
      "P-value > 0.15",
      "0.10 < P-value < 0.15",
      "0.05 < P-value < 0.10",
      "0.025 < P-value < 0.05",
      "0.01 < P-value < 0.025",
      "P-value < 0.01"
    )
    idx <- findInterval(V, thresholds) + 1
    p.value <- labels[idx]
  } else {
    idx <- which(alpha == kuiper.crits[, 1])
    p.val.threshold <- kuiper.crits[idx, 2]
    p.value <- p.val.threshold

    if (!quiet) {
      msg <- if (V > p.val.threshold) "Reject Null Hypothesis\n" else "Do Not Reject Null Hypothesis\n"
      message(msg)
    }
  }
  return(
    list(
      statistic = V,
      p.value = unname(p.value)
    )
  )
}

#' Watson's \eqn{U^2} Test of Circular Uniformity
#'
#' Watson's test statistic is a rotation-invariant Cramer - von Mises test
#'
#' @param x numeric vector. Values in degrees
#' @param alpha Significance level of the test. Valid levels are `0.01`, `0.05`,
#' and `0.1`.
#' This argument may be omitted (`NULL`, the default), in which case, a range
#' for the p-value will be returned.
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}-periodical
#' (`TRUE`, the default) or circular, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#' @param mu (optional) The specified mean direction (in degrees) in alternative
#'  hypothesis
#' @param dist Distribution to test for. The default, `"uniform"`, is the
#' uniform distribution. `"vonmises"` tests the von Mises distribution.
#' @param quiet logical. Prints the test's decision.
#'
#' @returns list containing the test statistic `statistic` and the significance
#' level `p.value`.
#'
#' @details
#' If `statistic > p.value`, the null hypothesis is rejected.
#' If not, randomness (uniform distribution) cannot be excluded.
#'
#' @references Mardia and Jupp (1999). Directional Statistics. John Wiley and
#' Sons.
#'
#' @export
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' watson_test(pidgeon_homing, alpha = .05)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' watson_test(sa.por$azi.PoR, alpha = .05)
#' watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
watson_test <- function(x, alpha = 0, dist = c("uniform", "vonmises"), axial = TRUE, mu = NULL, quiet = FALSE) {
  allowed_alphas <- c(0, 0.01, 0.025, 0.05, 0.1)
  if (!(alpha %in% allowed_alphas)) {
    stop("'alpha' must be one of: 0, 0.01, 0.025, 0.05, 0.1")
  }

  dist <- match.arg(dist)
  x <- x[!is.na(x)]
  n <- length(x)

  if (dist == "uniform") {
    # if (axial) {
    #   f <- 2
    # } else {
    #   f <- 1
    # }
    f <- 1
    x <- (x * f) %% 360

    # U2 Statistic:
    u <- sort(deg2rad(x)) / (2 * pi)
    if (is.null(mu)) {
      u.bar <- mean(u)
    } else {
      u.bar <- deg2rad(mu %% 360) / (2 * pi)
    }
    i <- seq_len(n)
    u2 <- sum((u - u.bar - (i - .5) / n + .5)^2) + 1 / (12 * n)
    statistic <- (u2 - 0.1 / n + 0.1 / (n^2)) * (1 + 0.8 / n)

    # P-value:
    crits <- c(99, 0.267, 0.221, 0.187, 0.152)
    if (n < 8) {
      p.value <- NA
      warning("Total Sample Size < 8:  Results are not valid")
    }

    if (alpha == 0) {
      thresholds <- c(rev(crits[-1]), Inf)

      messages <- c(
        "P-value > 0.10",
        "0.05 < P-value < 0.10",
        "0.025 < P-value < 0.05",
        "0.01 < P-value < 0.025",
        "P-value < 0.01"
      )
      p.value <- messages[findInterval(statistic, thresholds)]
      #
    } else {
      index <- (1:5)[alpha == allowed_alphas[1:5]]
      p.value <- crits[index]
      if (!quiet) {
        if (statistic > p.value) {
          message("Reject Null Hypothesis\n")
        } else {
          message("Do Not Reject Null Hypothesis\n")
        }
      }
    }
  } else {
    u2_crits <- cbind(
      c(0, 0.5, 1, 1.5, 2, 4, 100),
      c(0.052, 0.056, 0.066, 0.077, 0.084, 0.093, 0.096),
      c(0.061, 0.066, 0.079, 0.092, 0.101, 0.113, 0.117),
      c(0.081, 0.09, 0.11, 0.128, 0.142, 0.158, 0.164)
    )

    if (is.null(mu)) mu <- circular_mean(x, axial = axial, na.rm = FALSE)

    kappa_mle <- est.kappa(x, axial = axial)
    x <- x - mu
    x <- matrix(x, ncol = 1)
    z <- apply(x, 1, pvm, 0, kappa_mle)
    z <- sort(z)
    z.bar <- mean(z)
    i <- seq_len(n)
    sum.terms <- (z - (2 * i - 1) / (2 * n))^2
    statistic <- sum(sum.terms) - n * (z.bar - 0.5)^2 + 1 / (12 * n)

    row <- findInterval(kappa_mle, c(0.25, 0.75, 1.25, 1.75, 3, 5)) + 1

    if (alpha != 0) {
      col <- match(alpha, c(0.1, 0.05, 0.01)) + 1
      if (is.na(col)) {
        stop("'alpha' must be one of: 0.1, 0.05, 0.01")
      }

      p.value <- u2_crits[row, col]

      if (!quiet) {
        message(if (statistic > p.value) "Reject Null Hypothesis" else "Do Not Reject Null Hypothesis")
      }

    } else {
      breaks <- u2_crits[row, 2:4]
      labels <- c(
        "P-value > 0.10",
        "0.05 < P-value < 0.10",
        "0.01 < P-value < 0.05",
        "P-value < 0.01"
      )

      p.value <- labels[findInterval(statistic, breaks, rightmost.closed = TRUE) + 1]
    }
  }
  list(
    statistic = statistic,
    p.value = p.value
  )
}


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
        vm <- - (log(two_pi) + log_bessel + kappa) + kappa * cos(delta)
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
#' Computes the maximum likelihood estimate of \eqn{\kappa}, the concentration
#' parameter of a von Mises distribution, given a set of angular measurements.
#'
#' @param x numeric. angles in degrees
#' @param w numeric. weightings
#' @param bias logical parameter determining whether a bias correction is used
#' in the computation of the MLE. Default for bias is `FALSE` for no bias
#' correction.
#' @param axial logical. Whether the data are axial, i.e. pi-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}-periodical (`FALSE`).
#'
#' @returns numeric.
#' @export
#'
#' @examples
#' set.seed(123)
#' est.kappa(rvm(100, 90, 10), w = 1 / runif(100, 0, 10))
est.kappa <- function(x, w = NULL, bias = FALSE, axial = TRUE) {
  # Default weights
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)
  }

  f <- if(axial) 2 else 1
  x <- (x*f) %% 360

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
