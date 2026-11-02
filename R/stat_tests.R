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
#' \eqn{\sigma_\text{Hmax}}{SHmax} is obtained by the calculation of the average
#' azimuth and by a normalized \eqn{\chi^2}{chi-squared} test (Wdowsinki, 1998)
#'
#' @references Wdowinski, S., 1998, A theory of intraplate
#'   tectonics. *Journal of Geophysical Research: Solid Earth*, **103**,
#'   5037-5059, doi: 10.1029/97JB03390.
#'
#' @param prd Numeric vector containing the modeled azimuths of
#' \eqn{\sigma_\text{Hmax}}{SHmax}, i.e.
#' the return object from \code{model_shmax()}
#' @param obs Numeric vector containing the observed azimuth of
#' \eqn{\sigma_\text{Hmax}}{SHmax},
#' same length as \code{prd}
#' @param unc Uncertainty of observed \eqn{\sigma_\text{Hmax}}{SHmax}, either a
#' numeric vector or a number
#'
#' @returns Numeric vector
#'
#' @details
#' The normalized \eqn{\chi^2}{chi-squared} test is
#' \deqn{ \text{Norm} \chi^2_i = \frac{
#'    \sum^M_{i = 1} \left( \frac{\alpha_i - \alpha_{{predict}}}{\sigma_i}
#'    \right) ^2}
#'    {\sum^M_{i = 1} \left( \frac{90}{\sigma_i} \right) ^2 }}{
#'    (sum( ((obs-prd)/unc)^2 ) / sum( (90/unc)^2 )
#'    }
#' The value of the chi-squared test statistic is a number between 0 and 1
#' indicating the quality of the predicted \eqn{\sigma_\text{Hmax}}{SHmax}
#' directions. Low values
#' (\eqn{\le 0.15}) indicate good agreement,
#' high values (\eqn{> 0.7}) indicate a systematic misfit between predicted and
#' observed \eqn{\sigma_\text{Hmax}}{SHmax} directions.
#'
#' @family One-Sample-Tests
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
#' A test to determine whether a sample of circular or directional data is
#' evenly spread out or clustered around a single specific direction.
#' The test assesses the significance of the mean resultant length.
#' `rayleight_test_perm()` uses permutation to estimate p-values.
#'
#' @param x numeric vector. Values in degrees
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}{pi}-periodical
#' (`TRUE`, the default) or directional, i.e. \eqn{2 \pi}{2pi}-periodical (`FALSE`).
#' In case of axial data, the angles will be doubled for the test.
#' @param mu (optional) The specified or known mean direction (in degrees) in
#' alternative hypothesis.
#' @param quiet logical. Prints the test's decision.
#' @param n_perm integer. Number of permutations.
#' @param alpha Significance level of the test ranging between 0 and 1.
# #' Valid levels are `0.01`, `0.05`, and `0.1`.
# #' This argument may be omitted (`NULL`, the default), in which case, a range for the p-value will be returned.
#'
#' @details
#' ## Hypotheses
#' **Null Hypothesis** \eqn{H_0}{H0}: The population is distributed uniformly (randomly)
#' around the circle with no preferred direction.
#'
#' **Alternative Hypothesis** \eqn{H_1}{H1}: The population is not uniform and has a
#' unimodal (single-peaked) concentration in a preferred direction. When `mu`
#' is specified), angles are non-uniformly distributed around the specified direction.
#'
#'
#' *Mean Resultant Length* (\eqn{\bar{R}}{R̄} or R): A value between 0 and 1 that
#'  measures how concentrated the data points are.
#'
#'  * \eqn{\bar{R}}{R̄̄} = 0: The data is completely spread out around the circle
#'  * \eqn{\bar{R}}{R̄} = 1: All data points point in the exact same direction.
#'
#'  *p-value* The probability of seeing data this clustered purely by chance under the assumption of uniformity.
#'
#'  ## Interpretation
#'
#'  * Small p-value (p < 0.05): Reject the null hypothesis. The length of the
#'  mean resultant differs significantly from zero, and
#' the angles are not randomly distributed. You have strong evidence that the
#' data points share a significant preferred or mean direction (unimodal clustering).
#'  * Large p-value (p \eqn{\ge}{>=} 0.05): Fail to reject the null hypothesis. There is not
#'  enough evidence to claim a preferred direction, meaning the data looks random or uniform around the circle.
#'
#'
#' @note Although the Rayleigh test is consistent against (non-uniform)
#' von Mises alternatives, it is not consistent against alternatives with
#' `p = 0` (in particular, distributions with antipodal symmetry, i.e. axial
#' data). Tests of non-uniformity which are consistent against all alternatives
#' include Kuiper's test ([kuiper_test()]) and Watson's \eqn{U^2}{U2} test
#' ([watson_test()]).
#'
#' ## Limitations
#'
#'  * The test assumes a unimodal alternative (one main peak).
#'  * If your data has two opposite clusters (bimodal or axial data, like a bi-directional line trend),
#'  the Rayleigh test can yield a high/non-significant p-value because the opposing vectors cancel each other out.
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
#' @seealso [mean_resultant_length()], [circular_mean()]
#' @family One-Sample-Tests
#'
#' @name rayleigh-test
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' rayleigh_test(homing, axial = FALSE) # Do not reject null hypothesis.
#' # R = 0.22; stat = 0.497, p = 0.62
#' rayleigh_test_perm(homing, axial = FALSE)
#'
#' # Example data from Davis (1986), pp. 316
#' rayleigh_test(striae, axial = FALSE) # reject null hypothesis
#' rayleigh_test_perm(striae, axial = FALSE) # reject null hypothesis
#'
#' rayleigh_test(striae, mu = 105, axial = FALSE) # reject null hypothesis
#' rayleigh_test_perm(striae, mu = 105, axial = FALSE) # reject null hypothesis
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
#' rayleigh_test_perm(sa.por$azi.PoR, mu = 135, n_perm = 1e3) # reject null hypothesis
NULL

#' @rdname rayleigh-test
#' @export
rayleigh_test <- function(x, mu = NULL, axial = TRUE, alpha = 0.05, quiet = FALSE) {
  f <- if (isTRUE(axial)) 2 else 1
  if (is.null(mu)) {
    x <- x[!is.na(x)]
    xf <- (x * f) %% 360
    n <- length(x)
    R <- mean_resultant_length(xf, na.rm = FALSE)
    S <- 2 * n * R^2
    Z <- S / 2
    p.value <- rayleigh_p_value1(Z, n)

    result <- list(R = R, statistic = Z, p.value = p.value)

    if (isFALSE(quiet)) {
      message(if (p.value < alpha) "Reject Null Hypothesis\n" else "Do Not Reject Null Hypothesis\n")
    }
  } else {
    keep <- !is.na(x)
    x <- x[keep]
    x <- x * f
    mu <- mu * f
    xmu <- x - mu
    n <- length(x)
    C <- (sum(cosd(xmu))) / n
    s <- sqrt(2 * n) * C
    p.value <- rayleigh_p_value2(s, n)
    result <- list(C = C, statistic = s, p.value = p.value)
    if (isFALSE(quiet)) {
      message(if (p.value < alpha) "Reject Null Hypothesis\n" else "Do Not Reject Null Hypothesis\n")
    }
  }
  result
}

rayleigh_test2 <- function(x, axial = TRUE, alpha = 0.05) {
  dname <- deparse1(substitute(x))
  f <- if (isTRUE(axial)) 2 else 1
  x <- x[!is.na(x)]
  xf <- (x * f) %% 360
  n <- length(x)
  R <- mean_resultant_length(xf, na.rm = FALSE)
  p <- 2

  S     <- p * R^2 / n
  S_mod <- (1 - 1 / (2 * n)) * S + S^2 / (2 * n * (p + 2))

  structure(
    list(
      statistic   = c(S = S),
      parameter   = c(n = n, df = p),
      p.value     = stats::pchisq(S_mod, df = p, lower.tail = FALSE),
      estimate    = c(R = R, Rbar = R / n),
      method      = "Rayleigh test of uniformity\n\t(p-value from the modified statistic S*)",
      data.name   = dname,
      alternative = "unimodal distribution",
      modified    = S_mod,
      critical    = stats::qchisq(1-alpha, df = p),
      level       = 1-alpha,
      source      = "Mardia & Jupp (2000), eq. 10.4.4-10.4.6"
    ),
    class = "htest"
  )
}

#' @keywords internal
rayleigh_p_value1 <- function(K, n, wilkie = FALSE) {
  if (isFALSE(wilkie)) {
    P <- exp(-K)
    if (n < 50) {
      temp <- 1 +
        (2 * K - K^2) / (4 * n) -
        (24 * K - 132 * K^2 + 76 * K^3 - 9 * K^4) / (288 * n^2)
    } else {
      temp <- 1
    }
    min(max(P * temp, 0), 1)
  } else {
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


#' @rdname rayleigh-test
#' @export
rayleigh_test_perm <- function(x, mu = NULL, axial = TRUE, n_perm = 1000L){
  dname <- deparse1(substitute(x))
  f <- if (isTRUE(axial)) 2 else 1
  x <- stats::na.omit(x)
  xf <- (x * f) %% 360
  n <- length(x)

  if(is.null(mu)){
  stat <- mean_resultant_length(xf, na.rm = FALSE)
  names(stat) = "R"

  null_dist <- vapply(seq_len(n_perm), function(i) {
    rnd <- runif(n, 0, 360)
    mean_resultant_length(rnd)
  }, FUN.VALUE = numeric(1))

  } else {
    mu <- (mu * f) %% 360
    xmu <- xf - mu

    stat <- (sum(cosd(xmu))) / n
    names(stat) = "C"

    null_dist <- vapply(seq_len(n_perm), function(i) {
      rnd <- runif(n, 0, 360)
      (sum(cosd(rnd))) / n
    }, FUN.VALUE = numeric(1))
  }

  p.value <- (sum(null_dist >= stat) + 1L) / (n_perm + 1)

  structure(
    list(
      statistic   = stat,
      parameter   = c(n = n, df = 2, n_perm = n_perm),
      p.value     = p.value,
      estimate    = c(R = unname(stat) * n, Rbar = unname(stat)),
      method      = paste("Rayleigh test of uniformity\n\t(p-value from", n_perm, "permutations)"),
      data.name   = dname,
      alternative = "unimodal distribution"#,
      #modified    = S_mod,
      #critical    = stats::qchisq(1-alpha, df = 2),
      #level       = 1-alpha
    ),
    class = "htest"
  )
}



#' Weighted Goodness-of-fit Test for Circular Data
#'
#' Weighted version of the Rayleigh test (or V0-test) for uniformity against a
#' distribution with a priori expected von Mises concentration.
#' `weighted_rayleight_test_perm()` uses permutation to estimate p-values.
#'
#' @inheritParams rayleigh_test
#' @param w numeric vector weights of length `length(x)`. If `NULL`, the
#' non-weighted Rayleigh test is performed.
#' @param quiet logical. Prints the test's decision.
#' @param n_perm integer. Number of permutations
#'
#' @details
#' The Null hypothesis is uniformity (randomness). The alternative is a
#' distribution with a (specified) mean direction (`mu`).
#' If `statistic >= p.value`, the null hypothesis of randomness is rejected and
#' angles derive from a distribution with a (or the specified) mean direction.
#'
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
#' @name weighted-rayleigh-test
#' @family One-Sample-Tests
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
#'
#' weighted_rayleigh_perm(tibet.por$azi.PoR, mu = 90, w = 1 / tibet$unc)
#' weighted_rayleigh_perm(ice.por$azi.PoR, mu = 0, w = 1 / iceland$unc)
#' weighted_rayleigh_perm(sa.por$azi.PoR, mu = 135, w = 1 / san_andreas$unc)
NULL

#' @rdname weighted-rayleigh-test
#' @export
weighted_rayleigh <- function(x, mu = NULL, w = NULL, axial = TRUE, quiet = FALSE) {
  if (is.null(w)) {
    rayleigh_test(x, mu = mu, axial = axial)
  } else {
    # remove NA's
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]

    Z <- sum(w)
    n <- length(w)

    if (is.null(mu)) mu <- circular_mean(x, w, axial, na.rm = FALSE)

    d <- x - mu
    f <- if (isTRUE(axial)) 2 else 1

    m <- mean_SC(f * d, w = w, na.rm = FALSE)
    C <- as.numeric(m["C"])
    s <- sqrt(2 * n) * C
    p.value <- rayleigh_p_value2(s, n)

    result <- list(
      C = C,
      statistic = s,
      p.value = p.value
    )
    if (isFALSE(quiet)) {
      if (s >= p.value) {
        message("Reject Null Hypothesis\n")
      } else {
        message("Do Not Reject Null Hypothesis\n")
      }
    }
    return(result)
  }
}

#' @rdname weighted-rayleigh-test
#' @export
weighted_rayleigh_perm <- function(x, mu = NULL, w = NULL, axial = TRUE, n_perm = 1000L){
  if(is.null(w)) {
    rayleigh_test_perm(x, mu, axial, n_perm)
  } else {
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]

  Z <- sum(w)
  n <- length(w)

  if (is.null(mu)) mu <- circular_mean(x, w, axial, na.rm = FALSE)

  d <- x - mu
  f <- if (isTRUE(axial)) 2 else 1

  m <- mean_SC(f * d, w = w, na.rm = FALSE)
  stat <- as.numeric(m["C"])

  null_dist <- vapply(seq_len(n_perm), function(i) {
    rnd <- runif(n, 0, 360)
    m <- mean_SC(rnd)
    as.numeric(m["C"])
  }, FUN.VALUE = numeric(1))

  p.value <- (sum(null_dist >= stat) + 1L) / (n_perm + 1)

  return(list(statistic = stat, p.value = p.value))
  }
}


.uniform_scores <- function(theta) {
  # sorted theta / 2pi in [0, 1)
  sort((as.numeric(theta) %% (2 * pi)) / (2 * pi))
}

## Kuiper -----------------------------------------------------------------

.kuiper_sf <- function(v) {
  # survival function of the asymptotic Kuiper distribution
  if (v < 0.4) return(1)
  j <- 1:100
  p <- 2 * sum((4 * j^2 * v^2 - 1) * exp(-2 * j^2 * v^2))
  min(max(p, 0), 1)
}

.kuiper_isf <- function(alpha) {
  stats::uniroot(function(v) .kuiper_sf(v) - alpha, lower = 0.5, upper = 5,
          tol = .Machine$double.eps^0.75)$root
}

#' Kuiper Test of Circular Uniformity
#'
#' A statistical test used to determine whether a set of angular or circular
#' data points (such as times of day, compass directions, or degrees) are
#' spread out evenly around a circle or if they cluster in some way.
#'
#' @inheritParams rayleigh_test
#' @param quiet logical. Prints the test's decision.
#'
#' @details
#'  The **Null Hypothesis** (\eqn{H_0}{H0}): The data are distributed completely uniformly
#'  (randomly and evenly) around the circle.
#'
#'  The **Alternative Hypothesis** (\eqn{H_1}{H1}): The data are not uniform and
#'  show a preference, clustering, or pattern somewhere on the circle.
#'
#'  The Test Statistic (V or \eqn{D^{+} + D^{-}}{D⁺ + D⁻}): It measures the
#'  greatest positive and negative differences between your data's empirical
#'  cumulative distribution and a theoretical uniform distribution.
#'
#' ## Interpreting the Results
#'
#' * High Test Statistic / Low p-value (\eqn{p < \alpha}{p < alpha}, typically 0.05):
#' You reject the null hypothesis. This means your data are not uniform; they
#' have a significant preferred direction, grouping, or non-random pattern on the circle.
#'
#' * Low Test Statistic / High p-value (\eqn{p \ge 0.05}{p >= 0.05}): You fail to
#' reject the null hypothesis. There is no strong evidence to say the data are
#' different from a flat, uniform distribution. The points appear random across the circle.
#'
# #' If `statistic > p.value`, the null hypothesis is rejected.
# #' If not, randomness (uniform distribution) cannot be excluded.
#'
#' @note Kuiper's test statistic is a rotation-invariant Kolmogorov-type test statistic.
#' The critical values of a modified Kuiper's test statistic are used according
#' to the tabulation given in Stephens (1970).
#'
#' @returns A list with class `"htest"` containing the following components:
#' \describe{
#' \item{`statistic`}{the value of the D statistic}
#' \item{`parameter`}{the number of data}
#' \item{`p.value`}{the p-value for the test}
#' \item{`estimate`}{the D+ and D- values}
#' \item{`alternative`}{a character string describing the alternative hypothesis.}
#' \item{`method`}{a character string indicating what type of Kuiper test was performed.}
#' \item{`data.name`}{a character string giving the name(s) of the data.}
#' \item{`alpha`}{the significance level}
#' \item{`reject`}{logical. whether to reject the null-hypothesis based on the significance level.}
#' }
#'
#' @family One-Sample-Tests
#'
#' @export
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' kuiper_test(homing, alpha = .05)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' kuiper_test(sa.por$azi.PoR, alpha = .05)
kuiper_test <- function(x, alpha = 0.05, axial = TRUE, quiet = FALSE) {
  dname <- deparse1(substitute(x))
  # allowed_alphas <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15)
  # if (!(alpha %in% allowed_alphas)) {
  #   stop("'alpha' must be one of: 0, 0.01, 0.025, 0.05, 0.1, 0.15")
  # }
  #
  # thresholds <- c(1.537, 1.62, 1.747, 1.862, 2.001)
  #
  # kuiper.crits <- cbind(
  #   rev(allowed_alphas[-1]),
  #   thresholds
  # )
  f <- if (isTRUE(axial)) 2 else 1

  x <- x[!is.na(x)] # remove NA's
  x <- (x * f) %% 360
  n <- length(x)
  theta <- deg2rad(x)
  # u <- sort(deg2rad(x) %% (2 * pi)) / (2 * pi)
  # n <- length(x)
  # i <- 1:n
  # D.P <- max(i / n - u)
  # D.M <- max(u - (i - 1) / n)
  # sqrt_n <- sqrt(n)
  # V <- D.P + D.M
  # V <- V * (sqrt_n + 0.155 + 0.24 / sqrt_n)
  #
  # if (alpha == 0) {
  #   labels <- c(
  #     "P-value > 0.15",
  #     "0.10 < P-value < 0.15",
  #     "0.05 < P-value < 0.10",
  #     "0.025 < P-value < 0.05",
  #     "0.01 < P-value < 0.025",
  #     "P-value < 0.01"
  #   )
  #   idx <- findInterval(V, thresholds) + 1
  #   p.value <- labels[idx]
  # } else {
  #   idx <- which(alpha == kuiper.crits[, 1])
  #   p.val.threshold <- kuiper.crits[idx, 2]
  #   p.value <- p.val.threshold
  #
  #   if (isFALSE(quiet)) {
  #     msg <- if (V > p.val.threshold) "Reject Null Hypothesis\n" else "Do Not Reject Null Hypothesis\n"
  #     message(msg)
  #   }
  # }
  # return(
  #   list(
  #     statistic = V,
  #     p.value = unname(p.value)
  #   )
  # )


  u <- .uniform_scores(theta)
  n <- length(u)
  i <- seq_len(n)
  d_plus  <- max(i / n - u)
  d_minus <- max(u - (i - 1) / n)
  V <- (d_plus + d_minus) * (sqrt(n) + 0.155 + 0.24 / sqrt(n))

  structure(
  list(
    method      = "Kuiper's test of uniformity",
    alternative = "non-uniform distribution",
    data.name   = dname,
    statistic = c(V=V),
    p.value   = .kuiper_sf(V),
    critical  = .kuiper_isf(alpha),
    level     = 1 - alpha,
    parameter = c(n = n),
    source    = "Fisher (1993), sec. 4.3, Table A5",
    estimate     = c(`D+` = d_plus, `D-` = d_minus)
  ),
  class = "htest"
  )
}


## Watson U2 --------------------------------------------------------------

.watson_u2_sf <- function(u) {
  if (u <= 0.02) return(1)
  j <- 1:100
  p <- 2 * sum((-1)^(j - 1) * exp(-2 * j^2 * pi^2 * u))
  min(max(p, 0), 1)
}

.watson_u2_isf <- function(alpha, tol = .Machine$double.eps^0.75) {
  stats::uniroot(function(u) .watson_u2_sf(u) - alpha, lower = 0.03, upper = 2,
          tol = tol)$root
}

#' Watson's \eqn{U^2}{U²} Test for Goodness-of-Fit against known Distribution
#'
#' A non-parametric statistical test used for circular data to determine whether
#' a sample fits a specified theoretical distribution
#'
#' @inheritParams rayleigh_test
#' @param alpha Significance level of the test ranging between 0 and 1.
# #' Valid levels are `0.01`, `0.05`, and `0.1`.
# #' This argument may be omitted (`NULL`, the default), in which case, a range
# #' for the p-value will be returned.
#' @param dist Distribution to test for. The default, `"uniform"`, is the
#' circular uniform distribution. `"vonmises"` tests the von Mises distribution.
#'
#' @returns A list with class `"htest"` containing the following components:
#' \describe{
#' \item{`statistic`}{the value of the U2 statistic}
#' \item{`parameter`}{the number of data}
#' \item{`p.value`}{the p-value for the test}
#' \item{`estimate`}{the MLE of the kappa value for the von Mises distribution}
#' \item{`alternative`}{a character string describing the alternative hypothesis.}
#' \item{`method`}{a character string indicating what type of Watson test was performed.}
#' \item{`data.name`}{a character string giving the name(s) of the data.}
#' \item{`alpha`}{the significance level}
#' \item{`reject`}{logical. whether to reject the null-hypothesis based on the significance level.}
#' }
#'
#'
#' @details
#' ## Hypotheses
#' **Null Hypothesis** (\eqn{H_0}{H0}):  The circular sample comes from a specified
#' theoretical distribution (such as a uniform distribution or a specific von Mises distribution).
#'
#' **Alternative Hypothesis** (\eqn{H_1}{H1}): The circular sample does not follow the specified theoretical distribution.
#'
#' ## Interpretation
#' To interpret the output of Watson's \eqn{U^2}{U²}  test, compare your
#' calculated \eqn{U^2}{U²}  test statistic to the critical value from
#' Watson's goodness-of-fit/homogeneity tables at your chosen significance level (\eqn{\alpha}{alpha}, commonly set to 0.05),
#' or check the resulting p-value:
#'
#'  * If \eqn{U^2_{\text{calculated}} > U^2_{\text{critical}}} (or p < \eqn{\alpha}{alpha}):
#'  Reject the null hypothesis (\eqn{H_0}{H0}).
#'  Conclude that the data significantly deviates from the theoretical distribution.
#'
#'  * If \eqn{U^2_{\text{calculated}} \le U^2_{\text{critical}}} (or \eqn{p \ge \alpha}{p >= alpha}):
#'  Fail to reject the null hypothesis (\eqn{H_0}{H0}).
#'  There is not enough evidence to claim the data deviates from the expected model.
#'
#' @note Watson's test statistic is a rotation-invariant Cramer - von Mises test.
#' non-parametric, rank-based alternative to one-sample
#'
#' @references Mardia and Jupp (2000). Directional Statistics. John Wiley and
#' Sons.
#'
#' @export
#'
#' @family One-Sample-Tests
#' @seealso [vonmises] and [cunif]
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' watson_test(homing, axial = FALSE, alpha = .05)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' watson_test(sa.por$azi.PoR, alpha = .05)
#' watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
#' watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
watson_test <- function(x, alpha = 0.05, dist = c("uniform", "vonmises"), axial = TRUE, quiet = FALSE) {
  dname <- deparse1(substitute(x))
  # if (is.null(alpha)) alpha <- 0
  # allowed_alphas <- c(0, 0.01, 0.025, 0.05, 0.1)
  # if (!(alpha %in% allowed_alphas)) {
  #   stop("'alpha' must be one of: 0, 0.01, 0.025, 0.05, 0.1")
  # }
  dist <- match.arg(dist)
  # the vonmises critical-value table only has columns for alpha = 0.1/0.05/0.01 --
  # reject early instead of failing later after doing all the work
  if (dist == "vonmises" && alpha != 0 && !(alpha %in% c(0.1, 0.05, 0.01))) {
    stop("'alpha' must be one of: 0, 0.1, 0.05, 0.01 when dist = 'vonmises'")
  }

  x <- x[!is.na(x)]
  n <- length(x)
  f <- if (isTRUE(axial)) 2 else 1
  x2 <- (x * f) %% 360   # doubled-domain data used throughout, for both dist branches

  if (n < 8) {
    warning("Total Sample Size < 8:  Results are not valid")
  }

  if (dist == "uniform") {
    ## ---- U2 statistic ----
    # u <- sort(deg2rad(x2)) / (2 * pi)
    # u.bar <- mean(u)
    # i <- seq_len(n)
    # u2 <- sum((u - u.bar - (i - .5) / n + .5)^2) + 1 / (12 * n)
    # statistic <- (u2 - 0.1 / n + 0.1 / (n^2)) * (1 + 0.8 / n)
    #
    # ## crits correspond to alpha = 0.01, 0.025, 0.05, 0.1 respectively
    # crits       <- c(0.267, 0.221, 0.187, 0.152)
    # crit_alphas <- c(0.01, 0.025, 0.05, 0.1)
    #
    # if (n < 8) {
    #   p.value <- NA_character_
    #   critical.value <- NA_real_
    #   reject <- NA
    # } else if (alpha == 0) {
    #   thresholds <- c(rev(crits), Inf)                 # c(0.152, 0.187, 0.221, 0.267, Inf)
    #   messages <- c(
    #     "P-value > 0.10",
    #     "0.05 < P-value < 0.10",
    #     "0.025 < P-value < 0.05",
    #     "0.01 < P-value < 0.025",   # fixed: was "0.001" -- table only resolves to alpha=0.01
    #     "P-value < 0.01"            # fixed: was "0.001"
    #   )
    #   p.value <- messages[findInterval(statistic, thresholds) + 1]   # fixed: was missing "+ 1"
    #   critical.value <- NA_real_
    #   reject <- NA
    # } else {
    #   critical.value <- crits[crit_alphas == alpha]
    #   reject <- statistic > critical.value
    #   p.value <- NA_character_
    #   if (!quiet) {
    #     message(if (reject) "Reject Null Hypothesis" else "Do Not Reject Null Hypothesis")
    #   }
    # }

    theta <- deg2rad(x)
    u <- .uniform_scores(theta)
    n <- length(u)
    i <- seq_len(n)
    u2 <- sum((u - (2 * i - 1) / (2 * n))^2) - n * (mean(u) - 0.5)^2 +
      1 / (12 * n)
    u_mod <- (u2 - 0.1 / n + 0.1 / n^2) * (1 + 0.8 / n)

    name      = "Watson's U2 test against uniformity"
    alternative = "non-uniform distribution"
    statistic <- c("U2" = u_mod)
    p.value <- .watson_u2_sf(u_mod)
    critical.value <- .watson_u2_isf(alpha)
    reject <- statistic > critical.value
    extra     = NULL


  } else {
    ## ---- goodness-of-fit against von Mises ----

    if (is.null(alpha)) alpha <- 0
    allowed_alphas <- c(0, 0.01, 0.025, 0.05, 0.1)
    if (!(alpha %in% allowed_alphas)) {
      stop("'alpha' must be one of: 0, 0.01, 0.025, 0.05, 0.1")
    }

    u2_crits <- cbind(
      c(0, 0.5, 1, 1.5, 2, 4, 100),
      c(0.052, 0.056, 0.066, 0.077, 0.084, 0.093, 0.096),
      c(0.061, 0.066, 0.079, 0.092, 0.101, 0.113, 0.117),
      c(0.081, 0.09, 0.11, 0.128, 0.142, 0.158, 0.164)
    )
    kappa_mle <- est.kappa(x2)
    ## fixed: mean must be taken in the SAME (doubled) domain as x2/kappa_mle.
    ## the original computed circular_mean(x, axial=axial, ...) -- the mean of
    ## the UNdoubled data, in undoubled units -- then subtracted it from x2.
    mu <- circular_mean(x2, axial = FALSE, na.rm = FALSE)
    z <- pvm(x2 - mu, mean = 0, kappa = kappa_mle)   # vectorized directly; no apply()/matrix needed
    z <- sort(z)
    z.bar <- mean(z)
    i <- seq_len(n)
    statistic <- sum((z - (2 * i - 1) / (2 * n))^2) - n * (z.bar - 0.5)^2 + 1 / (12 * n)
    row <- findInterval(kappa_mle, c(0.25, 0.75, 1.25, 1.75, 3, 5)) + 1

    if (n < 8) {
      p.value <- NA_character_
      critical.value <- NA_real_
      reject <- NA
    } else if (alpha == 0) {
      breaks <- u2_crits[row, 2:4]
      labels <- c(
        "P-value > 0.10",
        "0.05 < P-value < 0.10",
        "0.01 < P-value < 0.05",
        "P-value < 0.01"
      )
      p.value <- labels[findInterval(statistic, breaks, rightmost.closed = TRUE) + 1]
      critical.value <- NA_real_
      reject <- NA
    } else {
      col <- match(alpha, c(0.1, 0.05, 0.01)) + 1
      critical.value <- u2_crits[row, col]
      reject <- statistic > critical.value
      p.value <- NA_character_
      # if (!quiet) {
      #   message(if (reject) "Reject Null Hypothesis" else "Do Not Reject Null Hypothesis")
      # }
    }
    name      = "Watson's U2 test against von Mises distribution"
    alternative = "Sample does not follow a von Mises distribution"
    extra <- c("Kappa MLE" = kappa_mle)
    names(statistic) = "U2"
  }

  structure(
    list(
      method      = name,
      alternative = alternative,
      statistic      = statistic,
      p_value        = p.value,        # descriptive bracket, only set when alpha == 0
      critical = critical.value, # tabulated threshold, only set when alpha != 0
      reject         = reject,         # logical decision, only set when alpha != 0
      level = 1- alpha,
      alpha          = alpha,
      parameter = c(n = n),
      #dist           = dist,
      estimate = extra,
      data.name = dname
    ),
    class = "htest"
  )
}

## Hodges-Ajne ------------------------------------------------------------

#' Hodges-Ajne Test
#'
#' Hodges-Ajne test of uniformity on the circle
#'
#' @inheritParams rayleigh_test
#'
#' @details  The statistic \eqn{m} is the smallest number of angles in a half-circle. Its
#' exact distribution is \deqn{P(m \le t) = 2^{(1 - n)} (n - 2t) C(n, t)} for \eqn{t < n / 3}
#' (Mardia and Jupp, 2000, sec. 6.3.4), which gives the p-value (1 above \eqn{n / 3}).
#' Small \eqn{m} indicates a preferred direction.
#'
#' @returns A list with class `"htest"` containing the following components:
#' \describe{
#' \item{`statistic`}{the value of the m statistic}
#' \item{`parameter`}{the number of data}
#' \item{`p.value`}{the p-value for the test}
#' \item{`alternative`}{a character string describing the alternative hypothesis.}
#' \item{`method`}{a character string indicating what type of test was performed.}
#' \item{`data.name`}{a character string giving the name(s) of the data.}
#' }
#'
#' @family One-Sample-Tests
#' @export
#'
#' @examples
#' hodges_ajne_test(rvm(100, mean = 0, kappa = 10))
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' hodges_ajne_test(sa.por$azi.PoR, alpha = .05)
hodges_ajne_test <- function(x, alpha = 0.05, axial = TRUE) {
  f <- if (isTRUE(axial)) 2L else 1L
  x0 <- x[!is.na(x)]
  n <- length(x0)
  f <- if (isTRUE(axial)) 2 else 1
  x2 <- x0 * f   # doubled-domain data used throughout, for both dist branches
  theta <- deg2rad(x2)

  t <- sort(as.numeric(theta) %% (2 * pi))
  n <- length(t)
  # number of angles in the half-circles starting at each angle
  ends   <- findInterval(t + pi, c(t, t + 2 * pi), left.open = TRUE)
  counts <- ends - (seq_len(n) - 1L)
  m <- min(counts, n - counts)
  p <- 1
  if (m < n / 3) {
    p <- min(1, exp((1 - n) * log(2) + log(n - 2 * m) + lchoose(n, m)))
  }
  structure(
  list(
    method      = "Hodges-Ajne test of uniformity",
    alternative = "unimodal distribution",
    statistic = c(m = m),
    p.value   = p,
    level     = 1 - alpha,
    parameter   = c(n = n),
    data.name   = deparse1(substitute(x)),
    source    = "Mardia & Jupp (2000), sec. 6.3.4"
  ),
  class = "htest"
  )
}

## Rao spacing ------------------------------------------------------------

#' Rao Spacing Test
#'
#' Rao's spacing test of uniformity on the circle
#'
#' @inheritParams rayleigh_test
#'
#' @details Test statistic is
#' \deqn{L = 1/2 \sum(|T_i - 360 / n|)} in degrees, with the spacings \eqn{T_i}
#' between consecutive angles (Mardia & Jupp 2000, eq. 6.3.45). The p-value is the proportion
#' of ``n_perm`` uniform samples with a larger \eqn{L}
#'
#' @returns A list with class `"htest"` containing the following components:
#' \describe{
#' \item{`statistic`}{the value of the L statistic}
#' \item{`parameter`}{the number of data and the number of permutations}
#' \item{`p.value`}{the p-value for the test}
#' \item{`alternative`}{a character string describing the alternative hypothesis.}
#' \item{`method`}{a character string indicating what type of test was performed.}
#' \item{`data.name`}{a character string giving the name(s) of the data.}
#' }
#'
#' @family One-Sample-Tests
#' @export
#'
#' @examples
#' set.seed(20250411)
#' rao_spacing_test(rvm(100, mean = 0, kappa = 10))
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' rao_spacing_test(sa.por$azi.PoR, alpha = .05)
rao_spacing_test <- function(x, alpha = 0.05, n_perm = 9999L, axial = TRUE) {
  f <- if (isTRUE(axial)) 2L else 1L
  x0 <- x[!is.na(x)]
  n <- length(x0)
  f <- if (isTRUE(axial)) 2 else 1
  x2 <- x0 * f   # doubled-domain data used throughout, for both dist branches
  theta <- deg2rad(x2)

  # spacing statistic L (degrees) for each row of an n_perm x n matrix
  spacing <- function(M) {
    M <- t(apply(M %% (2 * pi), 1L, sort))
    k <- ncol(M)
    gaps <- cbind(M[, -1L, drop = FALSE] - M[, -k, drop = FALSE],
                  M[, 1L] + 2 * pi - M[, k])
    0.5 * rowSums(abs(gaps * 180 / pi - 360 / k))
  }

  theta <- as.numeric(theta)
  n <- length(theta)
  L <- spacing(matrix(theta, nrow = 1L))
  sims <- spacing(matrix(runif(n_perm * n, 0, 2 * pi), n_perm, n))

  level <- 1 - alpha

  structure(
  list(
    method      = paste0("Rao's spacing test of uniformity with simulated ",
                         "p-value\n\t (based on ", n_perm, " replicates)"),
    statistic = c("L" = L),
    p.value   = (sum(sims >= L) + 1) / (n_perm + 1),
    critical  = stats::quantile(sims, probs = level, names = FALSE),
    level     = level,
    data.name   = deparse1(substitute(x)),
    alternative = "non-uniform distribution",
    parameter   = c(n = n, n_perm = n_perm),
    source    = "Mardia & Jupp (2000), eq. 6.3.45"
  ),
  class = 'htest'
  )
}



# Comparing samples ----------------------------------------------------------------

#' Watson's Two-Sample Test of Homogeneity
#'
#' Performs Watson's test for homogeneity on two samples of circular data.
#' `watson_two_test_perm()` uses permutation to estimate p-values.
#' Test whether samples of directions share the mean direction.
#'
#' @param x,y numeric vectors. Angles in degrees
#' @param n_perm integer. Number of permutations
#' @inheritParams watson_test
#'
#' @details A two-sample Watson's \eqn{U^2}{U²} permutation test determines whether two
#' independent groups of circular data (angles or directions) come from the same
#' underlying distribution.
#'
#' ## Hypotheses
#'
#' **Null Hypothesis** (\eqn{H_0}{H0}): The two samples come from the same circular
#' distribution (the two groups of angles share a common distribution around the circle).
#'
#' **Alternative Hypothesis** (\eqn{H_1}{H1}): The two samples come from different
#' circular distributions (the groups are oriented or dispersed differently around the circle).
#'
#' ## Interpretation
#'
#' The Test Statistic (\eqn{U^2}{U²}) measures the distance between the cumulative distribution
#' functions of the two circular samples. A larger \eqn{U^2}{U²} value means the
#' two sets of angles look more different from each other.
#'
#' The P-Value represents the probability of getting a U² value as large as (or
#' larger than) your observed value purely by chance, assuming the null hypothesis
#' is true. It is calculated by shuffling the group labels across your data many
#' times to build a reference permutation distribution.
#'
#' ### Making a Decision
#'
#' *  Low p-value (\eqn{p \le \alpha}{p <= alpha}, usually 0.05): Reject the null hypothesis. Conclude that
#' the two groups have significantly different circular distributions.
#' * High p-value (\eqn{p > \alpha}{p > alpha}): Fail to reject the null hypothesis. There is not
#' enough evidence to say the two groups are distributed differently around the circle.
#'
#' @note Critical values for the test statistic are obtained using the asymptotic
#' distribution of the test statistic. It is recommended to use the obtained
#' critical values and ranges for p-values only for combined sample sizes in
#' excess of 17. Tables are available for smaller sample sizes and can be found
#' in Mardia (1972) for instance.
#'
#' @name watson_two_sample
#' @family Two-Sample-Tests
#'
#' @returns A list with class `"htest"`
#'
#' @examples
#' set.seed(20250411)
#' x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
#' x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#' watson_two_test(x1, x2, axial = FALSE)
#' watson_two_test_perm(x1, x2, axial = FALSE)
#'
#' data1 <- rvm(n=20, mean = 0, kappa=3)
#' data2 <- rvm(n=20, mean = 90, kappa=2)
#' watson_two_test(data1, data2, axial = FALSE)
#' watson_two_test_perm(data1, data2, axial = FALSE)
#'
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' watson_two_test(sa.por$azi, 135, alpha = 0.05)
#' watson_two_test(sa.por$azi, rvm(100, 135, 10), alpha = 0.05)
#' watson_two_test_perm(sa.por$azi, rvm(100, 135, 10), alpha = 0.05)
NULL

#' @export
#' @rdname watson_two_sample
watson_two_test <- function(x, y, alpha = 0.05, axial = TRUE, quiet = FALSE) {
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  f <- if (isTRUE(axial)) 2L else 1L

  # # --- Preprocessing ---
  x <- stats::na.omit(x)
  y <- stats::na.omit(y)

  # if (is.null(alpha)) alpha <- 0
  # allowed_alphas <- c(0, 0.001, 0.01, 0.05, 0.1)
  # if (!(alpha %in% allowed_alphas)) {
  #   stop("'alpha' must be one of: 0, 0.001, 0.01, 0.05, 0.1")
  # }
  #
  #
  # f <- if (axial) 2L else 1L
  # x <- deg2rad((f * x) %% 360)
  # y <- deg2rad((f * y) %% 360)
  #
  # # --- U2 statistic ---
  # n1 <- length(x)
  # n2 <- length(y)
  # n <- n1 + n2
  #
  # labels <- c(rep(1L, n1), rep(2L, n2))
  # labels <- labels[order(c(x, y))]
  #
  # cx <- cumsum(labels == 1L)
  # cy <- cumsum(labels == 2L)
  #
  # d <- cy / n2 - cx / n1
  # dbar <- mean.default(d)
  # u2 <- (n1 * n2) / n^2 * sum((d - dbar)^2)
  #
  # if (n < 18) warning("Total sample size < 18: consult tabulated critical values")
  #
  # # --- Critical values (parallel to allowed_alphas) ---
  # crits <- c(NA_real_, 0.385, 0.268, 0.187, 0.152)
  #
  # if (alpha == 0) {
  #   # Return a descriptive p-value bracket
  #   thresholds <- c(0.152, 0.187, 0.268, 0.385)
  #   labels_p <- c(
  #     "P-value > 0.10",
  #     "0.05 < P-value < 0.10",
  #     "0.01 < P-value < 0.05",
  #     "0.001 < P-value < 0.01",
  #     "P-value < 0.001"
  #   )
  #   p.value <- labels_p[findInterval(u2, thresholds, rightmost.closed = TRUE) + 1L]
  # } else {
  #   crit <- crits[match(alpha, allowed_alphas)]
  #   if (!quiet) {
  #     message(if (u2 > crit) "Reject null hypothesis" else "Do not reject null hypothesis")
  #   }
  #   p.value <- crit
  # }
  #
  # list(statistic = u2, p.value = p.value)



  #level <-  1 - alpha
  theta <- list(x, y)
  theta <- lapply(theta, function(x){
    x2 <- x * f
    deg2rad(x2)
  })
  #if (length(theta) != 2L) stop("Watson's two-sample U2 compares two samples")
  n1 <- length(theta[[1L]])
  n2 <- length(theta[[2L]])
  n  <- n1 + n2
  pooled <- unlist(theta, use.names = FALSE) %% (2 * pi)
  first <- rep(c(1, 0), c(n1, n2))[order(pooled)]     # order() is stable
  a <- cumsum(first)
  d <- (seq_len(n) - a) / n2 - a / n1
  u2 <- n1 * n2 / n^2 * sum((d - mean(d))^2)

  structure(
    list(
      method      = "Watson two-sample U2 test",
      alternative = "Samples come from different distributions",
      parameter = c("n" = n),
      statistic = c("U2" = u2),
      p.value   = .watson_u2_sf(u2),
      critical  = .watson_u2_isf(alpha),
      level     = 1 - alpha,
      data.name = dname,
      source    = "Mardia & Jupp (2000), eq. 8.3.7"
    ),
    class = "htest"
  )
}


# watson_two_test2 <- function(x, mu, kappa, alpha = 0, axial = TRUE, quiet = FALSE){
#   if(is.null(alpha)) alpha <- 0
#   allowed_alphas <- c(0, 0.001, 0.01, 0.05, 0.1)
#   if (!(alpha %in% allowed_alphas)) {
#     stop("'alpha' must be one of: 0, 0.001, 0.01, 0.05, 0.1")
#   }
#
#   # --- Preprocessing ---
#   x <- x[!is.na(x)]
#
#   f <- if (axial) 2L else 1L
#   x <- deg2rad((f * x) %% 360)
#   mu <- (f * mu) %% 360
#
#   n <- length(x)
#   if (n < 18) warning("Total sample size < 18: consult tabulated critical values")
#
#   # --- Evaluate von Mises CDF at order statistics ---
#   x_sorted <- sort(x)
#   v <- pvm(rad2deg(x_sorted), mean = mu, kappa = kappa)
#
#   # --- U2 statistic ---
#   i    <- seq_len(n)
#   vbar <- mean.default(v)
#   u2   <- sum((v - vbar + (2*i - 1)/(2*n) - 0.5)^2) + 1/(12*n)
#
#   # --- Critical values (parallel to allowed_alphas) ---
#   crits <- c(NA_real_, 0.385, 0.268, 0.187, 0.152)
#
#   if (alpha == 0) {
#     # Return a descriptive p-value bracket
#     thresholds <- c(0.152, 0.187, 0.268, 0.385)
#     labels_p   <- c(
#       "P-value > 0.10",
#       "0.05 < P-value < 0.10",
#       "0.01 < P-value < 0.05",
#       "0.001 < P-value < 0.01",
#       "P-value < 0.001"
#     )
#     p.value <- labels_p[findInterval(u2, thresholds, rightmost.closed = TRUE) + 1L]
#   } else {
#     crit <- crits[match(alpha, allowed_alphas)]
#     if (!quiet) {
#       message(if (u2 > crit) "Reject null hypothesis" else "Do not reject null hypothesis")
#     }
#     p.value <- crit
#   }
#
#   list(statistic = u2, p.value = p.value)
# }
#

#' @export
#' @rdname watson_two_sample
watson_two_test_perm <- function(x, y, axial = TRUE, n_perm = 1000L, alpha = NULL) {
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  f <- if (isTRUE(axial)) 2L else 1L

  samp1 <- deg2rad((f * stats::na.omit(x) %% 360))
  samp2 <- deg2rad((f * stats::na.omit(y) %% 360))

  n1 <- length(samp1)
  n2 <- length(samp2)
  n <- n1 + n2
  combined <- c(samp1, samp2)

  # Bare U² computation — no validation, no output, just the number
  u2_stat <- function(s1, s2, n1, n2) {
    labels <- c(rep(1L, n1), rep(2L, n2))
    labels <- labels[order(c(s1, s2))]
    cx <- cumsum(labels == 1L)
    cy <- cumsum(labels == 2L)
    d <- cy / n2 - cx / n1
    dbar <- mean.default(d)
    (n1 * n2) / (n1 + n2)^2 * sum((d - dbar)^2)
  }

  Gstat <- u2_stat(samp1, samp2, n1, n2)

  nxtrm <- sum(vapply(seq_len(n_perm), function(r) {
    perm <- sample.int(n)
    Grand <- u2_stat(
      combined[perm[seq_len(n1)]],
      combined[perm[seq(n1 + 1L, n)]],
      n1, n2
    )
    Grand >= Gstat
  }, FUN.VALUE = logical(1))) + 1L

  pvalue <- nxtrm / (n_perm + 1)
  reject <- NULL
  if(!is.null(alpha)) {
    stopifnot(alpha >= 0  & alpha <= 1)
    reject <- pvalue < alpha
  }

structure(
  list(
    method      = paste("Watson two-sample U2 test based on", n_perm, "permutations"),
    alternative = "Samples come from different distributions",
    statistic = c("G" = Gstat),
    p.value = pvalue,
    parameter = c("n" = n, 'n_perm' = n_perm),
    alpha = alpha,
    data.name = dname,
    reject = reject),
  class = "htest"
)
}


#' Watson-Wheeler Test of Homogeneity of Means
#'
#' A a non-parametric statistical test used to determine whether two or more
#' independent samples of circular data (angles, directions, or periodic times)
#' come from the same underlying population distribution. The difference
#' between the samples can be in either the mean or the variance.
#'
#' @inheritParams watson_two_test
#' @importFrom stats na.omit
#'
#' @details
#' ## Hypotheses
#'
#' **Null Hypothesis (\eqn{H_0}{H0})** The samples come from identical populations (meaning
#' both the mean direction and the dispersion/variance are homogeneous across groups).
#'
#' **Alternative Hypothesis (\eqn{H_1}{H1})**: At least one sample comes from a different
#' population distribution, which can be due to a difference in the mean direction,
#' a difference in variance/concentration, or both.
#'
#' ## Interpretation
#'
#' **Test Statistic (W)**: This value follows an approximate \eqn{\chi^2}{chi-squared} distribution.
#' Higher values of W indicate larger discrepancies between the angular distributions of your groups.
#'
#' * If the p-value is less than your significance level (commonly \eqn{\alpha}{alpha} = 0.05),
#' you **reject** the null hypothesis.
#' This means you have strong evidence that the groups differ significantly in
#' their central direction or spread around the circle.
#' * If the p-value is greater than 0.05, you **fail to reject** the null hypothesis,
#' meaning there is no statistically significant evidence of difference among the groups.
#'
#' @note Important Considerations & Limitations:
#' * **Sensitivity to both mean and variance**: Because it detects differences in
#' either mean or variance, a significant result does not automatically mean
#' the mean angles are different; it could be driven entirely by differences
#' in concentration (variance).
#' * **Sample size requirement**: The chi-squared approximation requires each group
#' to have a minimum sample size (typically at least 10 elements per group) to remain valid.
#'
#' @family Two-Sample-Tests
#'
#' @returns A list with class `"htest"`
#' @export
#'
#' @examples
#' set.seed(20250411)
#' x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
#' x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#' watson_wheeler_test_perm(x1, x2, axial = FALSE)
#'
#' data1 <- rvm(n=20, mean = 0, kappa=3)
#' data2 <- rvm(n=20, mean = 90, kappa=2)
#' watson_wheeler_test_perm(data1, data2, axial = FALSE)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' watson_wheeler_test_perm(sa.por$azi.PoR, rvm(100, 135, 10))
watson_wheeler_test_perm <- function(x, y, axial = TRUE, n_perm = 1000L, alpha = NULL) {
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  f <- if (isTRUE(axial)) 2L else 1L

  samp1 <- deg2rad((f * stats::na.omit(x) %% 360))
  samp2 <- deg2rad((f * stats::na.omit(y) %% 360))

  n1 <- length(samp1)
  n2 <- length(samp2)
  n <- n1 + n2
  combined <- c(samp1, samp2)

  # Bare Watson-Wheeler statistic
  # W = (2/n) * (C1^2 + S1^2) where C1, S1 are summed uniform scores for group 1
  ww_stat <- function(s1, s2, n1, n2) {
    n <- n1 + n2
    rnks <- rank(c(s1, s2), ties.method = "random")
    cr <- rnks * 2 * pi / n
    C <- sum(cos(cr[seq_len(n1)]))
    S <- sum(sin(cr[seq_len(n1)]))
    2 * (n - 1) * (C^2 + S^2) / (n1 * n2)
  }

  Gstat <- ww_stat(samp1, samp2, n1, n2)

  nxtrm <- sum(vapply(seq_len(n_perm), function(r) {
    perm <- sample.int(n)
    Grand <- ww_stat(
      combined[perm[seq_len(n1)]],
      combined[perm[seq(n1 + 1L, n)]],
      n1, n2
    )
    Grand >= Gstat
  }, FUN.VALUE = logical(1))) + 1L

  pvalue <- nxtrm / (n_perm + 1)
  reject <- NULL
  if(!is.null(alpha)) {
    stopifnot(alpha >= 0  & alpha <= 1)
    reject <- pvalue < alpha
  }

  structure(
    list(
      method      = paste("Watson-Wheeler Test of Homogeneity of Means based on", n_perm, "permutations"),
      alternative = "Samples come from different distributions",
      statistic = c("G" = Gstat),
      p.value = pvalue,
      parameter = c("n" = n, 'n_perm' = n_perm),
      alpha = alpha,
      data.name = dname,
      reject = reject),
    class = "htest"
  )
}

.vonmises_kappa_mle <- function(Rbar, tol = .Machine$double.eps^0.75) {
  # solve A1(kappa) = I1(kappa) / I0(kappa) = Rbar
  if (Rbar < 1e-12) return(0)
  if (Rbar >= 1) return(Inf)
  A1 <- function(k) besselI(k, 1, expon.scaled = TRUE) /
    besselI(k, 0, expon.scaled = TRUE)
  stats::uniroot(function(k) A1(k) - Rbar, lower = 0, upper = 1 / (1 - Rbar),
          extendInt = "upX", tol = tol)$root
}

#' Watson-Williams Test of Homogeneity of Means
#'
#' Performs the Watson-Williams test for homogeneity of means between several
#' samples of circular data.
#'
#' @inheritParams uniform_scores_test
#'
#' @details
#' The Watson-Williams test for the homogeneity of means between two or more
#' groups is performed and the results are printed. The null hypothesis is that
#' means are equal across groups.
#'
#' The assumptions are that: (1) the samples are drawn from populations with a
#' von Mises distribution; (2) the parameter of concentration has the same
#' value in all populations; (3) this parameter is sufficiently large (i.e. > 1).
#' Assumptions 2 and 3 are checked and a warning is issued if they are not met.
#'
#' @family Two-Sample-Tests
#'
#' @returns A list with class `"htest"`
#' @export
#'
#' @examples
#' x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
#' x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#' watson_williams_test(list(x1, x2), axial = FALSE)
#'
#' data1 <- rvm(n=20, mean = 0, kappa=3)
#' data2 <- rvm(n=20, mean = 90, kappa=2)
#' watson_williams_test(list(data1, data2), axial = FALSE)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' watson_williams_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
watson_williams_test <- function(samples, alpha = 0.05, axial = TRUE) {
  level <- 1 - alpha
  if (!is.list(samples) || length(samples) < 2L) {
    stop("`samples` must be a list of at least two samples")
  }
  k  <- length(samples)
  nm <- names(samples)
  if (is.null(nm) || any(!nzchar(nm))) nm <- as.character(seq_len(k))
  nm <- make.unique(nm)

  f <- if (isTRUE(axial)) 2 else 1
  theta <- lapply(samples, function(x) {
    x <- as.numeric(x)
    f * deg2rad(x[is.finite(x)])                       # doubled if axial
  })
  n_i <- lengths(theta)
  if (any(n_i < 2L)) {
    stop("Each sample needs at least two observations (after removing missing values)")
  }

  C <- vapply(theta, function(t) sum(cos(t)), numeric(1L))
  S <- vapply(theta, function(t) sum(sin(t)), numeric(1L))
  n    <- sum(n_i)
  sumR <- sum(sqrt(C^2 + S^2))
  R    <- sqrt(sum(C)^2 + sum(S)^2)

  # pooled concentration and applicability checks
  kappa <- .vonmises_kappa_mle(sumR / n)
  if (kappa <= 0) stop("Pooled concentration is zero; the test is not applicable")
  thresh <- if (k == 2L) 2 else 1
  if (kappa < thresh) {
    warning(sprintf("Pooled concentration %.3g < %g; the test is probably not applicable",
                    kappa, thresh), call. = FALSE)
  }
  kt <- equal_kappa_test(samples, alpha = alpha,, axial = axial)
  kap_i <- if (!is.null(kt$kappa)) kt$kappa else kt$estimate
  if (kt$p.value < 1 - level) {
    warning(sprintf("Concentrations (%s) differ between samples (p = %.3g); the test might not be applicable",
                    paste(format(kap_i, digits = 3), collapse = ", "), kt$p.value),
            call. = FALSE)
  }

  F  <- (1 + 3 / (8 * kappa)) * (n - k) * (sumR - R) / ((k - 1) * (n - sumR))
  df <- c(df1 = k - 1, df2 = n - k)


  mu <- lapply(samples, function(x){
    circular_mean(x, axial = axial)
  }) |> unlist()

  structure(
    list(
      statistic  = c(F = F),
      parameter  = df,
      p.value    = stats::pf(F, df[[1L]], df[[2L]], lower.tail = FALSE),
      estimate   = setNames(mu, paste("mean direction", nm)),
      method     = paste("Watson-Williams test for homogeneity of",
                         if (f == 2) "mean axes" else "mean directions"),
      data.name  = paste(nm, collapse = " and "),
      kappa      = kappa,
      kappa_test = kt,
      level      = level
    ),
    class = "htest"
  )
}



#' Uniform scores test
#'
#' q-sample uniform-scores test tests the equality of given distributions, ie.
#' whether samples of directions share the mean direction.
#' The test is the distribution-free test of Mardia, Watson & Wheeler
#' on the uniform scores \eqn{2 \pi r_i / n} of the ranks in the pooled sample.
#'
#' @details \deqn{W = 2 \sum_{i=1}^q R_i^2 / n_i} with \eqn{2(q - 1)} degrees of freedom,
#' \eqn{R^* = 2 (n - 1) R_1^2 / (n_1 n_2)} for two samples
#'
#' @param samples list of samples of orientations
#' @inheritParams watson_two_test
#'
#' @returns A list with class `"htest"`
#' @family Two-Sample-Tests
#' @export
#'
#' @examples
#' set.seed(20250411)
#' x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
#' x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#' uniform_scores_test(list(x1, x2))
#'
#' #' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' uniform_scores_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
uniform_scores_test <- function(samples, alpha = 0.05, axial = TRUE) {
  level <-  1 - alpha
  dname <- paste(names(x), collapse = " and ")

  f <- if (isTRUE(axial)) 2L else 1L
  theta <- lapply(samples, function(x){
    deg2rad(x * f)
  })

  pooled <- unlist(theta, use.names = FALSE) %% (2 * pi)
  n <- length(pooled)
  groups <- rep(seq_along(theta), lengths(theta))
  # ties within a sample get consecutive ranks, ties between samples the
  # average rank
  ranks <- rank(pooled, ties.method = "first")
  for (value in unique(pooled[duplicated(pooled)])) {
    tied <- pooled == value
    if (length(unique(groups[tied])) > 1L) ranks[tied] <- mean(ranks[tied])
  }
  beta <- 2 * pi * ranks / n
  CS  <- rowsum(cbind(cos(beta), sin(beta)), groups)   # q x 2: C_i, S_i
  n_i <- tabulate(groups)
  W <- if (length(theta) == 2L) {
    2 * (n - 1) * sum(CS[1L, ]^2) / (n_i[1L] * n_i[2L])
  } else {
    2 * sum(rowSums(CS^2) / n_i)
  }

  df <- 2 * (length(theta) - 1)
  structure(
  list(
    method      = "Uniform scores test",
    alternative = "Samples come from different distributions",
    statistic = c(W = W),
    parameter = c(n = n),
    p.value   = stats::pchisq(W, df, lower.tail = FALSE),
    df        = df,
    critical  = stats::qchisq(level, df),
    level     = level,
    data.name = dname,
    source    = "Mardia & Jupp (2000), eq. 8.3.3-8.3.4, 8.4.3"
  ),
  class = "htest"
  )
}

#' Equal \eqn{\kappa} Test
#'
#' This function tests for the homogeneity of concentration parameters for
#' multiple samples of directional data.
#'
#' @inheritParams uniform_scores_test
#' @param method character. One of
#' \describe{
#' \item{`"f"`}{(default for two samples): \eqn{F = ((n_1 - R_1) / \nu_1) / ((n_2 - R_2) / \nu_2)},
#' two-sided, \eqn{F(\nu_1, \nu_2)} (Mardia and Jupp, 2000 eq. 10.5.3; Watson, 1956).}
#' \item{`"bartlett"`}{(default for more samples): Bartlett's test of the
#'       dispersions \eqn{(n_i - R_i) / \nu_i}, chi-square with \eqn{q - 1} degrees of
#'       freedom (Mardia and Jupp, 2000 eq. 10.6.21; McFadden & Lowes, 1981, eq. 64-66).}
#' }
#'
#' @details
#' The samples are assumed to have been drawn from von Mises populations.
#' The null hypothesis tested is that all populations sampled have the same concentration parameter, kappa \eqn{\kappa}.
#'
#' For concentrated data \eqn{2 \kappa (n_i - R_i)} is approximately chi-square
#' with \eqn{\nu_i = (p - 1)(n_i - 1)} degrees of freedom in 2 dimensions.
#'
#' @note Both need concentrated data (\eqn{\bar{R} \ge 0.7};
#' MJ eq. 10.6.24-10.6.27 and 7.3.24-7.3.25 give variance-stabilised tests for dispersed data.
#'
#' @references
#' Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 5.3, World Scientific Press, Singapore.
#'
#' Mardia, K. and Jupp, P. (1999). Directional Statistics, Section 7.4, John Wiley and Sons, England.
#'
#' @returns A list with class `"htest"`
#' @family Two-Sample-Tests
#' @export
#'
#' @examples
#' set.seed(20250411)
#' x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
#' x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#' equal_kappa_test(list(x1, x2))
#'
#' x1 <- rvm(100, mean = 90, kappa = 20)
#' x2 <- rvm(100, mean = 90, kappa = 10)
#' x3 <- rvm(100, mean = 90, kappa = 1)
#' equal_kappa_test(list(x1, x2, x3))
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' equal_kappa_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
equal_kappa_test <- function(samples, method = NULL, alpha = 0.05, axial = TRUE) {
  level <- 1 - alpha
  dname <- paste(names(x), collapse = " and ")

  f <- if (isTRUE(axial)) 2L else 1L
  samples <- lapply(samples, function(x){
    x * f
  })


  n_i <- lapply(samples, length) |> unlist()
  q   <- length(n_i)

  p <- 2
  nu  <- (p - 1) * (n_i - 1)

  R <- lapply(samples, mean_resultant_length) |> unlist()

  dev <- n_i - R
  # 2 kappa (n_i - R_i) ~ chi2(nu_i)  =>  kappa_i = nu_i / (2 (n_i - R_i))
  kappa <- setNames(nu / (2 * dev), paste0("kappa_", seq_len(q)))

  if (is.null(method)) method <- if (q == 2L) "f" else "bartlett"

  if (method == "f") {
    if (q != 2L) stop("The F test compares two samples")
    F <- (dev[1L] / nu[1L]) / (dev[2L] / nu[2L])
    p_value <- 2 * min(stats::pf(F, nu[1L], nu[2L]),
                       stats::pf(F, nu[1L], nu[2L], lower.tail = FALSE))
    return(structure(
      list(
        statistic   = c(F = F),
        parameter   = c(`num df` = nu[1L], `denom df` = nu[2L]),
        p.value     = min(p_value, 1),
        estimate    = c(`ratio of concentrations` = F),   # kappa_2 / kappa_1
        null.value  = c(`ratio of concentrations` = 1),
        alternative = "two.sided",
        method      = "F test of equal concentrations",
        data.name   = dname,
        kappa       = kappa,
        level       = level,
        source      = "Mardia & Jupp (2000), eq. 10.5.3"
      ),
      class = "htest"
    ))
  }

  total <- sum(nu)
  B <- total * log(sum(dev) / total) - sum(nu * log(dev / nu))
  B <- B / (1 + (sum(1 / nu) - 1 / total) / (3 * (q - 1)))
  structure(
    list(
      statistic   = c(`Bartlett's K-squared` = B),
      parameter   = c(df = q - 1),
      p.value     = stats::pchisq(B, q - 1, lower.tail = FALSE),
      estimate    = kappa,
      method      = "Bartlett test of equal concentrations",
      data.name   = dname,
      critical    = stats::qchisq(level, q - 1),
      level       = level,
      source      = "Mardia & Jupp (2000), eq. 10.6.21"
    ),
    class = "htest"
  )
}



## AR test ---------------------------------------------------------------------
ar_test_statistic <- function(s1, s2) {
  diffs <- outer(s1, s2, function(a, b) pi - abs(pi - abs(a - b)))
  sum(diffs)
}

#' Angular Randomisation Test of Homogeneity
#'
#' Performs Angular Randomisation Test for homogeneity on two samples of
#' circular data after Ruxton et al. (2023). Test whether samples of directions share the mean direction
#' P-values are estimated using permutation.
#'
#' @inheritParams watson_two_test
#' @param alpha (optional) numeric. Significance level of the test (values between 0 and 1).
#' @param axial logical. Whether the data are axial, i.e. \eqn{\pi}{pi}-periodical
#' (`TRUE`) or directional, i.e. \eqn{2 \pi}{2pi}-periodical (`FALSE`, the default).
#'
#' @details
#' **Null Hypothesis** (\eqn{H_0}{H₀}): The two circular samples share an identical underlying
#' probability distribution.
#'
#' **Alternative Hypothesis** (\eqn{H_{1}}{H1}): The two samples come from different distributions.
#'
#' ## Interpretation
#'
#' * Small p-value (\eqn{p < \alpha}, e.g., <0.05): Reject the null hypothesis.
#' This indicates strong evidence that the two samples come from different
#' circular distributions (differing in central tendency/mean direction or shape).
#'
#' * Large p-value (\eqn{p \ge \alpha}): Fail to reject the null hypothesis;
#' there is insufficient evidence to claim the two circular samples differ.
#'
#' @note
#' **Concentration Differences**: The test can suffer from markedly lower statistical
#' power if the underlying unimodal distributions differ by concentration
#' (dispersion/spread) rather than location—especially with small, uneven sample
#' sizes where the smaller sample comes from the more concentrated distribution.
#'
#' **Axial/Multimodal Data**: ART performs poorly and loses power when applied to
#' axially symmetric or symmetrically multimodal distributions.
#'
#'
#' @references Ruxton, G.D., Malkemper, E.P. & Landler, L. Evaluating the power
#' of a recent method for comparing two circular distributions: an alternative
#' to the Watson U2 test. Sci Rep 13, 10007 (2023). https://doi.org/10.1038/s41598-023-36960-1
#'
#'#' @returns A list with class `"htest"` containing the following components:
#' \describe{
#' \item{`statistic`}{the value of the G statistic}
#' \item{`parameter`}{the number of data in x and y, and the number of permutations}
#' \item{`p.value`}{the p-value for the test}
#' \item{`alternative`}{a character string describing the alternative hypothesis.}
#' \item{`method`}{a character string indicating what type of test was performed.}
#' \item{`alpha`}{the significance level}
#' \item{`reject`}{logical. whether to reject the null-hypothesis based on the significance level.}
#' \item{`data.name`}{a character string giving the name(s) of the data.}
#' }
#'
#' @export
#'
#' @family Two-Sample-Tests
#'
#' @examples
#' set.seed(20250411)
#' x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
#' x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#' ar_test(x1, x2)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' ar_test(sa.por$azi.PoR, rvm(100, 135, 10), axial = TRUE, alpha = 0.05)
ar_test <- function(x, y, n_perm = 1000L, axial = FALSE, alpha = NULL) {
  f <- if (isTRUE(axial)) 2 else 1

  samp1 <- deg2rad(((f * stats::na.omit(x)) %% 360))
  samp2 <- deg2rad(((f * stats::na.omit(y)) %% 360))

  n1 <- length(samp1)
  n2 <- length(samp2)
  n <- n1 + n2
  combined <- c(samp1, samp2)

  Gstat <- ar_test_statistic(samp1, samp2)

  # Permutation
  nxtrm <- sum(vapply(seq_len(n_perm), function(r) {
    perm <- sample.int(n)
    Grand <- ar_test_statistic(
      combined[perm[seq_len(n1)]],
      combined[perm[seq(n1 + 1, n)]]
    )
    Grand >= Gstat
  }, FUN.VALUE = logical(1))) + 1L

  pvalue <- nxtrm / (n_perm + 1)
  reject <- NULL
  if(!is.null(alpha)) {
    stopifnot(alpha >= 0  & alpha <= 1)
    reject <- pvalue < alpha
  }

  structure(
    list(
      method = "Angular Randomisation Test of Homogeneity",
      alternative = "Samples come from different distributions",
      statistic = c(G = Gstat),
      p.value = pvalue,
      alpha = alpha,
      parameter   = c(nx = n1, ny = n2, n_perm = n_perm),
      data.name = paste(deparse1(substitute(x)), "and", deparse1(substitute(y))),
      reject = reject
    ),
    class = 'htest'
  )
}

#' Wallraff Test of Angular Distances
#'
#' Performs the Wallraff test of angular distances or angular dispersion around the mean.
#'
#' @inheritParams equal_kappa_test
#' @param ref a vector of angles used as reference to compute the angular distances from, in each group.
#' It should contain as many elements as there are groups, *in the same order*.
#' The order is the order of the elements of the list `x`.
#'
#' If `ref` has less elements than the number of groups (typically one: a common reference for all groups), it is recycled to match the number of groups.
#'
#' If `ref` is `NULL` (the default), the mean angle of each group will be used as reference. In this situation, the Wallraff test becomes a comparison of angular dispersion around the mean
#'
#' @details
#' The Wallraff test of angular distances between two or more groups is performed
#' and the results are printed. The null hypothesis is that distances are equal across groups.
#'
#' The test proceeds by computing the angular distances from a reference angle,
#' in each group. The angular distance between two angles is the circular range.
#' Then the distances are compared with a usual rank sum test (Kruskal-Wallis, [stats::kruskal.test]).
#' When there are only two groups,
#' the Wilcoxon-Mann-Whitney test could be used but [stats::wilcox.test] without continuity
#' correction for the p-value is equivalent to [stats::kruskal.test] so only [stats::kruskal.test] is used here.
#'
#' The Wallraff test is most frequently used to compare angular dispersion
#' around the mean, between samples. In this case, the reference angle is the
#' mean angle of each sample. This is the default here, when no reference angles are provided.
#'
#'
#' @inherit circular::wallraff.test return
#' @inherit circular::wallraff.test references
#'
#' @family Two-Sample-Tests
#' @export
#'
#' @examples
#' set.seed(20250411)
#' x1 <- c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120)
#' x2 <- c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165)
#' walraff_test(list(x1, x2))
#'
#' x1 <- rvm(100, mean = 90, kappa = 20)
#' x2 <- rvm(100, mean = 90, kappa = 10)
#' x3 <- rvm(100, mean = 90, kappa = 1)
#' walraff_test(list(x1, x2, x3))
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' walraff_test(list(sa.por$azi.PoR, rvm(100, 135, 10)))
walraff_test <- function(samples, ref = NULL, axial = TRUE){
    if (!is.list(samples) || length(samples) < 2L) {
      stop("`samples` must be a list of at least two samples")
    }
    k  <- length(samples)
    nm <- names(samples)
    if (is.null(nm) || any(!nzchar(nm))) nm <- as.character(seq_len(k))
    nm <- make.unique(nm)

    samples <- lapply(samples, deg2rad)
    f <- if (isTRUE(axial)) 2 else 1
    x <- f * unlist(samples, use.names = FALSE)          # radians, doubled if axial
    g <- rep.int(seq_len(k), lengths(samples))
    ok <- is.finite(x)
    x <- x[ok]
    g <- g[ok]

    if (is.null(ref)) {
      SC <- rowsum(cbind(sin(x), cos(x)), g)             # k x 2 resultants
      if (any(rowSums(SC^2) < 1e-12)) {
        warning("Mean direction undefined for a sample with zero resultant")
      }
      centre <- atan2(SC[, 1L], SC[, 2L])
    } else {
      if (!length(ref) %in% c(1L, k)) {
        stop("`ref` must have length 1 or one value per sample")
      }
      centre <- rep_len(f * as.numeric(ref), k)
    }

    # angular distance to the centre, in [0, pi] (doubled) -> [0, pi/2] (axial)
    d <- abs(x - centre[g]) %% (2 * pi)
    d <- pmin(d, 2 * pi - d) / f

    res <- stats::kruskal.test(d, factor(g, levels = seq_len(k), labels = nm))
    res$method <- if (is.null(ref)) {
      "Wallraff rank sum test of angular distance to the sample mean directions"
    } else {
      "Wallraff rank sum test of angular distance to the reference direction"
    }
    res$data.name <- paste(nm, collapse = " and ")
    res
}
