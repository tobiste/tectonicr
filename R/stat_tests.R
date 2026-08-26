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
#' @family Tests
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
#' @param alpha Significance level of the test. Valid levels are `0.01`, `0.05`, and `0.1`.
#' This argument may be omitted (`NULL`, the default), in which case, a range for the p-value will be returned.
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
#'  * Large p-value (p ≥ 0.05): Fail to reject the null hypothesis. There is not
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
#' @family Tests
#'
#' @name rayleigh-test
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' rayleigh_test(pidgeon_homing, axial = FALSE) # Do not reject null hypothesis.
#' # R = 0.22; stat = 0.497, p = 0.62
#' rayleigh_test_perm(pidgeon_homing, axial = FALSE)
#'
#' # Example data from Davis (1986), pp. 316
#' finland_striae <- c(
#'   23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
#'   113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
#'   132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
#'   165, 171, 172, 179, 181, 186, 190, 212
#' )
#' rayleigh_test(finland_striae, axial = FALSE) # reject null hypothesis
#' rayleigh_test_perm(finland_striae, axial = FALSE) # reject null hypothesis
#'
#' rayleigh_test(finland_striae, mu = 105, axial = FALSE) # reject null hypothesis
#' rayleigh_test_perm(finland_striae, mu = 105, axial = FALSE) # reject null hypothesis
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
    min(max(P * temp, 0), 1)     # removed the dead "P * temp" line above this, which had no effect
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
  f <- if (isTRUE(axial)) 2 else 1
  x <- stats::na.omit(x)
  xf <- (x * f) %% 360
  n <- length(x)

  if(is.null(mu)){
  stat <- mean_resultant_length(xf, na.rm = FALSE)

  null_dist <- vapply(seq_len(n_perm), function(i) {
    rnd <- runif(n, 0, 360)
    mean_resultant_length(rnd)
  }, FUN.VALUE = numeric(1))

  } else {
    mu <- (mu * f) %% 360
    xmu <- xf - mu

    stat <- (sum(cosd(xmu))) / n

    null_dist <- vapply(seq_len(n_perm), function(i) {
      rnd <- runif(n, 0, 360)
      (sum(cosd(rnd))) / n
    }, FUN.VALUE = numeric(1))
  }

  p.value <- (sum(null_dist >= stat) + 1L) / (n_perm + 1)

  return(list(statistic = stat, p.value = p.value))
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
#' @family Tests
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
#'  The **Null Hypothesis** (\eqn{H_0}{H₀}): The data are distributed completely uniformly
#'  (randomly and evenly) around the circle.
#'
#'  The **Alternative Hypothesis** (\eqn{H_1}{H₁}): The data are not uniform and
#'  show a preference, clustering, or pattern somewhere on the circle.
#'
#'  The Test Statistic (V or \eqn{D^{+} + D^{-}}{D⁺ + D⁻}): It measures the
#'  greatest positive and negative differences between your data's empirical
#'  cumulative distribution and a theoretical uniform distribution.
#'
#' ## Interpreting the Results
#'
#' * High Test Statistic / Low p-value (\eqn{p < \alpha}{p < α}, typically 0.05):
#' You reject the null hypothesis. This means your data are not uniform; they
#' have a significant preferred direction, grouping, or non-random pattern on the circle.
#'
#' * Low Test Statistic / High p-value (\eqn{p \ge 0.05}{p ≥ 0.05}): You fail to
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
#' @returns list containing the test statistic `statistic` and the significance
#' level `p.value`.
#'
#' @family Tests
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
  f <- if (isTRUE(axial)) 2 else 1

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

    if (isFALSE(quiet)) {
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

#' Watson's \eqn{U^2}{U²} Test for Goodness-of-Fit against known Distribution
#'
#' A non-parametric statistical test used for circular data to determine whether
#' a sample fits a specified theoretical distribution
#'
#' @inheritParams rayleigh_test
#' @param alpha Significance level of the test. Valid levels are `0.01`, `0.05`,
#' and `0.1`.
#' This argument may be omitted (`NULL`, the default), in which case, a range
#' for the p-value will be returned.
#' @param dist Distribution to test for. The default, `"uniform"`, is the
#' circular uniform distribution. `"vonmises"` tests the von Mises distribution.
#'
#' @returns list containing the test statistic `statistic`, the significance
#' level `p.value`, the critical value `critical.value`, whether to `reject` the null hypothesis,
#' the significance level `alpha`, the tested distribution `dist`, and the number of data `n`
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
#' Watson's goodness-of-fit/homogeneity tables at your chosen significance level (\eqn{\alpha}{α}, commonly set to 0.05),
#' or check the resulting p-value:
#'
#'  * If \eqn{U^2_{\text{calculated}} > U^2_{\text{critical}}} (or p < \eqn{\alpha}{α}):
#'  Reject the null hypothesis (\eqn{H_0}{H0}).
#'  Conclude that the data significantly deviates from the theoretical distribution.
#'
#'  * If \eqn{U^2_{\text{calculated}} \le U^2_{\text{critical}}} (or \eqn{p \ge \alpha}{p ≥ α}):
#'  Fail to reject the null hypothesis (\eqn{H_0}{H0}).
#'  There is not enough evidence to claim the data deviates from the expected model.
#'
#' @note Watson's test statistic is a rotation-invariant Cramer - von Mises test.
#' non-parametric, rank-based alternative to one-sample
#'
#' @references Mardia and Jupp (1999). Directional Statistics. John Wiley and
#' Sons.
#'
#' @export
#'
#' @family Tests
#' @seealso [vonmises] and [cunif]
#'
#' @examples
#' # Example data from Mardia and Jupp (1999), pp. 93
#' pidgeon_homing <- c(55, 60, 65, 95, 100, 110, 260, 275, 285, 295)
#' watson_test(pidgeon_homing, axial = FALSE, alpha = .05)
#'
#' # San Andreas Fault Data:
#' data(san_andreas)
#' data("nuvel1")
#' PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
#' sa.por <- PoR_shmax(san_andreas, PoR, "right")
#' watson_test(sa.por$azi.PoR, alpha = .05)
#' watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
#' watson_test(sa.por$azi.PoR, alpha = .05, dist = "vonmises")
watson_test <- function(x, alpha = NULL, dist = c("uniform", "vonmises"), axial = TRUE, quiet = FALSE) {
  if (is.null(alpha)) alpha <- 0
  allowed_alphas <- c(0, 0.01, 0.025, 0.05, 0.1)
  if (!(alpha %in% allowed_alphas)) {
    stop("'alpha' must be one of: 0, 0.01, 0.025, 0.05, 0.1")
  }
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
    u <- sort(deg2rad(x2)) / (2 * pi)
    u.bar <- mean(u)
    i <- seq_len(n)
    u2 <- sum((u - u.bar - (i - .5) / n + .5)^2) + 1 / (12 * n)
    statistic <- (u2 - 0.1 / n + 0.1 / (n^2)) * (1 + 0.8 / n)

    ## crits correspond to alpha = 0.01, 0.025, 0.05, 0.1 respectively
    crits       <- c(0.267, 0.221, 0.187, 0.152)
    crit_alphas <- c(0.01, 0.025, 0.05, 0.1)

    if (n < 8) {
      p.value <- NA_character_
      critical.value <- NA_real_
      reject <- NA
    } else if (alpha == 0) {
      thresholds <- c(rev(crits), Inf)                 # c(0.152, 0.187, 0.221, 0.267, Inf)
      messages <- c(
        "P-value > 0.10",
        "0.05 < P-value < 0.10",
        "0.025 < P-value < 0.05",
        "0.01 < P-value < 0.025",   # fixed: was "0.001" -- table only resolves to alpha=0.01
        "P-value < 0.01"            # fixed: was "0.001"
      )
      p.value <- messages[findInterval(statistic, thresholds) + 1]   # fixed: was missing "+ 1"
      critical.value <- NA_real_
      reject <- NA
    } else {
      critical.value <- crits[crit_alphas == alpha]
      reject <- statistic > critical.value
      p.value <- NA_character_
      if (!quiet) {
        message(if (reject) "Reject Null Hypothesis" else "Do Not Reject Null Hypothesis")
      }
    }
  } else {
    ## ---- goodness-of-fit against von Mises ----
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
      if (!quiet) {
        message(if (reject) "Reject Null Hypothesis" else "Do Not Reject Null Hypothesis")
      }
    }
  }

  list(
    statistic      = statistic,
    p.value        = p.value,        # descriptive bracket, only set when alpha == 0
    critical.value = critical.value, # tabulated threshold, only set when alpha != 0
    reject         = reject,         # logical decision, only set when alpha != 0
    alpha          = alpha,
    dist           = dist,
    n              = n
  )
}

#' Watson's Two-Sample Test of Homogeneity
#'
#' Performs Watson's test for homogeneity on two samples of circular data.
#' `watson_two_test_perm()` uses permutation to estimate p-values.
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
#' *  Low p-value (\eqn{p \le \alpha}{p ≤ α}, usually 0.05): Reject the null hypothesis. Conclude that
#' the two groups have significantly different circular distributions.
#' * High p-value (\eqn{p > \alpha}{p > α}): Fail to reject the null hypothesis. There is not
#' enough evidence to say the two groups are distributed differently around the circle.
#'
#' @note Critical values for the test statistic are obtained using the asymptotic
#' distribution of the test statistic. It is recommended to use the obtained
#' critical values and ranges for p-values only for combined sample sizes in
#' excess of 17. Tables are available for smaller sample sizes and can be found
#' in Mardia (1972) for instance.
#'
#' @name watson_two_sample
#' @family Tests
#'
#' @returns list. Watson's two-sample test of homogeneity is performed, and the
#' results are printed. If alpha is specified and non-zero, the test statistic
#' is printed along with the critical value and decision. If alpha is omitted,
#' the test statistic is printed and a range for the p-value of the test is given.
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
#' watson_two_test_perm(sa.por$azi, rvm(100, 135, 10), alpha = 0.05)
NULL

#' @export
#' @rdname watson_two_sample
watson_two_test <- function(x, y, alpha = NULL, axial = TRUE, quiet = FALSE) {
  if (is.null(alpha)) alpha <- 0
  allowed_alphas <- c(0, 0.001, 0.01, 0.05, 0.1)
  if (!(alpha %in% allowed_alphas)) {
    stop("'alpha' must be one of: 0, 0.001, 0.01, 0.05, 0.1")
  }

  # --- Preprocessing ---
  x <- stats::na.omit(x)
  y <- stats::na.omit(y)

  f <- if (axial) 2L else 1L
  x <- deg2rad((f * x) %% 360)
  y <- deg2rad((f * y) %% 360)

  # --- U2 statistic ---
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  labels <- c(rep(1L, n1), rep(2L, n2))
  labels <- labels[order(c(x, y))]

  cx <- cumsum(labels == 1L)
  cy <- cumsum(labels == 2L)

  d <- cy / n2 - cx / n1
  dbar <- mean.default(d)
  u2 <- (n1 * n2) / n^2 * sum((d - dbar)^2)

  if (n < 18) warning("Total sample size < 18: consult tabulated critical values")

  # --- Critical values (parallel to allowed_alphas) ---
  crits <- c(NA_real_, 0.385, 0.268, 0.187, 0.152)

  if (alpha == 0) {
    # Return a descriptive p-value bracket
    thresholds <- c(0.152, 0.187, 0.268, 0.385)
    labels_p <- c(
      "P-value > 0.10",
      "0.05 < P-value < 0.10",
      "0.01 < P-value < 0.05",
      "0.001 < P-value < 0.01",
      "P-value < 0.001"
    )
    p.value <- labels_p[findInterval(u2, thresholds, rightmost.closed = TRUE) + 1L]
  } else {
    crit <- crits[match(alpha, allowed_alphas)]
    if (!quiet) {
      message(if (u2 > crit) "Reject null hypothesis" else "Do not reject null hypothesis")
    }
    p.value <- crit
  }

  list(statistic = u2, p.value = p.value)
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

  list(statistic = Gstat, p.value = pvalue, alpha = alpha, reject = reject)
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
#' **Null Hypothesis (H₀)** The samples come from identical populations (meaning
#' both the mean direction and the dispersion/variance are homogeneous across groups).
#'
#' **Alternative Hypothesis (H₁)**: At least one sample comes from a different
#' population distribution, which can be due to a difference in the mean direction,
#' a difference in variance/concentration, or both.
#'
#' ## Interpretation
#'
#' **Test Statistic (W)**: This value follows an approximate chi-squared (\eqn{\Chi^2}{χ²}) distribution.
#' Higher values of W indicate larger discrepancies between the angular distributions of your groups.
#'
#' * If the p-value is less than your significance level (commonly α = 0.05), you **reject** the null hypothesis.
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
#' @family Tests
#'
#' @returns list
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

  list(statistic = Gstat, p.value = pvalue, alpha = alpha, reject = reject)
}


ar_test_statistic <- function(s1, s2) {
  diffs <- outer(s1, s2, function(a, b) pi - abs(pi - abs(a - b)))
  sum(diffs)
}

#' Angular Randomisation Test of Homogeneity
#'
#' Performs Angular Randomisation Test for homogeneity on two samples of
#' circular data after Ruxton et al. (2023).
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
#' @returns list containing the test statistic, the p-value, the significance value `alpha` and a logical decision whether to reject the null hypothesis or not.
#'
#' @export
#'
#' @family Tests
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

  return(list(statistic = Gstat, p.value = pvalue, alpha = alpha, reject = reject))
}
