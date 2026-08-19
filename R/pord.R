# pord.R -- core: exact conditional test for P(y >= x) on a square table
# pord.table, pord.moments, pord.test
#
# Setting: n subjects each give a pair (x_i, y_i) on the same K-point ordinal
# scale.  Let T = #{i : y_i >= x_i}, i.e. the sum of the diagonal and upper
# triangle of the K x K table with rows = x and columns = y.
#
# Null hypothesis: x and y are independent.  Conditioning on BOTH margins
# (row sums r, column sums s) gives the Fisher-Yates multiple hypergeometric
# distribution
#     P(N = m) = ( prod_k r_k! * prod_l s_l! ) / ( n! * prod_{k,l} m_kl! ),
# the same conditioning as Fisher's exact test.  T is a linear functional of
# that table.  For K = 2 the table has one degree of freedom and the test
# reduces exactly to Fisher's exact test; for K > 2 it does not collapse to a
# single cell, so T is genuinely a different statistic from the
# Fisher-Freeman-Halton statistic (which uses the table probability itself and
# detects departures in every direction).
#
# Three computations are provided and agree with each other:
#   exact        dynamic programming over all tables with the given margins
#   permutation  shuffle y against x
#   normal       closed-form mean/variance + continuity correction

## ---------------------------------------------------------------------
## internal helpers
## ---------------------------------------------------------------------

.pord.pair <- function(x, y, K = NULL) {
  if (is.table(x) || (is.matrix(x) && is.null(y))) {
    M <- as.matrix(x)
    if (nrow(M) != ncol(M)) stop("a square table is required")
    K <- nrow(M)
    idx <- which(M > 0, arr.ind = TRUE)
    cnt <- M[idx]
    x <- rep(idx[, 1], cnt); y <- rep(idx[, 2], cnt)
  }
  if (length(x) != length(y)) stop("'x' and 'y' must have the same length")
  ok <- !is.na(x) & !is.na(y)
  x <- as.integer(x[ok]); y <- as.integer(y[ok])
  if (!length(x)) stop("no complete pairs")
  if (is.null(K)) K <- max(c(x, y))
  if (K < 2) stop("the scale must have at least two levels")
  if (min(c(x, y)) < 1 || max(c(x, y)) > K)
    stop("'x' and 'y' must take integer values in 1..K")
  list(x = x, y = y, K = K, n = length(x),
       r = as.integer(table(factor(x, levels = seq_len(K)))),
       s = as.integer(table(factor(y, levels = seq_len(K)))))
}

## exact null distribution of T by dynamic programming.
## Cells are filled in row-major order.  The state is the vector of remaining
## column capacities; the budget left in the current row is determined by their
## sum, because the sum at the start of row k is a constant.  Table weights are
## proportional to prod 1/m_kl!, rescaled at every sub-step to avoid underflow.
.pord.exact.pmf <- function(r, s) {
  K <- length(r); n <- sum(r)
  if (sum(s) != n) stop("row and column sums must agree")
  P <- s + 1
  enc <- function(C) {
    if (is.null(dim(C))) C <- matrix(C, nrow = 1)
    key <- C[, 1]
    if (K >= 2) for (l in 2:K) key <- key + prod(P[seq_len(l - 1)]) * C[, l]
    key
  }
  rowStart <- n - c(0, cumsum(r)[-K])
  cvec <- matrix(as.numeric(s), nrow = 1)
  W <- matrix(0, nrow = 1, ncol = n + 1); W[1, 1] <- 1
  lfac <- lgamma(seq_len(max(c(s, r)) + 2))   # lfac[m+1] = log(m!)

  for (k in seq_len(K)) {
    if (r[k] == 0) next
    for (l in seq_len(K)) {
      budget <- r[k] - rowStart[k] + rowSums(cvec)
      capl   <- cvec[, l]
      restl  <- if (l < K) rowSums(cvec[, (l + 1):K, drop = FALSE]) else rep(0, nrow(cvec))
      mlo <- pmax(0, budget - restl)
      mhi <- pmin(capl, budget)
      valid <- mhi >= mlo
      if (!any(valid)) stop("no reachable state; check the margins")
      aK <- NULL; aC <- NULL; aW <- NULL
      for (m in 0:max(mhi)) {
        sel <- which(valid & mlo <= m & m <= mhi)
        if (!length(sel)) next
        cc <- cvec[sel, , drop = FALSE]; cc[, l] <- cc[, l] - m
        Wm <- W[sel, , drop = FALSE] * exp(-lfac[m + 1])
        if (l >= k && m > 0) {
          Z <- matrix(0, nrow(Wm), n + 1)
          Z[, (m + 1):(n + 1)] <- Wm[, 1:(n + 1 - m), drop = FALSE]
          Wm <- Z
        }
        aK <- c(aK, enc(cc)); aC <- rbind(aC, cc); aW <- rbind(aW, Wm)
      }
      uk   <- sort(unique(aK))
      W    <- rowsum(aW, aK, reorder = TRUE)
      cvec <- aC[match(uk, aK), , drop = FALSE]
      mx <- max(W); if (mx > 0) W <- W / mx
    }
  }
  if (nrow(W) != 1) stop("internal error: the final state is not unique")
  W[1, ] / sum(W[1, ])
}

## number of K x K tables with the given margins, on the log10 scale.
## Used to decide whether the exact computation is affordable.
.pord.ntables.log10 <- function(r, s) {
  K <- length(r); n <- sum(r)
  P <- s + 1
  enc <- function(C) {
    if (is.null(dim(C))) C <- matrix(C, nrow = 1)
    key <- C[, 1]
    if (K >= 2) for (l in 2:K) key <- key + prod(P[seq_len(l - 1)]) * C[, l]
    key
  }
  rowStart <- n - c(0, cumsum(r)[-K])
  cvec <- matrix(as.numeric(s), nrow = 1); lw <- 0
  for (k in seq_len(K)) {
    if (r[k] == 0) next
    for (l in seq_len(K)) {
      budget <- r[k] - rowStart[k] + rowSums(cvec)
      capl   <- cvec[, l]
      restl  <- if (l < K) rowSums(cvec[, (l + 1):K, drop = FALSE]) else rep(0, nrow(cvec))
      mlo <- pmax(0, budget - restl); mhi <- pmin(capl, budget)
      aK <- NULL; aC <- NULL; aW <- NULL
      for (m in 0:max(mhi)) {
        sel <- which(mhi >= mlo & mlo <= m & m <= mhi)
        if (!length(sel)) next
        cc <- cvec[sel, , drop = FALSE]; cc[, l] <- cc[, l] - m
        aK <- c(aK, enc(cc)); aC <- rbind(aC, cc); aW <- c(aW, lw[sel])
      }
      uk <- sort(unique(aK))
      lw <- as.numeric(tapply(aW, aK, function(v) {
        mx <- max(v); mx + log10(sum(10^(v - mx)))
      })[as.character(uk)])
      cvec <- aC[match(uk, aK), , drop = FALSE]
    }
  }
  lw[1]
}

## ---------------------------------------------------------------------
## public
## ---------------------------------------------------------------------

#' @title Cross-tabulate a pair of ordinal ratings
#' @description
#' \code{pord.table} builds the K x K table with rows given by \code{x} and
#' columns by \code{y}, and splits it into the three regions used throughout
#' the package: \code{y < x} (lower triangle), \code{y == x} (diagonal) and
#' \code{y > x} (upper triangle).  Expected counts under independence with the
#' observed margins are reported alongside.
#'
#' @param x Integer vector of ratings in \code{1..K}, or a square table
#'   (in which case \code{y} is omitted).
#' @param y Integer vector of the same length as \code{x}.
#' @param K Number of scale points.  Taken from the data when \code{NULL}.
#'
#' @return A list with
#' \item{table}{the K x K table of observed counts.}
#' \item{expected}{expected counts under independence with the observed margins.}
#' \item{regions}{a data frame giving observed and expected counts for
#'   \code{y < x}, \code{y == x} and \code{y > x}.}
#' \item{n, K, r, s}{sample size, number of scale points, row and column sums.}
#' @seealso [pord.test()], [pord.sign()]
#' @examples
#' set.seed(1)
#' x <- sample(1:4, 60, TRUE, c(.05, .15, .35, .45))
#' y <- pmax(1, pmin(4, x - rbinom(60, 1, .4)))
#' pord.table(x, y)$regions
#' @export
pord.table <- function(x, y = NULL, K = NULL) {
  d <- .pord.pair(x, y, K)
  O <- table(x = factor(d$x, levels = seq_len(d$K)),
             y = factor(d$y, levels = seq_len(d$K)))
  E <- outer(d$r, d$s) / d$n
  reg <- function(M) c(sum(M[lower.tri(M)]), sum(diag(M)), sum(M[upper.tri(M)]))
  ro <- reg(as.matrix(O)); re <- reg(E)
  structure(list(
    table = O, expected = E,
    regions = data.frame(
      region   = c("y < x", "y == x", "y > x"),
      observed = ro, expected = re, difference = ro - re),
    n = d$n, K = d$K, r = d$r, s = d$s), class = "pord.table")
}

#' @title Closed-form null moments of the superiority count
#' @description
#' \code{pord.moments} returns the mean and variance of
#' \eqn{T = \#\{i : y_i \ge x_i\}} under independence with both margins fixed.
#' Both depend only on the two marginal distributions, and are exact -- they
#' agree with the exact enumeration to machine precision.  No simulation and no
#' enumeration is involved, so this is available for any \code{K} and any
#' sample size.
#'
#' With \eqn{p_k = P(x = k)}, \eqn{q_l = P(y = l)},
#' \eqn{G(k) = P(y \ge k)}, \eqn{F(l) = P(x \le l)} and
#' \eqn{\theta = \sum_k p_k G(k)},
#' \deqn{E[T/n] = \theta,}
#' \deqn{Var[T/n] = \frac{1}{n-1}\sum_{k,l} p_k q_l
#'       \left[1\{l \ge k\} - G(k) - F(l) + \theta\right]^2 .}
#'
#' Because both moments are functions of the margins alone, **the null
#' distribution is not shared between items with different margins**.
#'
#' @param x Integer vector of ratings in \code{1..K}, or a square table.
#' @param y Integer vector of the same length as \code{x}.
#' @param K Number of scale points.  Taken from the data when \code{NULL}.
#'
#' @return A list with \code{theta} (the null mean of \eqn{T/n}, i.e. the value
#'   expected under independence), \code{mean} and \code{sd} of \eqn{T} in
#'   counts, and \code{n}.
#' @details
#' The closed-form mean and variance follow from the finite-population
#' (permutation) variance formula for a sum scored over pairs, obtained by
#' double centering; the general framework is that of Wald and Wolfowitz
#' (1944), Hoeffding (1951) and Strasser and Weber (1999).  The formula
#' itself is derived in this package rather than taken from a published
#' source, and is verified against the exact null distribution in the
#' package tests.
#' @references
#' Wald, A. and Wolfowitz, J. (1944). Statistical tests based on permutations of
#'   the observations. \emph{Annals of Mathematical Statistics} \strong{15}, 358--372.
#'   \doi{10.1214/aoms/1177731207}
#'
#' Hoeffding, W. (1951). A combinatorial central limit theorem.
#'   \emph{Annals of Mathematical Statistics} \strong{22}, 558--566.
#'   \doi{10.1214/aoms/1177729545}
#'
#' Strasser, H. and Weber, C. (1999). On the asymptotic theory of permutation
#'   statistics. \emph{Mathematical Methods of Statistics} \strong{8}, 220--250.
#' @seealso [pord.test()]
#' @examples
#' set.seed(1)
#' x <- sample(1:4, 96, TRUE); y <- sample(1:4, 96, TRUE)
#' pord.moments(x, y)
#' @export
pord.moments <- function(x, y = NULL, K = NULL) {
  d <- .pord.pair(x, y, K)
  K <- d$K; n <- d$n
  p <- d$r / n; q <- d$s / n
  G <- vapply(seq_len(K), function(k) sum(q[k:K]), numeric(1))
  Fx <- vapply(seq_len(K), function(l) sum(p[1:l]), numeric(1))
  th <- sum(p * G)
  V <- 0
  for (k in seq_len(K)) for (l in seq_len(K))
    V <- V + p[k] * q[l] * ((l >= k) - G[k] - Fx[l] + th)^2
  V <- V / (n - 1)
  list(theta = th, mean = n * th, sd = n * sqrt(V), n = n)
}

#' @title Exact conditional test for paired ordinal superiority
#' @description
#' \code{pord.test} tests whether the number of pairs with \eqn{y \ge x} exceeds
#' the value expected when \code{x} and \code{y} are independent, conditioning
#' on both margins of the K x K table.
#'
#' This is a generalization of Fisher's exact test.  It uses the same
#' conditioning -- the Fisher-Yates multiple hypergeometric distribution over
#' all tables with the observed margins -- and reduces exactly to Fisher's exact
#' test when the scale is collapsed to two levels.  It differs only in the test
#' statistic: the Fisher-Freeman-Halton extension of Fisher's test to K x K
#' tables uses the probability of the table itself and detects departures from
#' independence in any direction, whereas \code{pord.test} counts one direction
#' only and therefore answers a narrower, more demanding question.
#'
#' The construction -- an exact conditional test of independence that orders
#' the tables with fixed margins by a statistic chosen for the question at
#' hand -- follows Agresti and Wackerly (1977); exact conditional inference
#' for ordered categories is treated by Agresti, Mehta and Patel (1990).
#'
#' @param x Integer vector of ratings in \code{1..K}, or a square table
#'   (in which case \code{y} is omitted).
#' @param y Integer vector of the same length as \code{x}.
#' @param K Number of scale points.  Taken from the data when \code{NULL}.
#' @param method One of \code{"auto"}, \code{"exact"}, \code{"permutation"},
#'   \code{"normal"}.  \code{"auto"} uses the exact computation when the number
#'   of tables with the given margins is below \code{max.tables}, and the
#'   permutation otherwise.
#' @param alternative \code{"greater"} (the default; more pairs with
#'   \eqn{y \ge x} than independence gives), \code{"less"}, or
#'   \code{"two.sided"}.  Note that \code{"less"} tests for negative
#'   association, \emph{not} for \eqn{y} falling short of \eqn{x}; use
#'   [pord.sign()] for that question.
#' @param B Number of permutations when \code{method = "permutation"}.
#' @param max.tables Threshold on \eqn{\log_{10}} of the number of tables, used
#'   by \code{method = "auto"}.
#' @param continuity Apply a continuity correction when
#'   \code{method = "normal"}.  Strongly recommended: the statistic is
#'   integer-valued and its lattice spacing is a substantial fraction of the
#'   null standard deviation, so the uncorrected normal tail understates the
#'   p-value by roughly a factor of 1.5.
#'
#' @return A list of class \code{"pord.test"} with the observed count and
#'   proportion, the value expected under independence (\code{theta}), the
#'   excess, the null mean and standard deviation, the p-value, the method
#'   actually used, and -- for \code{method = "exact"} -- the complete null
#'   probability mass function.
#' @references
#' Agresti, A. and Wackerly, D. (1977). Some exact conditional tests of
#'   independence for r x c cross-classification tables.
#'   \emph{Psychometrika} \strong{42}, 111--125.
#'
#' Agresti, A., Mehta, C. R. and Patel, N. R. (1990). Exact inference for
#'   contingency tables with ordered categories. \emph{Journal of the American
#'   Statistical Association} \strong{85}, 453--458.
#'
#' Fisher, R. A. (1935). The logic of inductive inference.
#'   \emph{Journal of the Royal Statistical Society} \strong{98}, 39--82.
#'
#' Agresti, A., Wackerly, D. and Boyett, J. M. (1979). Exact conditional tests
#'   for cross-classifications: approximation of attained significance levels.
#'   \emph{Psychometrika} \strong{44}, 75--83.
#'
#' Phipson, B. and Smyth, G. K. (2010). Permutation p-values should never be
#'   zero. \emph{Statistical Applications in Genetics and Molecular Biology}
#'   \strong{9}, Article 39. \doi{10.2202/1544-6115.1585}
#' @seealso [pord.moments()], [pord.sign()], [pord.group()], [pord.table()]
#' @examples
#' set.seed(1)
#' x <- sample(1:4, 60, TRUE, c(.05, .15, .35, .45))
#' y <- pmax(1, pmin(4, x - rbinom(60, 1, .35)))
#' pord.test(x, y)
#'
#' # collapsing to two levels reproduces Fisher's exact test
#' x2 <- as.integer(x == 4); y2 <- as.integer(y == 4)
#' pord.test(x2 + 1L, y2 + 1L)$p.value
#' @export
pord.test <- function(x, y = NULL, K = NULL,
                      method = c("auto", "exact", "permutation", "normal"),
                      alternative = c("greater", "less", "two.sided"),
                      B = 10000, max.tables = 9, continuity = TRUE) {
  method <- match.arg(method); alternative <- match.arg(alternative)
  d <- .pord.pair(x, y, K)
  obs <- sum(d$y >= d$x)
  mo  <- pord.moments(d$x, d$y, d$K)

  if (method == "auto") {
    lg <- try(.pord.ntables.log10(d$r, d$s), silent = TRUE)
    method <- if (!inherits(lg, "try-error") && is.finite(lg) && lg <= max.tables)
      "exact" else "permutation"
  }

  pmf <- NULL
  if (method == "exact") {
    pmf <- .pord.exact.pmf(d$r, d$s)
    up <- sum(pmf[(obs + 1):(d$n + 1)])
    lo <- sum(pmf[1:(obs + 1)])
    p <- switch(alternative, greater = up, less = lo,
                two.sided = min(1, 2 * min(up, lo)))
  } else if (method == "permutation") {
    sim <- replicate(B, sum(sample(d$y) >= d$x))
    up <- (sum(sim >= obs) + 1) / (B + 1)
    lo <- (sum(sim <= obs) + 1) / (B + 1)
    p <- switch(alternative, greater = up, less = lo,
                two.sided = min(1, 2 * min(up, lo)))
  } else {
    cc <- if (continuity) 0.5 else 0
    zu <- (obs - cc - mo$mean) / mo$sd
    zl <- (obs + cc - mo$mean) / mo$sd
    up <- stats::pnorm(zu, lower.tail = FALSE)
    lo <- stats::pnorm(zl)
    p <- switch(alternative, greater = up, less = lo,
                two.sided = min(1, 2 * min(up, lo)))
  }

  structure(list(
    statistic = obs, proportion = obs / d$n,
    theta = mo$theta, excess = obs / d$n - mo$theta,
    null.mean = mo$mean, null.sd = mo$sd,
    z = (obs - mo$mean) / mo$sd,
    p.value = p, alternative = alternative, method = method,
    B = if (method == "permutation") B else NA_integer_,
    n = d$n, K = d$K, r = d$r, s = d$s, null.pmf = pmf),
    class = "pord.test")
}

#' @export
print.pord.test <- function(x, ...) {
  cat("\n\tExact conditional test for paired ordinal superiority\n\n")
  cat(sprintf("n = %d, scale points = %d, method = %s%s\n", x$n, x$K, x$method,
              if (!is.na(x$B)) sprintf(" (B = %s)", format(x$B, big.mark = ",")) else ""))
  cat(sprintf("count(y >= x) = %d (%.1f%%); expected under independence %.1f%%; excess %+.1f points\n",
              x$statistic, 100 * x$proportion, 100 * x$theta, 100 * x$excess))
  cat(sprintf("null mean = %.2f, null sd = %.3f, z = %.2f\n", x$null.mean, x$null.sd, x$z))
  cat(sprintf("alternative: %s\np-value = %s\n\n", x$alternative, format.pval(x$p.value)))
  invisible(x)
}

#' @export
print.pord.table <- function(x, ...) {
  cat("\n\tPaired ordinal ratings (rows = x, columns = y)\n\n")
  print(x$table)
  cat("\n")
  print(x$regions, row.names = FALSE, digits = 4)
  cat("\n")
  invisible(x)
}
