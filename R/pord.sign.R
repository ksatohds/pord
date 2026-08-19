# pord.sign.R -- the other direction: is y falling short of x?
# pord.sign, pord.compare
#
# The primary test (pord.test) has independence as its null.  Its lower tail
# does NOT answer "is y falling short of x": under that null D = n - T, so a
# large D simply means negative association.  The natural null for the
# shortfall question is symmetry, P(y < x) = P(y > x).  Conditioning on the
# pair totals n_ij + n_ji removes the diagonal probabilities (which are free
# parameters under symmetry), and the lower-triangle count is then binomial
# with probability one half.  That is the sign test: discarding ties is a
# consequence of the conditioning, not an arbitrary choice.
#
# A ceiling in x is NOT a bias here.  Symmetry implies equal margins, so a
# concentration of x at the top with y below it is exactly the asymmetry the
# test is meant to detect.  (Versions up to 0.1.0 offered a `conditional`
# variant restricting to 1 < x < K; that selection is asymmetric in (x, y),
# so under symmetry the discordant probability is no longer one half and the
# reported p-values were wrong.  The argument has been removed.)

#' @title Sign test for paired ordinal shortfall
#' @description
#' \code{pord.sign} tests whether \eqn{y} falls below \eqn{x} more often than it
#' rises above it.  Conditioning on the number of discordant pairs leaves a
#' binomial distribution, so the diagonal (ties) drops out and the test is the
#' classical sign test.
#'
#' The null it actually tests is \emph{directional balance},
#' \eqn{P(Y < X) = P(Y > X)}.  Symmetry of the table (\eqn{\pi_{ij} = \pi_{ji}})
#' implies that null but is not implied by it: the test pools every off-diagonal
#' cell into two totals, so it is blind to a table in which one pair of cells
#' leans one way and another pair leans back.  With \eqn{n_{12} = 30, n_{21} = 0}
#' and \eqn{n_{34} = 0, n_{43} = 30} it returns \eqn{p = 1}, while Bowker's test
#' of symmetry gives \eqn{p < 10^{-13}}.  Use a Bowker or exact conditional
#' symmetry test when symmetry itself is the question; \code{pord.sign} answers
#' the directional one.
#'
#' A ceiling in \code{x} does not bias this test.  The null of symmetry
#' implies equal margins for \eqn{x} and \eqn{y}; when \eqn{x} piles up at the
#' top of the scale while \eqn{y} does not, that is the asymmetry being
#' tested, not an artefact.  A significant result on (nearly) every item is
#' therefore a finding; to rank items or domains by shortfall, use
#' [pord.compare()] or a heterogeneity test across items.
#'
#' The sign test is the method recommended for paired ordinal data by
#' Svensson (2001), who notes that the Wilcoxon signed-rank test is
#' \emph{not} appropriate because it ranks differences between ordinal ratings.
#'
#' @param x Integer vector of ratings in \code{1..K}, or a square table
#'   (in which case \code{y} is omitted).
#' @param y Integer vector of the same length as \code{x}.
#' @param K Number of scale points.  Taken from the data when \code{NULL}.
#' @param alternative \code{"greater"} (the default; \eqn{y} falls short more
#'   often than it exceeds), \code{"less"}, or \code{"two.sided"}.
#'
#' @return A list of class \code{"pord.sign"} with the counts below, tied and
#'   above, the proportion falling short, the proportion among discordant pairs,
#'   the exact binomial p-value, and diagnostics on how much of the sample sits
#'   at the ceiling.
#' @references
#' Svensson, E. (2001). Guidelines to statistical evaluation of data from rating
#'   scales and questionnaires. \emph{Journal of Rehabilitation Medicine}
#'   \strong{33}, 47--48.
#' @seealso [pord.test()], [pord.compare()], [pord.table()]
#' @examples
#' set.seed(1)
#' x <- sample(1:4, 96, TRUE, c(.05, .15, .35, .45))
#' y <- pmax(1, pmin(4, x - rbinom(96, 1, .35)))
#' pord.sign(x, y)
#' @export
pord.sign <- function(x, y = NULL, K = NULL,
                      alternative = c("greater", "less", "two.sided")) {
  alternative <- match.arg(alternative)
  d <- .pord.pair(x, y, K)
  lo <- sum(d$y < d$x); ti <- sum(d$y == d$x); hi <- sum(d$y > d$x)
  bt <- if (lo + hi > 0)
    stats::binom.test(lo, lo + hi, 0.5, alternative = alternative) else
      list(p.value = NA_real_)
  structure(list(
    below = lo, tied = ti, above = hi, n = d$n,
    prop.below = lo / d$n,
    prop.discordant = if (lo + hi > 0) lo / (lo + hi) else NA_real_,
    p.value = bt$p.value, alternative = alternative, K = d$K,
    prop.ceiling = mean(d$x == d$K), prop.floor = mean(d$x == 1)),
    class = "pord.sign")
}

#' @export
print.pord.sign <- function(x, ...) {
  cat("\n\tSign test for paired ordinal shortfall (exact test of symmetry)\n\n")
  cat(sprintf("n = %d\n", x$n))
  cat(sprintf("y < x: %d, y == x: %d, y > x: %d\n", x$below, x$tied, x$above))
  cat(sprintf("falling short: %.1f%% of all, %.1f%% of discordant pairs\n",
              100 * x$prop.below, 100 * x$prop.discordant))
  cat(sprintf("x at the ceiling: %.1f%%\n", 100 * x$prop.ceiling))
  cat(sprintf("alternative: %s\np-value = %s\n\n", x$alternative,
              format.pval(x$p.value)))
  invisible(x)
}

#' @title Pairwise comparison of shortfall between domains of items
#' @description
#' \code{pord.compare} compares, within respondents, how often \eqn{y} falls
#' short of \eqn{x} between domains of items -- for example between functional
#' domains made up of different numbers of items.  Each respondent contributes
#' the proportion of items in a domain for which \eqn{y < x}; domains are then
#' compared by the sign test on the direction of the within-respondent
#' difference, so no distance between ordinal categories is used.
#'
#' An overall Friedman test across domains, with Kendall's W, is reported
#' alongside the pairwise comparisons.
#'
#' @param X Matrix or data frame of \code{x} ratings, **respondents in rows and
#'   items in columns**.
#' @param Y Matrix or data frame of \code{y} ratings, same shape as \code{X}.
#' @param domain Vector of length \code{ncol(X)} assigning each item to a
#'   domain (subscale).  Note that this groups the *columns*; in the
#'   terminology of [stats::friedman.test()] the columns are its \code{groups}
#'   and the respondents in the rows are its \code{blocks}.
#' @param p.adjust.method Method passed to [stats::p.adjust()].
#'
#' @return A list of class \code{"pord.compare"} with the domain-level
#'   proportions and mean ranks, the Friedman test and Kendall's W, and a data
#'   frame of all pairwise sign tests.  Only respondents with complete data on
#'   every item are used, since the comparison is within respondents.
#' @seealso [pord.sign()], [pord.items()]
#' @examples
#' set.seed(1)
#' X <- matrix(sample(1:4, 40 * 6, TRUE, c(.05, .15, .35, .45)), 40, 6)
#' Y <- pmax(pmin(X - matrix(rbinom(40 * 6, 1, .4), 40, 6), 4), 1)
#' pord.compare(X, Y, domain = rep(c("A", "B"), each = 3))
#' @export
pord.compare <- function(X, Y, domain, p.adjust.method = "holm") {
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (!identical(dim(X), dim(Y))) stop("'X' and 'Y' must have the same shape")
  if (length(domain) != ncol(X)) stop("'domain' must have one entry per column")
  ok <- stats::complete.cases(X, Y)
  X <- X[ok, , drop = FALSE]; Y <- Y[ok, , drop = FALSE]
  n <- nrow(X)
  if (n < 2) stop("fewer than two complete respondents")
  lev <- unique(domain)
  Z <- vapply(lev, function(b) rowMeans((Y[, domain == b, drop = FALSE] <
                                         X[, domain == b, drop = FALSE]) * 1),
              numeric(n))
  colnames(Z) <- lev
  rk <- t(apply(Z, 1, rank))
  fr <- stats::friedman.test(Z)
  G <- length(lev)
  W <- 12 * sum((colSums(rk) - mean(colSums(rk)))^2) / (n^2 * (G^3 - G))
  pr <- utils::combn(G, 2)
  res <- do.call(rbind, lapply(seq_len(ncol(pr)), function(k) {
    a <- lev[pr[1, k]]; b <- lev[pr[2, k]]
    dd <- Z[, a] - Z[, b]
    hi <- sum(dd > 0); lo <- sum(dd < 0); ti <- sum(dd == 0)
    bt <- if (hi + lo > 0) stats::binom.test(hi, hi + lo, 0.5) else
      list(p.value = NA_real_)
    data.frame(comparison = paste(a, "vs", b), first.higher = hi,
               tied = ti, second.higher = lo,
               prop.discordant = if (hi + lo > 0) hi / (hi + lo) else NA_real_,
               p.value = bt$p.value, stringsAsFactors = FALSE)
  }))
  res$p.adjusted <- stats::p.adjust(res$p.value, p.adjust.method)
  res <- res[order(res$p.adjusted), ]
  structure(list(
    domains = data.frame(domain = lev, items = as.integer(table(domain)[lev]),
                        prop.short = colMeans(Z), mean.rank = colMeans(rk),
                        stringsAsFactors = FALSE),
    friedman = fr, kendall.W = W, pairwise = res, n = n,
    p.adjust.method = p.adjust.method), class = "pord.compare")
}

#' @export
print.pord.compare <- function(x, ...) {
  cat("\n\tBlockwise comparison of paired ordinal shortfall\n\n")
  cat(sprintf("n = %d complete respondents\n\n", x$n))
  print(x$domains[order(-x$domains$prop.short), ], row.names = FALSE, digits = 3)
  cat(sprintf("\nFriedman chi-squared = %.2f, df = %d, p = %s; Kendall's W = %.3f\n\n",
              x$friedman$statistic, x$friedman$parameter,
              format.pval(x$friedman$p.value), x$kendall.W))
  cat(sprintf("pairwise sign tests (%s-adjusted), %d of %d significant:\n",
              x$p.adjust.method, sum(x$pairwise$p.adjusted < 0.05),
              nrow(x$pairwise)))
  print(x$pairwise, row.names = FALSE, digits = 3)
  cat("\n")
  invisible(x)
}
