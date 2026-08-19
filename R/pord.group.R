# pord.group.R -- group comparison and multi-item wrappers
# pord.group, pord.items
#
# Group comparison permutes the group label across respondents while keeping
# each respondent's whole vector of items intact.  The permutation is therefore
# unaffected by within-respondent and between-item correlation.
#
# Two statistics answer different questions and are reported separately:
#   "excess"     per group, the observed proportion minus the value expected
#                under independence using that group's OWN margins.  This
#                adjusts for differences in the marginal distributions and asks
#                whether the strength of the association is common to the groups.
#   "proportion" the raw proportion, compared across groups item by item.  This
#                does not adjust for the margins, and is the quantity a plot of
#                group proportions actually displays.
#
# Pooling the groups for the primary analysis is justified by the first; the
# second is an exploratory scan for items whose level differs by group.

.pord.excess <- function(x, y, K) {
  n <- length(x)
  p <- as.numeric(table(factor(x, levels = seq_len(K)))) / n
  q <- as.numeric(table(factor(y, levels = seq_len(K)))) / n
  G <- vapply(seq_len(K), function(k) sum(q[k:K]), numeric(1))
  mean(y >= x) - sum(p * G)
}

.pord.chi2 <- function(nn, xx) {
  N <- sum(nn); X <- sum(xx)
  if (N == 0 || X == 0 || X == N) return(0)
  p <- X / N; e1 <- nn * p; e0 <- nn * (1 - p); s <- 0
  ok <- e1 > 0; s <- s + sum((xx[ok] - e1[ok])^2 / e1[ok])
  ok <- e0 > 0; s <- s + sum(((nn - xx)[ok] - e0[ok])^2 / e0[ok])
  s
}

#' @title Does the paired ordinal pattern differ between groups?
#' @description
#' \code{pord.group} tests whether the relation between \code{X} and \code{Y}
#' differs between groups of respondents, by permuting the group label while
#' keeping each respondent's whole vector of item responses intact.  Because
#' whole vectors are moved, the test is unaffected by within-respondent and
#' between-item correlation.
#'
#' Two statistics are computed, because they answer different questions.
#' \code{"excess"} subtracts, for each group, the value expected under
#' independence using that group's own margins, and so asks whether the
#' \emph{strength of the association} is common to the groups; this is what
#' justifies pooling the groups for [pord.test()].  \code{"proportion"}
#' compares the raw proportions item by item without adjusting for the margins,
#' and so corresponds to what a plot of group proportions displays.  A group can
#' differ on the second while agreeing on the first, when its marginal
#' distributions differ.
#'
#' @param X Matrix or data frame of \code{x} ratings, **respondents in rows and
#'   items in columns**.
#' @param Y Matrix or data frame of \code{y} ratings, same shape as \code{X}.
#' @param group Vector of group labels, **one per respondent** (i.e. one per
#'   row of \code{X}), such as occupation or site.  This follows the
#'   psychometric convention in which \code{group} is a person-level variable
#'   (as in \code{mirt::multipleGroup} and \code{lavaan::cfa}).  Note that it is
#'   *not* the \code{groups} of [stats::friedman.test()], which refers to the
#'   columns.
#' @param direction \code{"ge"} for \eqn{y \ge x} (the default) or \code{"lt"}
#'   for \eqn{y < x}.
#' @param K Number of scale points.  Taken from the data when \code{NULL}.
#' @param B Number of permutations.
#' @param p.adjust.method Method passed to [stats::p.adjust()] for the item-wise
#'   tests.
#'
#' @return A list of class \code{"pord.group"} with the group-level
#'   proportions and excesses, the global permutation p-value for each
#'   statistic, a max-statistic p-value that absorbs the multiplicity across
#'   items, and a per-item table.
#' @seealso [pord.test()], [pord.items()]
#' @examples
#' set.seed(1)
#' X <- matrix(sample(1:4, 60 * 5, TRUE, c(.05, .15, .35, .45)), 60, 5)
#' Y <- pmax(pmin(X - matrix(rbinom(60 * 5, 1, .4), 60, 5), 4), 1)
#' g <- rep(c("a", "b", "c"), each = 20)
#' pord.group(X, Y, g, B = 500)
#' @export
pord.group <- function(X, Y, group, direction = c("ge", "lt"), K = NULL,
                       B = 10000, p.adjust.method = "holm") {
  direction <- match.arg(direction)
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (!identical(dim(X), dim(Y))) stop("'X' and 'Y' must have the same shape")
  if (length(group) != nrow(X)) stop("'group' must have one entry per row")
  if (is.null(K)) K <- max(c(X, Y), na.rm = TRUE)
  OK <- !is.na(X) & !is.na(Y)
  HIT <- if (direction == "ge") (Y >= X) & OK else (Y < X) & OK
  g0 <- as.integer(factor(group)); lev <- levels(factor(group)); G <- length(lev)
  if (G < 2) stop("at least two groups are required")
  J <- ncol(X)

  prop <- function(g) vapply(seq_len(G), function(j) {
    i <- which(g == j); sum(HIT[i, ]) / sum(OK[i, ]) }, numeric(1))
  exc <- function(g) vapply(seq_len(G), function(j) {
    i <- which(g == j)
    xs <- X[i, ][OK[i, ]]; ys <- Y[i, ][OK[i, ]]
    e <- .pord.excess(xs, ys, K)
    if (direction == "ge") e else -e }, numeric(1))
  cells <- function(g) {
    nn <- xx <- matrix(0, J, G)
    for (j in seq_len(G)) { i <- which(g == j)
      nn[, j] <- colSums(OK[i, , drop = FALSE])
      xx[, j] <- colSums(HIT[i, , drop = FALSE]) }
    list(n = nn, x = xx)
  }

  o.prop <- diff(range(prop(g0)))
  o.exc  <- diff(range(exc(g0)))
  c0 <- cells(g0)
  o.chi <- vapply(seq_len(J), function(q) .pord.chi2(c0$n[q, ], c0$x[q, ]), numeric(1))

  sim <- replicate(B, {
    g <- sample(g0); cl <- cells(g)
    c(diff(range(prop(g))), diff(range(exc(g))),
      vapply(seq_len(J), function(q) .pord.chi2(cl$n[q, ], cl$x[q, ]), numeric(1)))
  })
  p.prop <- (sum(sim[1, ] >= o.prop) + 1) / (B + 1)
  p.exc  <- (sum(sim[2, ] >= o.exc) + 1) / (B + 1)
  p.item <- vapply(seq_len(J), function(q)
    (sum(sim[2 + q, ] >= o.chi[q]) + 1) / (B + 1), numeric(1))
  p.max  <- (sum(apply(sim[3:(2 + J), , drop = FALSE], 2, max) >= max(o.chi)) + 1) / (B + 1)

  pm <- c0$x / c0$n; colnames(pm) <- lev
  items <- data.frame(item = if (!is.null(colnames(X))) colnames(X) else
                        paste0("item", seq_len(J)),
                      pm, range = apply(pm, 1, function(v) diff(range(v))),
                      chisq = o.chi, p.value = p.item,
                      p.adjusted = stats::p.adjust(p.item, p.adjust.method),
                      check.names = FALSE, stringsAsFactors = FALSE)
  items <- items[order(-items$chisq), ]

  structure(list(
    groups = data.frame(group = lev, proportion = prop(g0), excess = exc(g0),
                        stringsAsFactors = FALSE),
    direction = direction,
    p.proportion = p.prop, p.excess = p.exc, p.max = p.max,
    items = items, B = B, p.adjust.method = p.adjust.method),
    class = "pord.group")
}

#' @export
print.pord.group <- function(x, ...) {
  lab <- if (x$direction == "ge") "y >= x" else "y < x"
  cat(sprintf("\n\tGroup comparison of paired ordinal pattern (%s)\n\n", lab))
  print(x$groups, row.names = FALSE, digits = 3)
  cat(sprintf("\npermutation of group labels, B = %s\n",
              format(x$B, big.mark = ",")))
  cat(sprintf("  strength of association (excess) : p = %.4f\n", x$p.excess))
  cat(sprintf("  raw proportion                   : p = %.4f\n", x$p.proportion))
  cat(sprintf("  item-wise, max statistic         : p = %.4f\n", x$p.max))
  cat(sprintf("  item-wise, %s-adjusted: %d of %d significant\n\n",
              x$p.adjust.method, sum(x$items$p.adjusted < 0.05), nrow(x$items)))
  print(utils::head(x$items, 5), row.names = FALSE, digits = 3)
  cat("\n")
  invisible(x)
}

#' @title Apply the paired ordinal tests across many items
#' @description
#' \code{pord.items} runs [pord.test()] (and optionally [pord.sign()]) on every
#' column of a pair of matrices, and adjusts for multiplicity across items.
#'
#' @param X Matrix or data frame of \code{x} ratings, **respondents in rows and
#'   items in columns**.
#' @param Y Matrix or data frame of \code{y} ratings, same shape as \code{X}.
#' @param K Number of scale points.  Taken from the data when \code{NULL}.
#' @param sign Also run the sign test for \eqn{y < x} on each item.
#' @param p.adjust.method Method passed to [stats::p.adjust()].
#' @param ... Passed to [pord.test()], e.g. \code{method} or \code{B}.
#'
#' @return A data frame with one row per item: sample size, observed
#'   proportion, the value expected under independence, the excess, the p-value
#'   and its adjusted version, and -- when \code{sign = TRUE} -- the
#'   corresponding sign-test columns.
#' @seealso [pord.test()], [pord.sign()], [pord.group()]
#' @examples
#' set.seed(1)
#' X <- matrix(sample(1:4, 60 * 4, TRUE, c(.05, .15, .35, .45)), 60, 4)
#' Y <- pmax(pmin(X - matrix(rbinom(60 * 4, 1, .4), 60, 4), 4), 1)
#' pord.items(X, Y, sign = TRUE)
#' @export
pord.items <- function(X, Y, K = NULL, sign = FALSE,
                       p.adjust.method = "holm", ...) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (!identical(dim(X), dim(Y))) stop("'X' and 'Y' must have the same shape")
  if (is.null(K)) K <- max(c(X, Y), na.rm = TRUE)
  nm <- if (!is.null(colnames(X))) colnames(X) else paste0("item", seq_len(ncol(X)))
  out <- do.call(rbind, lapply(seq_len(ncol(X)), function(j) {
    tt <- pord.test(X[, j], Y[, j], K = K, ...)
    row <- data.frame(item = nm[j], n = tt$n, proportion = tt$proportion,
                      expected = tt$theta, excess = tt$excess,
                      z = tt$z, p.value = tt$p.value, method = tt$method,
                      stringsAsFactors = FALSE)
    if (sign) {
      ss <- pord.sign(X[, j], Y[, j], K = K)
      row$prop.ceiling <- ss$prop.ceiling
      row$sign.n <- ss$n
      row$sign.below <- ss$below; row$sign.tied <- ss$tied; row$sign.above <- ss$above
      row$sign.p <- ss$p.value
    }
    row
  }))
  out$p.adjusted <- stats::p.adjust(out$p.value, p.adjust.method)
  if (sign) out$sign.p.adjusted <- stats::p.adjust(out$sign.p, p.adjust.method)
  out[order(-out$excess), ]
}
