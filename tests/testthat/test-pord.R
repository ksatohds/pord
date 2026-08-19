# Tests use synthetic data only.  The empirical study that motivated the
# package is not distributed; its reproduction script lives in
# inst/reproduce/ and is excluded from the build.

make_pair <- function(n = 120, K = 4, seed = 1, drop = 0.35) {
  set.seed(seed)
  x <- sample(seq_len(K), n, TRUE, prob = seq_len(K))
  y <- pmax(1L, pmin(K, x - stats::rbinom(n, 1, drop)))
  list(x = x, y = y, K = K)
}

## ---------------------------------------------------------------------
## The exact null is a dynamic program over all tables with the observed
## margins, and its cost climbs steeply in n: on one machine, n = 40 took
## 0.08s, n = 80 took 7s and n = 120 took 92s.  The identities being checked
## hold at every n, so the routine tests use a small n and run everywhere;
## the large-n versions carry skip_on_cran(), which keeps CRAN's check short
## without leaving these properties unverified there.
N_FAST <- 40
N_SLOW <- 120

test_that("the closed-form moments equal the moments of the exact null", {
  for (sd in 1:4) {
    d <- make_pair(n = N_FAST, seed = sd)
    mo <- pord.moments(d$x, d$y, K = d$K)
    tt <- pord.test(d$x, d$y, K = d$K, method = "exact")
    v  <- seq_along(tt$null.pmf) - 1
    mu <- sum(v * tt$null.pmf)
    s  <- sqrt(sum((v - mu)^2 * tt$null.pmf))
    expect_equal(mo$mean, mu, tolerance = 1e-9)
    expect_equal(mo$sd,   s,  tolerance = 1e-9)
  }
})

test_that("the closed-form moments still match at a larger sample", {
  skip_on_cran()
  d <- make_pair(n = N_SLOW, seed = 1)
  mo <- pord.moments(d$x, d$y, K = d$K)
  tt <- pord.test(d$x, d$y, K = d$K, method = "exact")
  v  <- seq_along(tt$null.pmf) - 1
  mu <- sum(v * tt$null.pmf)
  expect_equal(mo$mean, mu, tolerance = 1e-9)
  expect_equal(mo$sd, sqrt(sum((v - mu)^2 * tt$null.pmf)), tolerance = 1e-9)
})

test_that("the exact null distribution is a proper distribution", {
  d <- make_pair(n = N_FAST, seed = 7)
  tt <- pord.test(d$x, d$y, K = d$K, method = "exact")
  expect_equal(sum(tt$null.pmf), 1, tolerance = 1e-10)
  expect_true(all(tt$null.pmf >= -1e-12))
})

test_that("the exact null is still a proper distribution at a larger sample", {
  skip_on_cran()
  d <- make_pair(n = N_SLOW, seed = 7)
  tt <- pord.test(d$x, d$y, K = d$K, method = "exact")
  expect_equal(sum(tt$null.pmf), 1, tolerance = 1e-10)
  expect_true(all(tt$null.pmf >= -1e-12))
})

## ---------------------------------------------------------------------
test_that("with two scale points the test reduces to Fisher's exact test", {
  for (sd in 1:5) {
    set.seed(sd)
    x <- sample(1:2, 80, TRUE); y <- sample(1:2, 80, TRUE)
    got <- pord.test(x, y, K = 2, method = "exact", alternative = "less")$p.value
    ref <- stats::fisher.test(table(factor(x, 1:2), factor(y, 1:2)),
                              alternative = "less")$p.value
    expect_equal(got, ref, tolerance = 1e-10)
  }
})

## ---------------------------------------------------------------------
test_that("exact, permutation and normal agree well enough", {
  d <- make_pair(n = N_FAST, seed = 3)
  e <- pord.test(d$x, d$y, K = d$K, method = "exact")$p.value
  set.seed(11)
  p <- pord.test(d$x, d$y, K = d$K, method = "permutation", B = 5000)$p.value
  n <- pord.test(d$x, d$y, K = d$K, method = "normal")$p.value
  expect_equal(p, e, tolerance = 0.03)          # Monte Carlo error, B = 5000
  expect_gt(n / e, 0.5); expect_lt(n / e, 2)    # continuity-corrected normal
})

test_that("the three methods still agree at the sample size a user would have", {
  skip_on_cran()
  d <- make_pair(n = 150, seed = 3)
  e <- pord.test(d$x, d$y, K = d$K, method = "exact")$p.value
  set.seed(11)
  p <- pord.test(d$x, d$y, K = d$K, method = "permutation", B = 20000)$p.value
  n <- pord.test(d$x, d$y, K = d$K, method = "normal")$p.value
  expect_equal(p, e, tolerance = 0.02)
  expect_gt(n / e, 0.5); expect_lt(n / e, 2)
})

test_that("the continuity correction matters", {
  d <- make_pair(n = 150, seed = 5)
  with.cc <- pord.test(d$x, d$y, K = d$K, method = "normal", continuity = TRUE)$p.value
  no.cc   <- pord.test(d$x, d$y, K = d$K, method = "normal", continuity = FALSE)$p.value
  expect_gt(with.cc, no.cc)                     # the correction is conservative
})

## ---------------------------------------------------------------------
test_that("standardizing does not change the permutation test", {
  d <- make_pair(seed = 9)
  set.seed(2); obs <- sum(d$y >= d$x)
  sim <- replicate(5000, sum(sample(d$y) >= d$x))
  m0 <- mean(sim); s0 <- stats::sd(sim)
  pT <- (sum(sim >= obs) + 1) / 5001
  pZ <- (sum((sim - m0) / s0 >= (obs - m0) / s0) + 1) / 5001
  expect_identical(pT, pZ)
})

## ---------------------------------------------------------------------
test_that("pord.table splits the square table into three regions", {
  d <- make_pair(seed = 4)
  tb <- pord.table(d$x, d$y, K = d$K)
  expect_equal(sum(tb$regions$observed), tb$n)
  expect_equal(sum(tb$regions$expected), tb$n, tolerance = 1e-10)
  expect_equal(tb$regions$observed[2] + tb$regions$observed[3],
               sum(d$y >= d$x))
  # a table can be passed instead of two vectors
  tb2 <- pord.table(tb$table)
  expect_equal(tb$regions$observed, tb2$regions$observed)
})

## ---------------------------------------------------------------------
test_that("the sign test drops ties and matches the classical binomial", {
  d <- make_pair(seed = 6)
  s1 <- pord.sign(d$x, d$y, K = d$K)
  expect_equal(s1$below + s1$tied + s1$above, length(d$x))
  ref <- stats::binom.test(s1$below, s1$below + s1$above, 0.5,
                           alternative = "greater")$p.value
  expect_equal(s1$p.value, ref, tolerance = 1e-12)
})

test_that("the sign test holds its level under symmetry, either way up", {
  # X, Y iid: the joint is symmetric, so the rejection rate must sit AT the
  # nominal level -- not merely below it.  The removed `conditional` variant
  # failed by being far too conservative here (its selection on 1 < x < K
  # made the discordant probability 0.218 rather than 0.5), which an
  # upper-bound-only assertion would have passed.  Both a top-heavy and a
  # bottom-heavy scale are checked, since the old error changed sign with the
  # shape of the margin.
  for (pr in list(top = c(.05, .15, .35, .45), bottom = c(.45, .35, .15, .05))) {
    set.seed(22)
    rej <- mean(replicate(2000, {
      x <- sample(1:4, 60, TRUE, pr); y <- sample(1:4, 60, TRUE, pr)
      pord.sign(x, y, K = 4, alternative = "two.sided")$p.value < 0.05
    }))
    expect_gt(rej, 0.02)   # not conservative: it must actually reject at 5%
    expect_lt(rej, 0.08)   # not anticonservative
  }
})

test_that("the sign test has power against a real downward shift", {
  # The mirror of the level test: when y genuinely falls short, it must reject.
  set.seed(23)
  rej <- mean(replicate(400, {
    x <- sample(1:4, 60, TRUE, c(.05, .15, .35, .45))
    y <- pmax(1L, x - stats::rbinom(60, 1, 0.35))
    pord.sign(x, y, K = 4)$p.value < 0.05
  }))
  expect_gt(rej, 0.95)
})

## ---------------------------------------------------------------------
test_that("Kendall's W shares the tie correction of the Friedman chi-squared", {
  # A corrected chi-squared paired with an uncorrected W is the bug fixed in
  # 0.2.3: a reader computing chi2 / (n (G - 1)) must recover the reported W.
  set.seed(31)
  X <- matrix(sample(1:4, 80 * 6, TRUE, prob = 1:4), 80, 6)
  Y <- pmax(pmin(X - matrix(stats::rbinom(80 * 6, 1, .4), 80, 6), 4), 1)
  cp <- pord.compare(X, Y, domain = rep(c("A", "B", "C"), each = 2))
  n <- cp$n; G <- nrow(cp$domains)
  expect_equal(cp$kendall.W,
               as.numeric(cp$friedman$statistic) / (n * (G - 1)),
               tolerance = 1e-12)
  # and it equals the tie-corrected closed form
  Z <- vapply(c("A", "B", "C"), function(b) {
    j <- rep(c("A", "B", "C"), each = 2) == b
    rowMeans((Y[, j, drop = FALSE] < X[, j, drop = FALSE]) * 1)
  }, numeric(nrow(X)))
  rk <- t(apply(Z, 1, rank))
  S <- sum((colSums(rk) - n * (G + 1) / 2)^2)
  TT <- sum(apply(rk, 1, function(r) { u <- as.numeric(table(r)); sum(u^3 - u) }))
  expect_equal(cp$kendall.W,
               12 * S / (n^2 * (G^3 - G) - n * TT), tolerance = 1e-10)
  expect_gte(cp$kendall.W, 0); expect_lte(cp$kendall.W, 1)
})

test_that("pord.compare returns one row per pair of domains", {
  set.seed(8)
  X <- matrix(sample(1:4, 60 * 6, TRUE, prob = 1:4), 60, 6)
  Y <- pmax(pmin(X - matrix(stats::rbinom(60 * 6, 1, .4), 60, 6), 4), 1)
  cp <- pord.compare(X, Y, domain = rep(c("A", "B", "C"), each = 2))
  expect_equal(nrow(cp$domains), 3L)
  expect_equal(nrow(cp$pairwise), 3L)            # choose(3, 2)
  expect_true(all(cp$pairwise$p.adjusted >= cp$pairwise$p.value))
})

## ---------------------------------------------------------------------
test_that("pord.group is invariant to relabelling the groups", {
  set.seed(12)
  X <- matrix(sample(1:4, 60 * 5, TRUE, prob = 1:4), 60, 5)
  Y <- pmax(pmin(X - matrix(stats::rbinom(60 * 5, 1, .4), 60, 5), 4), 1)
  g <- rep(c("a", "b", "c"), each = 20)
  set.seed(1); r1 <- pord.group(X, Y, g, B = 300)
  set.seed(1); r2 <- pord.group(X, Y, factor(g, labels = c("x", "y", "z")), B = 300)
  expect_equal(r1$p.excess, r2$p.excess)
  expect_equal(r1$groups$proportion, r2$groups$proportion)
})

test_that("the two directions of pord.group are complementary up to ties", {
  set.seed(13)
  X <- matrix(sample(1:4, 60 * 5, TRUE, prob = 1:4), 60, 5)
  Y <- pmax(pmin(X - matrix(stats::rbinom(60 * 5, 1, .4), 60, 5), 4), 1)
  g <- rep(c("a", "b"), each = 30)
  set.seed(1); a <- pord.group(X, Y, g, direction = "ge", B = 300)
  set.seed(1); b <- pord.group(X, Y, g, direction = "lt", B = 300)
  expect_equal(a$groups$proportion + b$groups$proportion,
               rep(1, 2), tolerance = 1e-12)
})

## ---------------------------------------------------------------------
test_that("pord.items adjusts across items and orders by excess", {
  set.seed(15)
  X <- matrix(sample(1:4, 60 * 4, TRUE, prob = 1:4), 60, 4)
  Y <- pmax(pmin(X - matrix(stats::rbinom(60 * 4, 1, .4), 60, 4), 4), 1)
  colnames(X) <- colnames(Y) <- paste0("Q", 1:4)
  out <- pord.items(X, Y, sign = TRUE)
  expect_equal(nrow(out), 4L)
  expect_true(all(out$p.adjusted >= out$p.value))
  expect_true(!is.unsorted(rev(out$excess)))
  expect_true(all(c("sign.below", "sign.p.adjusted") %in% names(out)))
})

## ---------------------------------------------------------------------
test_that("missing values and bad input are handled", {
  d <- make_pair(seed = 2)
  x <- d$x; y <- d$y; x[1:5] <- NA
  expect_equal(pord.test(x, y, K = 4)$n, length(d$x) - 5L)
  expect_error(pord.test(1:10, 1:9))
  expect_error(pord.test(c(0L, 1L), c(1L, 1L), K = 2))   # outside 1..K
  expect_error(pord.table(matrix(1:6, 2, 3)))            # not square
})
