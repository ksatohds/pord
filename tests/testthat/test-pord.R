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
test_that("the closed-form moments equal the moments of the exact null", {
  for (sd in 1:4) {
    d <- make_pair(seed = sd)
    mo <- pord.moments(d$x, d$y, K = d$K)
    tt <- pord.test(d$x, d$y, K = d$K, method = "exact")
    v  <- seq_along(tt$null.pmf) - 1
    mu <- sum(v * tt$null.pmf)
    s  <- sqrt(sum((v - mu)^2 * tt$null.pmf))
    expect_equal(mo$mean, mu, tolerance = 1e-9)
    expect_equal(mo$sd,   s,  tolerance = 1e-9)
  }
})

test_that("the exact null distribution is a proper distribution", {
  d <- make_pair(seed = 7)
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
  d <- make_pair(n = 150, seed = 3)
  e <- pord.test(d$x, d$y, K = d$K, method = "exact")$p.value
  set.seed(11)
  p <- pord.test(d$x, d$y, K = d$K, method = "permutation", B = 20000)$p.value
  n <- pord.test(d$x, d$y, K = d$K, method = "normal")$p.value
  expect_equal(p, e, tolerance = 0.02)          # Monte Carlo error
  expect_gt(n / e, 0.5); expect_lt(n / e, 2)    # continuity-corrected normal
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
test_that("the sign test drops ties and the conditional variant drops both ends", {
  d <- make_pair(seed = 6)
  s1 <- pord.sign(d$x, d$y, K = d$K)
  expect_equal(s1$below + s1$tied + s1$above, length(d$x))
  s2 <- pord.sign(d$x, d$y, K = d$K, conditional = TRUE)
  expect_equal(s2$n, sum(d$x > 1 & d$x < d$K))
  expect_lte(s2$n, s1$n)
  # the conditional variant is less biased when the ceiling is heavy
  set.seed(21)
  x <- rep(4L, 100); y <- sample(1:4, 100, TRUE)
  expect_lt(pord.sign(x, y, K = 4)$p.value, 1e-10)   # naive: rejects trivially
  expect_error(pord.sign(x, y, K = 4, conditional = TRUE))  # nothing left
})

## ---------------------------------------------------------------------
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
  out <- pord.items(X, Y, sign = TRUE, conditional = TRUE)
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
