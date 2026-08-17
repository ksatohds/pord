# pord

Exact conditional tests for paired ordinal data on a square table.

Two ratings on the same K-point scale are collected from the same respondent —
for example how much is *expected* of a task and how much is *achieved*. The
question is whether the second reaches the first more often than independence
would give, and, separately, which items fall short.

## The test

Let `T = #{i : y_i >= x_i}`, the sum of the diagonal and upper triangle of the
K × K table with rows `x` and columns `y`. `pord.test()` conditions on **both
margins** and computes the exact null distribution of `T`.

This is a generalization of Fisher's exact test. It uses the same conditioning
— the Fisher–Yates multiple hypergeometric distribution over all tables with the
observed margins — and **reduces exactly to Fisher's exact test when the scale is
collapsed to two levels**. It differs only in the statistic: the
Fisher–Freeman–Halton extension of Fisher's test to K × K tables uses the
probability of the table itself and detects departures in every direction,
whereas `pord.test()` counts one direction only and so answers a narrower,
more demanding question.

The null distribution is obtained **exactly**, by dynamic programming over all
tables with the given margins — no simulation. The null mean and variance have
closed forms depending only on the two margins, so they are available for any K
at no cost. Because both depend on the margins, **the null distribution is not
shared between items**.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ksatohds/pord")
```

## Usage

```r
library(pord)

set.seed(1)
x <- sample(1:4, 96, TRUE, c(.05, .15, .35, .45))   # expectation
y <- pmax(1, pmin(4, x - rbinom(96, 1, .35)))       # achievement

pord.test(x, y)            # does achievement reach expectation?
pord.sign(x, y, conditional = TRUE)   # where does it fall short?
pord.table(x, y)           # the 4 x 4 table split into three regions
```

For many items at once, and for comparing groups of respondents:

```r
pord.items(X, Y, sign = TRUE)              # one row per item, Holm-adjusted
pord.group(X, Y, group = occupation)       # does the pattern differ by group?
pord.compare(X, Y, domain = subscale)      # which domains fall short more?
```

`X` and `Y` have **respondents in rows and items in columns**. `group` is a
person-level variable, one entry per row, following the convention of
`mirt::multipleGroup()` and `lavaan::cfa()`; `domain` groups the *columns*.
Note that this is not the terminology of `stats::friedman.test()`, where
`groups` refers to the columns and `blocks` to the rows.

## Two questions, two nulls

`pord.test()` has **independence** as its null. Its lower tail does *not* ask
whether `y` falls short of `x`: under that null the shortfall count is `n - T`,
so a large value simply means negative association.

The natural null for the shortfall question is **symmetry**,
`P(y < x) = P(y > x)`. Conditioning on the off-diagonal pair totals removes the
diagonal (ties drop out as a consequence of the conditioning, not by choice) and
leaves a binomial — the sign test, provided by `pord.sign()`.

When `x` sits at an end of the scale one direction is impossible, so the
unconditional sign test is biased. `conditional = TRUE` restricts to respondents
with `1 < x < K`, for whom both directions remain available.

## Choosing a method

| method | cost | accuracy |
|---|---|---|
| `"exact"` | seconds for small K | exact |
| `"permutation"` | `B` shuffles | exact in principle; resolution limited to `1/(B+1)` |
| `"normal"` | instant | closed-form moments; **needs the continuity correction**, which is on by default |

`method = "auto"` counts the tables with the given margins and picks the exact
computation when that is affordable.

The statistic is integer-valued and its lattice spacing is a substantial
fraction of the null standard deviation, so an uncorrected normal tail
understates the p-value by roughly a factor of 1.5. The correction is applied
unless `continuity = FALSE`.

## Notes

Power depends on the margins. The excess over independence collapses to zero
when `x` takes a single value, whatever that value is, so items whose first
rating is concentrated carry little information regardless of the second.
`pord.sign()` reports the proportion at the ceiling for this reason.

The statistic is **not** monotone in the positive-quadrant-dependence ordering:
a margin-preserving transfer that increases concordance can decrease
`P(y >= x)`. The test therefore targets a specific direction and should not be
described as a test for positive dependence.

## References

- Fisher, R. A. (1935). The logic of inductive inference. *JRSS* **98**, 39–82.
- Wald, A. & Wolfowitz, J. (1944). Statistical tests based on permutations of the observations. *Ann. Math. Statist.* **15**, 358–372.
- Hoeffding, W. (1951). A combinatorial central limit theorem. *Ann. Math. Statist.* **22**, 558–566.
- Strasser, H. & Weber, C. (1999). On the asymptotic theory of permutation statistics. *Math. Methods Statist.* **8**, 220–250.
- Svensson, E. (2001). Guidelines to statistical evaluation of data from rating scales and questionnaires. *J. Rehabil. Med.* **33**, 47–48.
- Phipson, B. & Smyth, G. K. (2010). Permutation p-values should never be zero. *Stat. Appl. Genet. Mol. Biol.* **9**, Article 39.

## License

MIT © Kenichi Satoh
