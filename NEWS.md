# pord 0.2.1 (2026-08-19)

## Wording and test coverage

No change to any computed value.

- `pord.sign()` was described as "the exact conditional test of symmetry".
  It is not: the null it tests is directional balance, `P(Y < X) = P(Y > X)`.
  Symmetry implies that null but is not implied by it, because the test pools
  every off-diagonal cell into two totals.  With `n_12 = 30, n_21 = 0` and
  `n_34 = 0, n_43 = 30` the test returns `p = 1` while Bowker's test of
  symmetry gives `p < 1e-13`.  The help page now says so and points to a
  Bowker or exact conditional symmetry test when symmetry is the question.
- `pord.test()` and `pord.moments()` were titled in terms of "paired ordinal
  superiority", which reads as stochastic dominance or a comparison of
  margins.  Both margins are conditioned on, so a marginal shift is not what
  the test can see; what varies is the pairing.  Retitled around the quantity
  actually used: an excess of `y >= x` pairs over independence.
- The level test for `pord.sign()` asserted only that the rejection rate was
  below 0.07, which a test that never rejects would also pass -- precisely
  the failure mode of the `conditional` variant removed in 0.2.0.  It now
  brackets the rate from both sides, on a bottom-heavy scale as well as a
  top-heavy one, and a companion test checks power against a real downward
  shift.

# pord 0.2.0 (2026-08-19)

## Breaking change: `conditional` removed from `pord.sign()` and `pord.items()`

The `conditional = TRUE` variant (restricting the sign test to respondents
with `1 < x < K`) reported incorrect p-values and has been removed.

The restriction is asymmetric in `(x, y)`, so under the null of symmetry the
probability that a discordant pair has `y < x` is not one half.  For example,
with `x`, `y` iid on `P(1,2,3,4) = (.05,.15,.35,.45)` -- a fully symmetric
joint -- that probability is 0.218 among respondents with `x` in `{2, 3}`;
comparing against 0.5 makes the test conservative there and anticonservative
when the scale is bottom-heavy.

The motivation for the variant was also mistaken: a ceiling in `x` does not
bias the ordinary sign test.  Symmetry implies equal margins, so `x`
concentrating at the top of the scale while `y` sits below it is exactly the
asymmetry the test is built to detect.  Use the ordinary `pord.sign()`;
to rank items or domains by shortfall, use `pord.compare()`.

# pord 0.1.0 (2026-08-16)

First public release.
