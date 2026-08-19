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
