# Accuracy

This document records how closely the polyfill reproduces two references. The
first is the R, SciPy, and mpmath values that the php_stats extension is
validated against. The second is the extension itself. Every number below was
measured, not estimated. A function is shipped only when it passes the
extension's own fixture tolerances and the parity suite. Anything that cannot
meet them stays undefined and is listed under [Not provided](#not-provided).

Reference: php-stats commit pinned in `tests/fixtures/REFERENCE`. The fixtures
are copied from that commit by `bin/fetch-fixtures` and never edited.

## How the numbers are produced

- `bin/accuracy-report` runs the polyfill on every case in the R fixtures and
  prints the tables in [Against R](#against-r). Each comparison uses the
  tolerance of the extension's own suite for that fixture, |a − b| ≤ abs + rel·|b|.
  "Tolerance used" is the largest |a − b| / (abs + rel·|b|) over all entries, so
  100% would be a failure. "Max relative error" is |a − b| / |b| over entries
  that are not rounding noise. It skips entries whose expected magnitude is
  below the absolute tolerance or below 10⁻³ of the largest value in the same
  result entry, such as the residuals of an exact fit.
- The parity suites in `tests/Parity` compare the polyfill with the loaded
  extension on the fixtures and on seeded random and adversarial inputs. Run
  them with `PARITY_REPORT=parity.json` to record the largest differences.
- Both are part of the test suite: `composer test` (polyfill only) and
  `composer test:ext` (with the extension loaded).

## Descriptive statistics

The descriptive functions are ports of the extension's algorithms, including
its two-pass and compensated summations. Against the extension they agree
exactly (largest difference 0.0) on the fixtures and on seeded random inputs.
The random inputs include values offset by 10⁹ and 10¹² (catastrophic
cancellation for a naive variance), integers, tiny values, n from 1 to 500,
and null/NaN under both missing-value policies. Exception codes are identical
for every invalid input tested.

## Regression models

`stats_ols` solves least squares by Householder QR of the column-scaled design,
then a one-sided Jacobi SVD of the small R factor. The singular values are
those of the design, so the extension's rank rule (s ≤ ε·max(n, p)·s_max)
makes the same decisions. `stats_logistic`, `stats_probit`, and
`stats_poisson` port the extension's damped Newton / Fisher-scoring loop step
for step. This includes the starting values, line search, separation
certificate, convergence rule, compensated likelihood sums, and statistics. The
statrs special functions behind p-values and likelihood constants (erfc, ln Γ,
the regularized incomplete beta, and the inverse error function) are ported
from the same crate version.

### Against R

| Function | Reference | Comparisons | Max relative error | Tolerance used | Tested domain |
|---|---|---:|---:|---:|---|
| `stats_cor` | R (descriptive-r.json) | 3 | 1.3e-16 | 0% | correlated, negative, offset |
| `stats_cov` | R (descriptive-r.json) | 6 | 2.1e-16 | 0% | correlated, negative, offset |
| `stats_logistic` | R (logistic-r.json) | 78 | 7.4e-9 | 6% | n=80, p=3, n=80, p=2, n=80, p=1 |
| `stats_logistic` | R (r-compatibility.json) | 52 | 3.5e-6 | 20% | logistic_intercept_individual, logistic_intercept_grouped, logistic_no_intercept_individual, logistic_no_intercept_grouped |
| `stats_logistic` | R (stress-r.json) | 44 | 2.5e-12 | 0% | logistic_rare, logistic_offset_-1000, logistic_offset_1000, logistic_extreme_known_predictors |
| `stats_mean` | R (descriptive-r.json) | 8 | 1.5e-16 | 0% | small, integers, offset, tiny, and 4 more |
| `stats_median` | R (descriptive-r.json) | 8 | 0.0e+0 | 0% | small, integers, offset, tiny, and 4 more |
| `stats_ols` | R (ols-r.json) | 108 | 3.5e-13 | 0% | n=8, p=2, n=8, p=3, n=5, p=2, n=10, p=3, and 1 more |
| `stats_ols` | R (r-compatibility.json) | 30 | 3.8e-12 | 0% | ols_intercept_individual, ols_no_intercept_individual |
| `stats_ols` | R (stress-r.json) | 54 | 1.6e-8 | 0% | nearly_collinear, poorly_scaled, small_outcome, large_outcome, and 2 more |
| `stats_poisson` | R (poisson-r.json) | 78 | 8.8e-11 | 0% | n=120, p=3, n=120, p=2, n=120, p=1 |
| `stats_poisson` | R (r-compatibility.json) | 26 | 8.0e-7 | 39% | poisson_intercept_individual, poisson_no_intercept_individual |
| `stats_poisson` | R (stress-r.json) | 55 | 1.7e-8 | 21% | poisson_tiny_known_means, poisson_sparse, poisson_large, poisson_offset_-1000, and 1 more |
| `stats_probit` | R (probit-r.json) | 78 | 6.3e-8 | 17% | n=120, p=3, n=120, p=2, n=120, p=1 |
| `stats_probit` | R (r-compatibility.json) | 52 | 1.8e-5 | 20% | probit_intercept_individual, probit_intercept_grouped, probit_no_intercept_individual, probit_no_intercept_grouped |
| `stats_probit` | R (stress-r.json) | 44 | 7.4e-11 | 0% | probit_rare, probit_offset_-1000, probit_offset_1000, probit_extreme_known_predictors |
| `stats_quantile_array` | R (descriptive-r.json) | 8 | 0.0e+0 | 0% | small, integers, offset, tiny, and 4 more |
| `stats_quantile` | R (descriptive-r.json) | 80 | 0.0e+0 | 0% | small, integers, offset, tiny, and 4 more |
| `stats_sd` | R (descriptive-r.json) | 8 | 3.9e-16 | 0% | small, integers, offset, tiny, and 4 more |
| `stats_var` | R (descriptive-r.json) | 16 | 5.2e-16 | 0% | small, integers, offset, tiny, and 4 more |
| `stats_weighted_mean` | R (descriptive-r.json) | 2 | 2.1e-16 | 0% | analytic, counts |
| `stats_weighted_var` | R (descriptive-r.json) | 6 | 1.9e-16 | 0% | analytic, counts |

The GLM differences from R (up to about 2·10⁻⁵ relative in the probit
r-compatibility case) come from R's own stopping rule (`glm` stops when the
relative deviance change is below 10⁻⁸). The extension shows the same
differences. Descriptive statistics are compared at a relative tolerance of
10⁻¹².

### Against the extension

The parity suite compares every result key. Each value is judged on the scale
of its quantity. Coefficients and intervals are compared relative to
|b| + standard error, t/z statistics and R² relative to |b| + 1, p-values on the
log scale, and fitted values and residuals relative to |b| + max |y|. Other keys
are compared relative to |b| + 10⁻³·max |b|. The required agreement is 10⁻¹⁰
for OLS on the fixtures, 10⁻⁹ for OLS on random data, and 10⁻⁷ for the
iterative models. The iterative models stop once a full step is below 10⁻⁸
relative, so two implementations agree to about that level and not to
rounding.

| Group | Inputs | OLS | Logistic | Probit | Poisson |
|---|---|---:|---:|---:|---:|
| fixtures | all R fixtures, with and without the automatic intercept | 2.7e-10 | 2.2e-12 | 1.1e-12 | 2.3e-13 |
| random | 40 seeded designs per model: n 12–400, 1–5 predictors, weights, offsets, missing rows | 7.7e-13 | 6.1e-12 | 7.9e-9 | 6.8e-14 |
| scaling | predictors and outcomes scaled by 10⁻⁸ to 10⁸ | 1.4e-13 | 2.3e-13 | 3.5e-14 | 1.8e-14 |
| collinear | condition numbers about 10⁶ to 10¹⁰ | 1.2e-5 | 8.8e-9 | (rejected identically) | 3.1e-9 |
| separation | complete, quasi-complete, rare events (1 in 2000), grouped near-boundary proportions | | 9.7e-9 | 2.5e-14 | |
| offsets / counts | binary offsets −50 to 40; Poisson counts up to 10¹², offsets −700 to 300, sparse counts | | 5.4e-14 | 7.5e-14 | 1.2e-8 |
| errors | invalid shapes, types, domains, weights, policies, rank deficiency | 1.3e-13 | same codes | same codes | 6.7e-8 |

Every invalid input in the error group raises the same exception code in both.
In the collinear OLS group, two backward-stable solvers agree to about
cond(X)·ε. At a condition number of 10¹⁰, 1.2·10⁻⁵ of a standard error is the
expected agreement. Probit fits on the collinear designs are rejected by both
implementations with the same code.

The stopping iteration can differ between the two. On the fixtures and random
data it differs by at most one. In the stress groups the likelihood can be flat
to rounding near the optimum (for example, Poisson offsets near 300 make the
linear predictor a difference of numbers near 300). There the stopping
iteration differed by up to two while the estimates still agreed to 10⁻⁸.

## Distributions

All 13 families in the extension are provided: normal, t, chi-square, F,
exponential, logistic, log-normal, gamma, beta, binomial, Poisson, negative
binomial, and negative binomial by mean. Each family has 14 functions. These
are the density or mass function and its logarithm, cdf, logcdf, sf, logsf,
and quantile, each with an `_array` variant. That makes 182 functions in all.

The implementation ports the extension's distribution layer and the parts of
the cdflib crate it calls. The cdflib parts are TOMS 708 (incomplete beta),
TOMS 654 (incomplete gamma), Cody's normal CDF, and the `dinvnr`/`dt1`
quantile starts. The extension's log-scale tails, Loader's saddle-point
densities, and its safeguarded Newton quantile search are ported as well.

### Against the references

The fixture comparisons use the extension's own tolerances. These are
1e-12 relative for probabilities and densities (with a 1e-300 floor),
1e-12·max(1, |value|) for log-scale values, and 1e-15 + 1e-12·|value| for
quantiles. The table gives the largest error observed against each reference,
measured the same way. R and SciPy are measured only at points where that
reference itself agrees with mpmath, and the bracketed count is the number of
points excluded because it did not.

| Family | mpmath | R | SciPy | Tested domain |
|---|---:|---:|---:|---|
| normal | 2.9e-14 | 1.8e-14 | 1.1e-13 | sd 1e-8 to 1e3, \|x\| to 1e13 |
| t | 1.8e-13 | 2.6e-13 [4] | 6.6e-13 [40] | df 0.5 to 1e6, \|x\| to 1e10 |
| chi-square | 8.3e-13 | 8.3e-13 [3] | 8.9e-13 [64] | df 0.5 to 1e6, x to 1e7 |
| F | 1.6e-13 | 5.2e-13 [15] | 3.4e-13 [46] | df1 0.5 to 100, df2 0.5 to 1e6, x to 1e20 |
| exponential | 1.5e-14 | 2.1e-16 | 2.1e-16 | rate 1e-5 to 1e5 |
| logistic | 1.5e-15 | 2.0e-16 | 2.0e-16 | scale 1e-5 to 2 |
| log-normal | 9.1e-14 | 2.1e-14 [1] | 1.1e-13 | sdlog 0.5 to 3 |
| gamma | 7.5e-13 | 7.5e-13 | 5.6e-13 [40] | shape 0.01 to 1e4, scale 0.01 to 1e5 |
| beta | 1.5e-13 | 5.7e-14 | 2.6e-13 [15] | shapes 0.1 to 1e3 |
| binomial | 3.3e-14 | 5.7e-14 [3] | 9.9e-13 [31] | size up to 1e6 |
| Poisson | 3.7e-13 | 3.7e-13 [6] | 3.7e-13 [57] | lambda 0 to 1e6 |
| negative binomial | 6.8e-14 | 2.5e-14 [4] | 4.5e-14 [35] | size 0.5 to 1e5 |
| negative binomial (mu) | 5.1e-14 | 3.8e-13 [10] | 9.5e-13 [70] | mu 0 to 1e5, size 0.5 to 1e6 |

Probabilities are tested down to about 1e-300, and the log-scale functions
beyond the double range. `tests/Unit/DistributionsStressTest.php` adds
comparisons with closed forms and summed series in the extreme tails. These
cover the normal asymptotics to z = 1e150, Cauchy and t with df = 2 to 1e300,
chi-square with df 1 and 2, the incomplete gamma for integer shapes up to 1e6,
beta with a unit shape, F with df 2, and the support edges of the binomial and
negative binomial. The suite also checks densities at extreme parameters and
quantile/cdf round trips for p from 1e-300 to 1 − 2⁻⁵².

### Against the extension

Every family is bit-identical to the extension: the largest difference is 0 on
about 8,600 fixture points and on seeded sweeps over tails, ±INF evaluation
points, p in {0, 1, 5e-324, 1e-300, 1 − 2⁻⁵²}, and extreme parameters.
Exception codes are the same for every invalid input tested.

### Limits shared with the extension

For gamma, chi-square, and Poisson, log tails below 1e-280 lose accuracy at
shape or lambda of about 1e6 and above. The relative error is 2e-12 at 1e6,
7e-12 at 1e7, and 5e-10 at 1e8, from cancellation in the extension's
formula. The polyfill reproduces the extension's values exactly there, and a
stress test records the bound.

### Extreme parameters

cdflib's incomplete beta and gamma routines never return for some extreme
arguments and panic in the deep asymptotic regime of the incomplete gamma.
Since php-stats bd98bf7 the extension returns a failure up front for those
regimes, and the polyfill applies the same conditions. They are an overflowing
a + b, a beta parameter above 1e300 with x or 1 − x below 1e-300, and a gamma
shape or point above 1e300. A panic inside either routine also becomes a
failure. The functions report every such failure as exception 1006, and a
quantile search stops with 1006 as soon as a tail evaluation fails.
`tests/Parity/DistributionsExtremeParityTest.php` compares both
implementations on a grid of parameters from 1e-300 to `PHP_FLOAT_MAX` for F,
beta, gamma, chi-square, t, Poisson, and the negative binomial by mean. The
values or exception codes agree on every point.

For `stats_t_quantile` at very small df, the quantile search used to step from
a very large starting point by ±2 and could take seconds per call. Both the
extension (since 0aafc26) and the polyfill now search that run of points by
exponential probing and bisection. The visited points are the same, so the
results equal those of the original walk.

## Formula interface

`stats_model_matrix`, `stats_lm`, and `stats_glm` port the extension's formula
parser, term expansion, evaluation, and coding. The fits reuse the polyfill's
OLS and GLM code described above.

### Against R

The 62 cases in `formula-r.json` cover 33 model matrices, 14 `lm` fits, and 15
`glm` fits. The tolerances are those of the extension's suite: 1e-12 for
design matrices (1e-9 + 1e-7·|b| for the orthogonal `poly()` case) and
1e-9 + 1e-7·|b| for fits. Column names, terms, and retained rows must match
exactly.

| Function | Comparisons | Max relative error | Tolerance used | Datasets |
|---|---:|---:|---:|---|
| `stats_model_matrix` | 266 | 6.7e-16 | 0% | housing, schools, cbpp, synthetic, small |
| `stats_lm` | 238 | 2.7e-12 | 0% | housing, schools, synthetic |
| `stats_glm` | 240 | 5.9e-10 | 1% | learning, cbpp, clinic, synthetic |

### Against the extension

`tests/Parity/FormulaParityTest.php` compares the two on every fixture case.
It also covers about 170 syntax forms, about 120 malformed formulas, factor
level and collation edge cases, both missing-value policies, about 28 weight
shapes, 20 family and link combinations, about 50 data shapes, and 30 random
designs. Design matrices, column names, terms, `assign`, row accounting, key
order, and warning text are identical, and every invalid input raises the same
exception code with the same message. The largest scale-aware difference in
the fits is 4.4e-13 on the fixtures and 7.9e-13 across the other cases.

## Not provided

The following functions need the native extension and are left undefined:
`stats_lmm`, `stats_lmm_fit`, `stats_logistic_mixed`, `stats_lmer`,
`stats_glmer`, every `*_rng` function,
the `stats_draw_*` functions, and `stats_hello`.
`PhpStats\Polyfill\Polyfill::supports()` reports which functions are
available.

## Known differences from the extension

- The listwise-deletion warning is raised with `E_USER_WARNING`, because PHP
  code cannot raise `E_WARNING`. The message text is identical, including the
  `stats_ols(): ` prefix.
- A scalar argument of the wrong type follows PHP's rules for a userland
  function. In files without `strict_types`, numeric strings and booleans are
  coerced, and other values raise `TypeError`. The extension raises an
  `\Exception` with the message "Invalid value given for argument" instead.
- The reported `iterations` of an iterative fit can differ, as described above.
