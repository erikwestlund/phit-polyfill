# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project uses [Semantic Versioning](https://semver.org/).

Each release states the php-stats commit it reproduces. The same commit is
pinned in `tests/fixtures/REFERENCE`, and the reference fixtures are copied
from it.

## [0.1.0] - Unreleased

Reproduces php-stats commit
[`91ee538`](https://github.com/erikwestlund/php-stats/commit/91ee538253e3ecf7b3f086441305c3f7237200fc).
php-stats has no tagged release yet.

### Added

- Polyfill bootstrap that defines each `stats_*` function only when the
  php_stats extension does not, generated from the extension's API contract.
  Signatures, parameter names, and defaults match the extension.
- `PhpStats\Polyfill\Polyfill::supports()`, `implementation()`,
  `functions()`, and `nativeLoaded()`.
- Descriptive statistics: `stats_mean`, `stats_var`, `stats_sd`,
  `stats_median`, `stats_quantile`, `stats_quantile_array`, `stats_cov`,
  `stats_cor`, `stats_weighted_mean`, and `stats_weighted_var`.
- Regression: `stats_ols`, `stats_logistic`, `stats_probit`, and
  `stats_poisson`, with weights, offsets, the automatic intercept, listwise
  deletion, and the extension's result keys, diagnostics, and exception codes.
- Distributions: 13 families (normal, t, chi-square, F, exponential, logistic,
  log-normal, gamma, beta, binomial, Poisson, negative binomial, and negative
  binomial by mean). Each has the density or mass function, cdf, survival
  function, their logarithms, and the quantile, with `_array` variants (182
  functions).
- Formula interface: `stats_model_matrix`, `stats_lm`, and `stats_glm`.
- Seeded sampling: the 13 `*_rng` functions, with the extension's ChaCha20
  generator and inverse-CDF sampling, so the same seed gives the same draws.
- Posterior-draw kernels: `stats_draw_moments`, `stats_draw_covariance_map`,
  and `stats_draw_affine`, including the fused multiply-add accumulation of
  the extension's matrix product.
- `stats_hello`.
- A randomized differential parity suite that compares the polyfill with the
  loaded extension on seeded random inputs.
- Test suites ported from the extension and run against the R, SciPy, and
  mpmath reference fixtures. A parity suite compares the polyfill with the
  loaded extension.
- `bin/fetch-fixtures`, `bin/generate-bootstrap`, `bin/accuracy-report`, and
  `benchmarks/run.php` with recorded results in `benchmarks/results.json`.
- ACCURACY.md with the largest observed errors per function and reference.

### Extension-only

These functions are left undefined and need the native extension:
`stats_lmm`, `stats_lmm_fit`, `stats_logistic_mixed`, `stats_lmer`,
`stats_glmer`, and `stats_read_csv`.

### Known differences from the extension

- The listwise-deletion warning is `E_USER_WARNING`, because PHP code cannot
  raise `E_WARNING`. The text is identical.
- Scalar arguments of the wrong type follow PHP's rules for userland
  functions.
- The `iterations` of an iterative fit can differ slightly where the
  likelihood is flat to rounding. The estimates still agree.

[0.1.0]: https://github.com/erikwestlund/php-stats-polyfill/releases/tag/v0.1.0
