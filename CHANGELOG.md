# Changelog

All notable changes to this project are documented in this file. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project uses [Semantic Versioning](https://semver.org/).

Each release states the php-stats commit it reproduces. The same commit is
pinned in `tests/fixtures/REFERENCE`, and the reference fixtures are copied
from it.

## [0.1.0] - Unreleased

Reproduces php-stats commit
[`8e6a61e`](https://github.com/erikwestlund/php-stats/commit/8e6a61ee10f0f0986fd51624bfca8912f99fe933).
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
  Results are identical to the extension's to the bit.
- Distributions: 13 families (normal, t, chi-square, F, exponential, logistic,
  log-normal, gamma, beta, binomial, Poisson, negative binomial, and negative
  binomial by mean). Each has the density or mass function, cdf, survival
  function, their logarithms, and the quantile, with `_array` variants (182
  functions).
- Formula interface: `stats_model_matrix`, `stats_lm`, and `stats_glm`, with
  the fitted factor levels and `poly()` coefficients in a `predvars` key.
- `stats_predict` for `stats_lm`, `stats_glm`, `stats_lmer`, and
  `stats_glmer` fits, as R's `predict()`: link or response scale, delta-method
  standard errors, confidence and prediction intervals for linear models, and
  random-effect or population predictions for mixed models.
- Seeded sampling: the 13 `*_rng` functions, with the extension's ChaCha20
  generator and inverse-CDF sampling, so the same seed gives the same draws.
- Posterior-draw kernels: `stats_draw_moments`, `stats_draw_covariance_map`,
  and `stats_draw_affine`, including the fused multiply-add accumulation of
  the extension's matrix product.
- Mixed models: `stats_lmm`, `stats_lmm_fit`, `stats_logistic_mixed`,
  `stats_lmer`, and `stats_glmer`, with exact PHP ports of the extension's
  BOBYQA (bobyqa 0.3.2), Nelder–Mead (argmin 0.11.0), and nalgebra SVD,
  symmetric eigendecomposition, and matrix products.
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

`stats_read_csv` is left undefined and needs the native extension.

### Known differences from the extension

- The listwise-deletion warning is `E_USER_WARNING`, because PHP code cannot
  raise `E_WARNING`. The text is identical.
- Scalar arguments of the wrong type follow PHP's rules for userland
  functions.

[0.1.0]: https://github.com/erikwestlund/php-stats-polyfill/releases/tag/v0.1.0
