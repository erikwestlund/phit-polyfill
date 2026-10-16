# php-stats-polyfill

A pure-PHP polyfill for the [php_stats](https://github.com/erikwestlund/php-stats)
native extension. It defines the extension's `stats_*` functions only when the
extension is not loaded, with the same names, parameter names, order and
defaults, result keys, warnings, and exception codes. Code written against
php_stats therefore runs unchanged on hosts where the extension cannot be
installed. When the extension is loaded, it is used and the polyfill does
nothing.

The package is part of the [phpstats](https://phpstats.org) family, alongside
phibble, phplyr, phpolars, and cmdstanphp. It follows the pattern of
`symfony/polyfill-mbstring`.

## Versions

Each polyfill release reproduces one php-stats commit, which is pinned in
`tests/fixtures/REFERENCE`. The reference fixtures are copied from it.

| Polyfill | php-stats | Notes |
|---|---|---|
| 0.1.0 | [`91ee538`](https://github.com/erikwestlund/php-stats/commit/91ee538253e3ecf7b3f086441305c3f7237200fc) | php-stats has no tagged release yet |

See [CHANGELOG.md](CHANGELOG.md) for what each release contains.

## Installation

```bash
composer require erikwestlund/php-stats-polyfill
```

PHP 8.5 or later is required. Composer loads `bootstrap.php` automatically. For
each function the polyfill implements, it declares the global function only if
`function_exists()` reports that it is missing.

```php
$fit = stats_ols(
    x: [[1.2, 0.4], [2.3, 1.1], [3.1, 0.9], [4.8, 2.2], [5.0, 1.7]],
    y: [2.1, 4.2, 5.8, 8.1, 10.2],
);

$fit['coefficients'];  // intercept first, as with the extension
$fit['p_values'];
```

To check at runtime whether a function is available without the extension,
use `PhpStats\Polyfill\Polyfill`:

```php
use PhpStats\Polyfill\Polyfill;

Polyfill::supports('stats_logistic');  // true
Polyfill::supports('stats_read_csv');  // false: needs the extension
Polyfill::nativeLoaded();              // whether php_stats is loaded
```

## Coverage

| Area | Functions | Polyfill |
|---|---|---|
| Descriptive statistics | `stats_mean`, `stats_var`, `stats_sd`, `stats_median`, `stats_quantile`, `stats_quantile_array`, `stats_cov`, `stats_cor`, `stats_weighted_mean`, `stats_weighted_var` | Provided |
| Regression | `stats_ols`, `stats_logistic`, `stats_probit`, `stats_poisson`, with weights, offsets, the automatic intercept, and both missing-value policies | Provided |
| Formula interface | `stats_model_matrix`, `stats_lm`, `stats_glm` | Provided |
| Distributions | 13 families × 14 functions; see [Distributions](#distributions) | Provided |
| Mixed models | `stats_lmm`, `stats_lmm_fit`, `stats_logistic_mixed`, `stats_lmer`, `stats_glmer` | Provided (slower than the extension; see [ACCURACY.md](ACCURACY.md#speed)) |
| Random numbers | 13 seeded `*_rng` functions, one per distribution family | Provided |
| Posterior draws | `stats_draw_moments`, `stats_draw_covariance_map`, `stats_draw_affine` | Provided |
| CSV reading | `stats_read_csv` (phibble's PHP reader is the fallback) | Extension only |
| Other | `stats_hello` | Provided |

Functions marked "Extension only" stay undefined without the extension, so
`function_exists()` reports them correctly.

### Distributions

All 13 distribution families are provided, each with the full set of
operations and their `_array` variants:

| Family | Prefix | Parameters | Operations |
|---|---|---|---|
| Normal | `stats_normal_` | mean, sd | pdf, logpdf, cdf, logcdf, sf, logsf, quantile |
| Student's t | `stats_t_` | df | same |
| Chi-square | `stats_chisquare_` | df | same |
| F | `stats_f_` | df1, df2 | same |
| Exponential | `stats_exponential_` | rate | same |
| Logistic | `stats_logistic_` | location, scale | same |
| Log-normal | `stats_lognormal_` | meanlog, sdlog | same |
| Gamma | `stats_gamma_` | shape, scale | same |
| Beta | `stats_beta_` | alpha, beta | same |
| Binomial | `stats_binomial_` | size, prob | pmf, logpmf, cdf, logcdf, sf, logsf, quantile |
| Poisson | `stats_poisson_` | lambda | same |
| Negative binomial | `stats_nbinom_` | size, prob | same |
| Negative binomial by mean | `stats_nbinom_mu_` | size, mu | same |

For example, `stats_t_sf(2.5, df: 12)` and
`stats_normal_quantile_array([0.025, 0.975])`.

Each family also has a seeded sampler, such as `stats_normal_rng(1000, seed: 42)`.
The polyfill uses the extension's generator (ChaCha20 seeded as in rand_core,
with inverse-CDF sampling), so a seed gives the same draws, bit for bit, with
or without the extension.

### Formula interface

`stats_model_matrix`, `stats_lm`, and `stats_glm` accept an R-style formula
over named columns. They port the extension's parser, term expansion, and
coding, so column names, term order, contrasts, and missing-value handling
match the extension. Where the extension differs from R, the polyfill follows
the extension (see the php-stats formula documentation).

```php
$data = [
    'rent' => [1879.25, 1805.03, 2022.29, 1267.85, 1509.24, 2240.10],
    'area' => [74.4, 57.4, 66.4, 49.1, 45.7, 92.4],
    'district' => ['north', 'south', 'north', 'east', 'east', 'south'],
];

$fit = stats_lm('rent ~ log(area) + district', $data);
$fit['named'];    // coefficients keyed by column name

$design = stats_model_matrix('~ area * district', $data);
$design['columns'];

$logit = stats_glm('y ~ x + offset(log(t))', $other, family: 'binomial');
```

Supported syntax includes `+ - : * / %in% ^`, `0`/`1`/`-1`, `.`, `offset()`,
`I()`, the common transforms (`log`, `log2`, `log10`, `log1p`, `exp`, `sqrt`,
`abs`, and others), `poly()`, `factor()`, `ordered()`, and `cbind()` responses
for grouped binomial data. Families are gaussian, binomial (logit or probit),
and poisson. Formulas with random-effect bars such as `(1 | g)` are accepted
by `stats_model_matrix`, which returns the random-effect design. `stats_lm`
and `stats_glm` reject them with code 1007, because `stats_lmer` and
`stats_glmer` need the extension.

## How it works

The implementations are ports of the extension's Rust code, not independent
reimplementations, so results agree to rounding in most cases.

- **Least squares.** `stats_ols` scales the design columns to unit norm and
  factors the design by Householder QR, using the same operation order as the
  extension's linear-algebra library. The singular values of the small R
  factor (computed by a Jacobi sweep) equal the design's, so the extension's
  rank rule (s ≤ ε·max(n, p)·s_max) makes the same rank decisions.
  Coefficients and covariances come from R by back-substitution. Weighted fits
  factor rows in decreasing-weight order, as the extension does.
- **Generalized linear models.** `stats_logistic` and `stats_poisson` use
  damped Newton IRLS and `stats_probit` uses Fisher scoring. Each iteration
  solves a weighted least-squares step with the same QR, halves the step
  until the likelihood does not decrease, and stops only after a small full
  step. The separation certificate, boundary handling, starting values for
  extreme offsets, and exception codes match the extension.
- **Distributions.** The extension's distribution layer is ported together
  with the parts of cdflib it calls: TOMS 708 (incomplete beta), TOMS 654
  (incomplete gamma), Cody's normal CDF, and the `dinvnr`/`dt1` quantile
  starts. Log-scale tails, saddle-point densities, and the safeguarded quantile
  search are ported as well.
- **Special functions.** The statrs routines that the extension's estimators
  call (erfc, ln Γ, the regularized incomplete beta, and the inverse error
  function) are ported from the same crate version.
- **Validation.** Inputs are validated in the extension's order, so an input
  that is invalid in several ways raises the same exception code.

Exceptions are plain `\Exception` instances whose code identifies the category,
as in the extension:

| Code | Meaning |
|---:|---|
| 1001 | Invalid input |
| 1002 | Dimension mismatch |
| 1003 | Non-finite value or overflow |
| 1004 | Rank deficiency |
| 1005 | Factorization failure |
| 1006 | Non-convergence (including separation) |
| 1007 | Unsupported input size |
| 1008 | Internal error |

## Accuracy

Every provided function passes the extension's R, SciPy, and mpmath fixtures
at the extension's own tolerances, and a parity suite that compares the
polyfill with the loaded extension. [ACCURACY.md](ACCURACY.md) lists, per
function, the largest observed error against each reference and the tested
domain.

In summary:

- Descriptive statistics agree with the extension exactly on every input
  tested, including values offset by 10¹².
- OLS, logistic, probit, and Poisson regression return the same bits as the
  extension on every tested input, including designs with condition numbers
  up to 10¹⁰, extreme prior weights, separation, and extreme offsets.
- Distribution functions are bit-identical to the extension on every input
  tested, and within 10⁻¹² of mpmath down to probabilities near 10⁻³⁰⁰.
- The formula functions produce identical design matrices, column names, and
  terms, and identical fits.
- Seeded samplers, posterior-draw kernels, and the five mixed-model functions
  are bit-identical to the extension in the randomized parity suite.

### Known differences

- The listwise-deletion warning uses `E_USER_WARNING` because PHP code cannot
  raise `E_WARNING`. The text is identical.
- Scalar arguments of the wrong type follow PHP's userland rules (coercion in
  non-strict files, otherwise `TypeError`). The extension raises an
  `\Exception` instead.
- Matrix products above 5×5 follow the extension's aarch64 build, which uses
  fused multiply-adds. On an x86-64 CPU without FMA the extension's products
  can differ in the last bit.

## Performance

Benchmarks compare the polyfill with the extension in the same process on an
Apple M5 Ultra with PHP 8.5.10. Each design has 10 standard normal predictors
plus an intercept. Times are the median of five runs with the tracing JIT on
(`-d opcache.enable_cli=1 -d opcache.jit=tracing -d opcache.jit_buffer_size=64M`).

| Model | n | Polyfill | Extension | Ratio |
|---|---:|---:|---:|---:|
| `stats_ols` | 10,000 | 7.2 ms | 2.6 ms | 2.8× |
| `stats_ols` | 100,000 | 81 ms | 28 ms | 2.9× |
| `stats_logistic` | 10,000 | 43 ms | 6.7 ms | 6.5× |
| `stats_logistic` | 100,000 | 446 ms | 69 ms | 6.4× |

The logistic fit includes the null model that every intercept model fits for
`null_deviance`. Without the JIT the polyfill is about three times slower:
247 ms for OLS and 1.46 s for logistic regression at n = 100,000.

For comparison, solving the normal equations (X'X)b = X'y in plain PHP took
20 ms at n = 100,000. That baseline returns coefficients only, without
standard errors, rank detection, or the other results. Forming X'X also
squares the condition number, so it loses accuracy on ill-conditioned designs.

The polyfill's coefficients were identical to the extension's in every case.
The recorded run is in [benchmarks/results.json](benchmarks/results.json),
with its full context: CPU, platform, PHP and JIT settings, median and minimum
times, the plain-PHP baseline for OLS, and the agreement check. Reproduce it
with:

```bash
php -d extension=/path/to/libphp_stats.so -d opcache.enable_cli=1 \
    -d opcache.jit=tracing -d opcache.jit_buffer_size=64M \
    benchmarks/run.php --sizes=10000,100000 --repeat=5 --seed=5 \
    --json=benchmarks/results.json
```

## Development

```bash
composer install
composer test        # polyfill only
composer test:ext    # with the extension loaded from .ext/libphp_stats.dylib
```

The composer scripts use Herd's `php85`. Run `vendor/bin/pest` with any PHP 8.5
binary otherwise.

- `bin/fetch-fixtures [checkout]` copies the reference fixtures from the
  php-stats commit pinned in `tests/fixtures/REFERENCE`. The copies are
  committed and never edited by hand.
- `bin/generate-bootstrap` regenerates `bootstrap.php` and
  `src/FunctionMap.php` from the extension's API contract
  (`tests/fixtures/php-api.json`). A function is declared only when its
  implementation exists.
- `bin/accuracy-report` prints the tables in ACCURACY.md.
- `tests/Unit` holds the ports of the extension's suites. They call the
  polyfill directly, so they test it whether or not the extension is loaded.
  `tests/Parity` runs only with the extension and compares the two.

## License

MIT. See [LICENSE](LICENSE).
