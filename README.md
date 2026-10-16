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
Polyfill::supports('stats_lmm');       // false: needs the extension
Polyfill::nativeLoaded();              // whether php_stats is loaded
```

## Coverage

| Area | Functions | Polyfill |
|---|---|---|
| Descriptive statistics | `stats_mean`, `stats_var`, `stats_sd`, `stats_median`, `stats_quantile`, `stats_quantile_array`, `stats_cov`, `stats_cor`, `stats_weighted_mean`, `stats_weighted_var` | Provided |
| Regression | `stats_ols`, `stats_logistic`, `stats_probit`, `stats_poisson`, with weights, offsets, the automatic intercept, and both missing-value policies | Provided |
| Distributions | See [Distributions](#distributions) | Provided where accurate |
| Mixed models | `stats_lmm`, `stats_lmm_fit`, `stats_logistic_mixed` | Extension only |
| Random numbers and draws | every `*_rng` function, `stats_draw_affine`, `stats_draw_covariance_map`, `stats_draw_moments` | Extension only |
| Other | `stats_hello` | Extension only |

Functions marked "Extension only" stay undefined without the extension, so
`function_exists()` reports them correctly.

### Distributions

Distribution functions are in progress and will be listed here.

## How it works

The implementations are ports of the extension's Rust code, not independent
reimplementations, so results agree to rounding in most cases.

- **Least squares.** `stats_ols` scales the design columns to unit norm,
  factors the design by Householder QR, and takes a one-sided Jacobi SVD of the
  small R factor. Those singular values equal the design's, so the extension's
  rank rule (s ≤ ε·max(n, p)·s_max) makes the same rank decisions.
- **Generalized linear models.** `stats_logistic` and `stats_poisson` use
  damped Newton IRLS and `stats_probit` uses Fisher scoring. Each iteration
  solves a weighted least-squares step with the same QR and SVD, halves the step
  until the likelihood does not decrease, and stops only after a small full
  step. The separation certificate, boundary handling, starting values for
  extreme offsets, and exception codes match the extension.
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
- OLS agrees with the extension to about 10⁻¹⁰ relative on the reference
  fixtures. On designs with condition numbers up to 10¹⁰ it agrees to about
  cond(X)·ε, which is what two backward-stable solvers can promise.
- Logistic, probit, and Poisson regression agree with the extension to about
  10⁻⁸, the level of their shared convergence rule, and usually to 10⁻¹².

### Known differences

- The listwise-deletion warning uses `E_USER_WARNING` because PHP code cannot
  raise `E_WARNING`. The text is identical.
- Scalar arguments of the wrong type follow PHP's userland rules (coercion in
  non-strict files, otherwise `TypeError`). The extension raises an
  `\Exception` instead.
- The `iterations` reported by an iterative fit can differ by one, and by a
  little more when the likelihood is flat to rounding near the optimum. The
  estimates still agree.

## Performance

Benchmarks compare the polyfill with the extension in the same process on an
Apple M5 Ultra with PHP 8.5.10. Each design has 10 standard normal predictors
plus an intercept. Times are the median of five runs with the tracing JIT on
(`-d opcache.enable_cli=1 -d opcache.jit=tracing -d opcache.jit_buffer_size=64M`).

| Model | n | Polyfill | Extension | Ratio |
|---|---:|---:|---:|---:|
| `stats_ols` | 10,000 | 6.7 ms | 2.8 ms | 2.4× |
| `stats_ols` | 100,000 | 72 ms | 28 ms | 2.6× |
| `stats_logistic` | 10,000 | 38 ms | 6.9 ms | 5.5× |
| `stats_logistic` | 100,000 | 392 ms | 70 ms | 5.6× |

The logistic fit includes the null model that every intercept model fits for
`null_deviance`. Without the JIT the polyfill is about three times slower
(226 ms for OLS and 1.25 s for logistic regression at n = 100,000).
For comparison, solving the normal equations (X'X)b = X'y in plain PHP took
20 ms at n = 100,000. That baseline returns coefficients only, without
standard errors, rank detection, or the other results. Forming X'X also
squares the condition number, so it loses accuracy on ill-conditioned designs.
The polyfill's coefficients agreed with the extension's to within 2·10⁻¹² standard
errors in every case. The recorded run, with its full context (CPU, platform,
PHP and JIT settings, median and minimum times, a plain-PHP normal-equations
baseline for OLS, and the agreement check), is in
[benchmarks/results.json](benchmarks/results.json). Reproduce it with:

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
