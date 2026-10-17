# OPcache folds negated float comparisons incorrectly for NaN

Recorded 2026-10-16 while porting the php_stats distribution quantiles. Nothing
has been reported upstream.

## Summary

With OPcache enabled, PHP's optimizer rewrites a negated comparison such as
`!($x > 0.0)` into the inverted comparison `$x <= 0.0`. The two are equivalent
for ordinary numbers but not for NaN. Every comparison with NaN is false, so
`!(NAN > 0.0)` is `true` while `NAN <= 0.0` is `false`. A script therefore
gives different results depending on whether OPcache optimized the file.

The rewrite happens in optimizer pass 5 (the CFG-based block pass, bit `0x10`
of `opcache.optimization_level`). The JIT is not involved. The result is the
same with the JIT disabled, in function mode, and in tracing mode.

## Environment

- PHP 8.5.10 (cli, NTS, clang 15.0.0, built 2026-09-24), from Laravel Herd.
- PHP 8.4.25 from the same distribution behaves the same way.
- macOS 27.0 on arm64 (Apple M5 Ultra).
- No `php.ini` is loaded. Herd's binary defaults are `opcache.enable_cli=On`,
  `opcache.optimization_level=0x7FFEBFFF`, `opcache.jit=disable`,
  `opcache.file_update_protection=2`, and no `opcache.file_cache`.

## Reproduction

```php
<?php
function not_greater(float $x): bool { return ! ($x > 0.0); }
function not_less(float $x): bool { return ! ($x < 0.0); }
function not_ge(float $x): bool { return ! ($x >= 0.0); }
function not_le(float $x): bool { return ! ($x <= 0.0); }
function untyped($x) { return ! ($x > 0.0); }
function branch(float $x): bool { if (! ($x > 0.0)) { return true; } return false; }
function local(): bool { $x = NAN; return ! ($x > 0.0); }
function eq(float $x): bool { return ! ($x == 0.0); }

echo json_encode([not_greater(NAN), not_less(NAN), not_ge(NAN), not_le(NAN),
    untyped(NAN), branch(NAN), local(), eq(NAN)]), "\n";
```

Every value should be `true`. Save the script and wait at least two seconds
before running it. With `opcache.file_update_protection=2`, OPcache does not
cache or optimize a file modified less than two seconds earlier, so a script
run immediately after it is written does not show the problem.

| Settings | Output |
|---|---|
| `-d opcache.enable_cli=0` | `[true,true,true,true,true,true,true,true]` |
| defaults (`opcache.enable_cli=On`) | `[false,false,false,false,false,false,false,true]` |
| `-d opcache.optimization_level=0x10` (pass 5 only) | same as defaults |
| `-d opcache.optimization_level=0x7FFEBFEF` (all but pass 5) | all `true` |
| `-d opcache.jit=function` or `tracing`, `-d opcache.jit_buffer_size=64M` | same as defaults |
| tracing JIT with `opcache.optimization_level=0x7FFEBFEF` | all `true` |
| PHP 8.4.25, `-d opcache.enable_cli=1` | same as PHP 8.5 defaults |

Equality is not affected: `!($x == 0.0)` stays `true`, presumably because the
rewrite to `!=` is also correct for NaN.

Each optimization bit was tested on its own (`opcache.optimization_level=1 << n`).
Only bit 4 (`0x10`) reproduces the problem.

## Opcodes

`-d opcache.opt_debug_level=0x10000` (before the optimizer) and `0x20000`
(after) show the rewrite for `not_greater()`:

```
before                                   after
0001 T1 = IS_SMALLER float(0) CV0($x)    0001 T1 = IS_SMALLER_OR_EQUAL CV0($x) float(0)
0002 T2 = BOOL_NOT T1                    0002 RETURN T1
0003 VERIFY_RETURN_TYPE T2
0004 RETURN T2
```

`BOOL_NOT(0 < x)` becomes `x <= 0`. The transformation assumes that `<` is a
total order, which does not hold for floats when NaN is involved.

## How it showed up here

The polyfill's continuous quantile search ports the extension's bracket
expansion. In php-stats bd98bf7 the extension stops the search when a tail
evaluation returns NaN. The port first expressed the stop condition as
`!($g($v) > 0.0)`. Under OPcache that became `$g($v) <= 0.0`, which is false
for NaN. As a result the search kept walking to −∞ and returned `-INF`, while
the extension raised exception 1006. The case was
`stats_t_quantile(1 - 0.999999999999, 1e299)`. With OPcache disabled the
polyfill matched the extension.

The problem first looked intermittent, because tests that ran within two
seconds of editing a file used the unoptimized code.

## What the polyfill does

- Code under `src/` does not apply `!` to a comparison. NaN-sensitive
  conditions test `is_nan()` explicitly, and range checks are written as
  `is_nan($x) || $x < $lo || $x > $hi`. Conditions whose negation matters for
  NaN use an `if`/`else` instead of `!`.
- `tests/Unit/NanSafetyTest.php` fails if a `!` applied to a parenthesised
  comparison appears in `src/`.
- The full test suite is run with OPcache disabled, with OPcache enabled, and
  with the tracing JIT. Source files must be at least two seconds old for the
  OPcache runs to exercise the optimizer.

A separate function-JIT problem found while setting up CI is described in
[function-jit.md](function-jit.md).
