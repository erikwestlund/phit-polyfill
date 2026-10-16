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

## CI note: setup-php enables the function JIT

setup-php's php.ini sets `opcache.jit=1235`, so a CI job that only passes
`-d opcache.enable_cli=1` runs the function JIT, not the optimizer alone.
With that setting PHP 8.5.11 prints `ENTRY n (block b start s) - live var v`
while compiling large functions, and on 2026-10-16 two such jobs on commit
04110a9 hung after `NanSafetyTest`. The hang did not recur on the next commit,
and it could not be reproduced locally with the official Docker image, the
ondrej/php build, or Homebrew's 8.5.10, with or without `opcache.jit=1235`.
The workflow now disables the JIT in the `opcache` job, tests `opcache.jit=1235`
in a separate `function-jit` job, and times every job out after 20 minutes.

The `function-jit` job on commit 298cdb0 hung at the same point, after printing
`ENTRY 102 (block 580 start 3111) - live var 2930`. The message appears to come
from the JIT's register allocator, which counts IR instructions, so the
function being compiled is probably one of the four large BOBYQA routines
(`Bobyqb::solve`, `Rescue::rescue`, `TrustRegion::trsbox`, and
`Geometry::geostep`, each with 1,000 to 1,400 opcodes). With
`opcache.jit=1235` the JIT compiles functions once they become hot, which
would explain why the hang is intermittent. The optimizer alone
(`opcache.jit=disable`) and the tracing JIT pass. The job is non-blocking until
the cause is confirmed. Nothing has been reported upstream.
