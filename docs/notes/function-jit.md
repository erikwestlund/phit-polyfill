# The function JIT miscompiles code entered partway through a function

Recorded 2026-10-16 on PHP 8.5.10 (macOS arm64) and 8.5.11 (GitHub's Ubuntu
x86-64 and macOS runners, and the official `php:8.5-cli` image under x86-64
emulation). The upstream issue is
[php-src #24214](https://github.com/php/php-src/issues/24214), opened
2026-10-09 and still open. Nothing from this project has been added to it.

PHP 8.4.25 (Herd, macOS arm64) has the same bug. With the code from before the
fixes below, `opcache.jit=1235` sent the NelderMead test closure into the same
endless loop in three of four hot-counter settings and printed a `live var`
line each time. The fixed code passes on PHP 8.4 under every JIT setting.
PHP 8.3 has a different JIT with different problems; see
[php83-jit.md](php83-jit.md).

## Recommendation

On PHP 8.4 and 8.5, run the polyfill without the JIT, with the tracing JIT
(`opcache.jit=tracing`), or with PHP's `function` mode (`opcache.jit=function`,
which is `1205`). Avoid a function JIT with hot-counter triggers, such as
`opcache.jit=1235`, the value that setup-php's php.ini sets. Under that setting
the polyfill has produced wrong singular values and an endless loop in test
code, as described below.

## The bug

The function JIT compiles a whole function at once. The compiled code can be
entered from the interpreter partway through the function:

- at a loop header, when a loop becomes hot before its function does (an
  on-stack-replacement, or OSR, entry). Only the hot-counter trigger
  (`opcache.jit` digit T = 3, as in `1235`) creates these.
- after a call that may have run in the interpreter, such as a call resolved
  by name at run time (a continuation entry).

At such an entry the JIT can restore only values that are bound to a slot of
the PHP call frame. Some values live only in a CPU register. In some of these
cases PHP prints

```
ENTRY n (block b start s) - live var v
```

to stderr while compiling, where `n` is the opcode index of the entry point.
In others it prints nothing. In both cases the compiled code reads whatever
the register holds. The maintainer's analysis in #24214 traces one cause to an
assignment such as `$n = count($a)` and the spill-conflict check added for
GH-16821.

## What happened here

CI jobs that enabled OPcache ran `opcache.jit=1235`, because setup-php sets it.

1. **An endless loop in a test closure.** The cost function in
   `tests/Unit/NelderMeadTest.php` computed `$n = count($y)` and looped
   `for ($i = 0; $i < $n; $i++)`. Entered through an OSR entry, `$n` was
   garbage, the loop ran past the array, and the test emitted warnings until
   memory ran out. This was the "hang" after `NanSafetyTest`. It depended on
   hot-counter timing and address-space layout, so it was intermittent.

2. **Wrong singular values, silently.** `Svd::reflectRows` applies a
   Householder reflection in nested loops. When an inner loop over rows
   became hot first, a value of the outer loop was not restored and one entry
   of Vᵀ came out as zero. No message was printed. `DenseTest` failed in
   about one run in ten with the default hot counters, and deterministically
   with `opcache.jit_hot_loop=10`.

3. **Harmless entries in `Cdflib::betaFrac`.** The JIT reported five
   continuation entries after the by-name calls `fdiv()`, `abs()`, and
   `is_nan()`. Those calls resolve to internal functions, which never return
   through the interpreter, so the entries were not taken and results were
   correct.

Fully qualifying every internal-function call in `src/` (`\abs()`) was tried as
a general fix and reverted (commit e902431, reverted in 1b7b95f). It removed
the by-name entries, but it changed the generated code everywhere. On arm64
it exposed three more miscompiles (in `betaFrac`, `Dense::multiply`, and
`LeastSquares::singularValues`), and on x86-64 the tracing JIT and
`opcache.jit=1205` then computed wrong incomplete-beta values in
`Cdflib::betaGrat`, which the unqualified code does not.

## What the polyfill does

- The NelderMead test closure computes its loop bound without assigning a
  temporary and iterates the quadratic case with `foreach`.
- `Svd::reflectRows` moves its loops over rows into separate functions, so an
  OSR entry there has only parameters and local variables live. The arithmetic
  is unchanged.
- `betaFrac` calls `\fdiv()`, `\abs()`, and `\is_nan()` fully qualified and
  stores the new ratio in a fresh variable before copying the old one
  (`$next = \fdiv(...); $r0 = $r; $r = $next;`). Without the reordering the
  qualified version computed wrong values on arm64, because the JIT let the
  copy `$r0` see the call result written into `$r`.
- CI runs the suite without OPcache, with the optimizer alone, with the
  tracing JIT, and with `opcache.jit=1205` on Linux and macOS. These jobs fail
  on a test failure and on any `- live var` line. A job with
  `opcache.jit=1235` runs as well, but it does not block, because further OSR
  miscompiles cannot be ruled out.

## Finding affected code

- `opcache.jit_debug=0x1000000` prints each function's IR before compiling it,
  headed `JIT$<function>: ; after folding`. The most recent header before an
  `ENTRY ... live var` line names the function. The printing changes timing,
  so it does not always reproduce what a plain run does.
- `opcache.jit_bisect_limit=N` compiles only the first `N` functions and names
  the first one it skips. A binary search on `N` finds a function that
  computes a wrong result without any message. The compile order varies
  between runs of the full suite, so the search needs a standalone script
  whose order is fixed.
- `opcache.jit_hot_loop` set to small values (1, 3, or 10) changes which loops
  become hot first and turns some intermittent OSR failures into
  deterministic ones.

## Verification

The suite (552 tests, and 656 with the extension's parity tests) was run on
arm64 and on x86-64 under emulation with OPcache off, the optimizer alone, the
tracing JIT, `opcache.jit=1205`, `1215`, and `1235` with default, low, and
OSR-first hot counters. Every run passed without a `live var` line, including
the `1235` runs. A sweep of 80 large inputs (regressions up to 120 columns,
formula and mixed models, prediction, distribution grids, samplers, and draw
kernels) gave results identical to the run without OPcache under all of those
settings on arm64, and under the tracing JIT and `1205` on x86-64. Passing
`1235` runs do not rule out other OSR miscompiles, which depend on which loops
become hot first.
