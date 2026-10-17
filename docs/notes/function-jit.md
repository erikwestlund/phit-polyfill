# The function JIT loses values at entry points

Recorded 2026-10-16 on PHP 8.5.10 (macOS arm64) and 8.5.11 (GitHub's Ubuntu
x86-64 and macOS runners). The upstream issue is
[php-src #24214](https://github.com/php/php-src/issues/24214), opened
2026-10-09 and still open. Nothing from this project has been added to it.

## Summary

PHP's function JIT (`opcache.jit=1235`, and other `opcache.jit` values with a
trigger digit from 0 to 3) compiles a whole function at once. The compiled code
can be entered from the interpreter partway through the function:

- at a loop header, when a loop becomes hot before its function does (an
  on-stack-replacement, or OSR, entry; trigger 3 only);
- after a call that may have run in the interpreter, such as a call to a user
  function or a call resolved by name at run time (a continuation entry).

At such an entry the JIT can restore only values that are bound to a slot of
the PHP call frame. If its register allocator keeps a live value only in a CPU
register, PHP prints

```
ENTRY n (block b start s) - live var v
```

to stderr and emits code that reads whatever the register holds. `n` is the
opcode index of the entry point. The maintainer's analysis in #24214 traces
one cause to `$n = count($a)` and the spill-conflict check added for GH-16821.

The tracing JIT (`opcache.jit=tracing`) and OPcache without the JIT were not
affected in any test below.

## What happened here

setup-php's php.ini sets `opcache.jit=1235`, so CI jobs that enabled OPcache
ran the function JIT. On commit cb0395c the only wrong behavior was problem 1,
in test code. The JIT also reported five entries in `Cdflib::betaFrac`, but
those followed by-name calls to internal functions, which never return
through the interpreter, and a large-input sweep matched the run without
OPcache. Problems 2 to 4 appeared after the internal-function calls in `src/`
were fully qualified, which changes the code the JIT generates, and were fixed
before that change was committed.

1. **A hang in a test closure.** The cost function in
   `tests/Unit/NelderMeadTest.php` computed `$n = count($y)` and looped
   `for ($i = 0; $i < $n; $i++)`. Entered through an OSR entry at the loop
   header, `$n` was garbage, the loop ran past the array, and the test emitted
   warnings until memory ran out. This was the "hang" after `NanSafetyTest`.
   Whether it happened depended on hot-counter timing and address-space
   layout, so it was intermittent.

2. **A copy that sees a later call result.** In `Cdflib::betaFrac` the
   continued-fraction loop does `$r0 = $r; $r = \fdiv($anp1, $bnp1);` and then
   tests `abs($r - $r0)`. When the call is compiled as a direct internal call,
   the JIT lets `$r0` share `$r`'s value and then writes the call result into
   `$r`'s frame slot, so `$r0` sees the new value, the difference is zero, and
   the series stops after one step. The beta, binomial, negative binomial, and
   F distributions returned values off by up to 10⁻⁴ relative, with no message
   on stderr. The same happens with a user-function call. A standalone
   reproducer gives the wrong result on the first call with
   `opcache.jit=1205` (compile everything at load).

3. **`\count()` at an OSR entry.** A fully qualified `\count($a)` compiles to
   the COUNT opcode. The JIT inlines it as a 32-bit load followed by a zero
   extension and binds only the extended value to the variable's slot. At an
   OSR entry the 32-bit value is reported as a live var and is not restored.
   In `Dense::multiply`, which computes `$k = \count($a)` and
   `$n = \count($b)` before its loops, this made matrix products loop
   forever. An unqualified `count()` in a namespace is an ordinary call and is
   not affected.

4. **A value of an outer loop at an inner loop's OSR entry.** The one-sided
   Jacobi sweep in `LeastSquares::singularValues` has three nested loops.
   When an inner loop over `$k` became hot first, two values of the outer
   loops were reported as live vars and were not restored. With 70 or more
   columns the singular values differed from the run without the JIT. The
   unit tests use smaller matrices and did not show it.

## Finding affected code

- `opcache.jit_debug=0x1000000` prints each function's IR before compiling it,
  headed `JIT$<function>: ; after folding`. The most recent header before an
  `ENTRY ... live var` line names the affected function. This printing changes
  timing, so it does not always reproduce what a plain run does.
- `opcache.jit_bisect_limit=N` compiles only the first `N` functions and names
  the first one it skips. A binary search on `N` finds a function that
  computes a wrong result without any message.
- `opcache.jit_hot_loop=1` with `opcache.jit_hot_func=255` compiles nearly
  every function through an OSR entry at its first loop, which turns the
  intermittent OSR failures into deterministic ones.
- `opcache.opt_debug_level=0x20000` dumps the optimized opcodes. Problem 2
  shows up there as `CVx = QM_ASSIGN CVy` followed by `CVy = DO_ICALL` or
  `DO_UCALL` while `CVx` is still read. A scan of `src/` for that shape found
  `betaFrac` and two sites in `TrustRegion::trsbox`.

## What the polyfill does

- Calls to internal functions in `src/` are fully qualified (`\abs($x)`), which
  removes the by-name continuation entries. `count()` is the exception and
  stays unqualified because of problem 3. `tests/Unit/JitSafetyTest.php`
  enforces both rules.
- `betaFrac` and `trsbox` store a call result in a fresh variable before
  copying the old value, for example
  `$next = \fdiv(...); $r0 = $r; $r = $next;`.
- The inner loops of `singularValues` are separate functions
  (`jacobiDots` and `jacobiRotate`), so their OSR entries have only
  parameters and local accumulators live.
- The NelderMead test closure computes its loop bound without assigning a
  temporary and iterates the quadratic case with `foreach`.
- CI runs the suite under `opcache.jit=1235` on Linux and macOS and under
  `opcache.jit=1205`. These jobs fail on a test failure and also on any
  `- live var` line in the output.

Besides the test and parity suites, a sweep of 80 large inputs (regressions
up to 120 columns, formula models, mixed models, prediction, dense
distribution grids, samplers, and draw kernels) gave results identical to the
run without OPcache under `opcache.jit=1235` with default, low, and OSR-first
hot counters, under `opcache.jit=1205`, and under the tracing JIT, with no
`live var` lines.

These changes cover the code paths that the suites and the sweep exercise. Other code shapes can trigger the
same bug, so a user who sees a `live var` line under the function JIT, or a
result that differs from the tracing JIT, should report it here.
