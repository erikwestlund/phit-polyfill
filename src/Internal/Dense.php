<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

/**
 * Dense matrix kernels of nalgebra 0.35.0 for dynamically sized f64 matrices,
 * reproduced operation for operation so results match the extension's bits.
 *
 * Matrices are column-major lists of columns: $a[$j][$i] is row i of column j.
 * Vectors are plain lists.
 *
 * @internal
 */
final class Dense
{
    /** nalgebra's gemm uses matrixmultiply only when every dimension exceeds this. */
    private const int SMALL_DIM = 5;

    /** matrixmultiply's depth block for f64 on aarch64 (the NEON kernel uses S_KC = 256). */
    private const int KC = 256;

    /** Below about 2^-850, Fma::fma's error term may lose bits to subnormals. */
    private const float TINY_PRODUCT = 7.0e-257;

    /** Fma::fma's limit for exact splitting, also kept clear of overflow in its sums. */
    private const float HUGE = 1e300;

    /**
     * nalgebra's dotx for dynamic vectors over indices [$start, $end): eight
     * interleaved partial sums over blocks of eight, combined pairwise, then
     * the remainder in order.
     *
     * @param  list<float>  $x
     * @param  list<float>  $y
     */
    public static function dot(array $x, array $y, int $start = 0, ?int $end = null): float
    {
        $end ??= count($x);
        $acc0 = $acc1 = $acc2 = $acc3 = $acc4 = $acc5 = $acc6 = $acc7 = 0.0;
        $i = $start;
        while ($end - $i >= 8) {
            $acc0 += $x[$i] * $y[$i];
            $acc1 += $x[$i + 1] * $y[$i + 1];
            $acc2 += $x[$i + 2] * $y[$i + 2];
            $acc3 += $x[$i + 3] * $y[$i + 3];
            $acc4 += $x[$i + 4] * $y[$i + 4];
            $acc5 += $x[$i + 5] * $y[$i + 5];
            $acc6 += $x[$i + 6] * $y[$i + 6];
            $acc7 += $x[$i + 7] * $y[$i + 7];
            $i += 8;
        }
        $result = 0.0;
        $result += $acc0 + $acc4;
        $result += $acc1 + $acc5;
        $result += $acc2 + $acc6;
        $result += $acc3 + $acc7;
        for (; $i < $end; $i++) {
            $result += $x[$i] * $y[$i];
        }

        return $result;
    }

    /**
     * nalgebra's `a * b` (gemm with alpha 1, beta 0). When any of the result's
     * rows, its columns, or the inner dimension is at most 5 it accumulates
     * columns gemv-style (y = a0 * b0, then y = ak * bk + y, unfused);
     * otherwise each entry is matrixmultiply's NEON kernel: a fused
     * multiply-add chain over k from +0.0 in blocks of 256, with later blocks
     * added to the stored partial result.
     *
     * @param  list<list<float>>  $a  m x k, as columns
     * @param  list<list<float>>  $b  k x n, as columns
     * @param  int|null  $rows  m, needed only when $a has no columns
     * @return list<list<float>> m x n, as columns
     */
    public static function multiply(array $a, array $b, ?int $rows = null): array
    {
        $k = count($a);
        $n = count($b);
        $m = $rows ?? ($k > 0 ? count($a[0]) : 0);
        if ($m > self::SMALL_DIM && $n > self::SMALL_DIM && $k > self::SMALL_DIM) {
            return self::multiplyFused($a, $b, $m, $k, $n);
        }
        $c = [];
        for ($j = 0; $j < $n; $j++) {
            $bj = $b[$j];
            if ($k === 0) {
                $c[$j] = array_fill(0, $m, 0.0);

                continue;
            }
            $y = [];
            $column = $a[0];
            $value = $bj[0];
            for ($i = 0; $i < $m; $i++) {
                $y[$i] = $column[$i] * $value;
            }
            for ($p = 1; $p < $k; $p++) {
                $column = $a[$p];
                $value = $bj[$p];
                for ($i = 0; $i < $m; $i++) {
                    $y[$i] = $column[$i] * $value + $y[$i];
                }
            }
            $c[$j] = $y;
        }

        return $c;
    }

    /**
     * A term whose product is exactly zero (one factor zero, the other finite)
     * leaves a fused accumulator unchanged, except that it turns an
     * accumulated -0.0 into +0.0. Such terms are skipped, and an entry whose
     * nonzero terms cancel to zero is recomputed with the full chain, so the
     * result is the kernel's to the bit. Sparse factors (block-diagonal Lambda, the
     * trapezoidal R) then cost only their nonzero terms.
     *
     * @return list<list<float>>
     */
    private static function multiplyFused(array $a, array $b, int $m, int $k, int $n): array
    {
        $finiteColumn = [];
        foreach ($a as $p => $column) {
            $finite = true;
            foreach ($column as $v) {
                $finite = $finite && is_finite($v);
            }
            $finiteColumn[$p] = $finite;
        }
        $c = [];
        for ($j = 0; $j < $n; $j++) {
            $bj = $b[$j];
            $terms = [];
            for ($p = 0; $p < $k; $p++) {
                if ($bj[$p] != 0.0 || $finiteColumn[$p] === false) {
                    $terms[] = $p;
                }
            }
            $y = [];
            for ($i = 0; $i < $m; $i++) {
                $total = 0.0;
                $contributed = false;
                for ($start = 0; $start < $k; $start += self::KC) {
                    $stop = min($k, $start + self::KC);
                    $acc = 0.0;
                    foreach ($terms as $p) {
                        if ($p < $start || $p >= $stop) {
                            continue;
                        }
                        $x = $a[$p][$i];
                        $v = $bj[$p];
                        if ($x == 0.0 && is_finite($v)) {
                            continue;
                        }
                        $contributed = true;
                        $acc = self::fma($x, $v, $acc);
                    }
                    $total = $start === 0 ? $acc : $acc + $total;
                }
                // With no contributing term every step adds a zero product to
                // +0.0, which stays +0.0; only a cancelled sum needs the full chain.
                if ($contributed && $total == 0.0) {
                    $total = self::fusedEntry($a, $bj, $i, $k);
                }
                $y[$i] = $total;
            }
            $c[$j] = $y;
        }

        return $c;
    }

    /** One entry of the fused kernel with every term, zeros included. */
    private static function fusedEntry(array $a, array $bj, int $i, int $k): float
    {
        $total = 0.0;
        for ($start = 0; $start < $k; $start += self::KC) {
            $stop = min($k, $start + self::KC);
            $acc = 0.0;
            for ($p = $start; $p < $stop; $p++) {
                $acc = self::fma($a[$p][$i], $bj[$p], $acc);
            }
            $total = $start === 0 ? $acc : $acc + $total;
        }

        return $total;
    }

    /**
     * A correctly rounded fused multiply-add over all finite operands. Inside
     * Fma::fma's exact range (operands and |c| at most 1e300, |a * b| at least
     * about 2^-850) it is Fma::fma; there its fallbacks would drop the
     * product's rounding error. Outside it: with c = ±0 the result is the
     * rounded product (keeping an underflowed product's sign); with |c| far
     * above the product it is c; otherwise a and b are scaled by powers of two
     * into [1, 2), c alongside, and the result scaled back, settling a
     * subnormal result's second rounding from the exact remainder.
     */
    public static function fma(float $a, float $b, float $c): float
    {
        $product = $a * $b;
        if ($a == 0.0 || $b == 0.0 || is_nan($product) || is_nan($c) || is_infinite($c)
            || is_infinite($a) || is_infinite($b)) {
            return Fma::fma($a, $b, $c);
        }
        $magnitude = abs($product);
        if ($magnitude >= self::TINY_PRODUCT && $magnitude <= self::HUGE && abs($a) <= self::HUGE
            && abs($b) <= self::HUGE && abs($c) <= self::HUGE) {
            return Fma::fma($a, $b, $c);
        }
        if ($c == 0.0) {
            return $product;
        }
        [$aUnit, $aShift] = self::unit($a);
        [$bUnit, $bShift] = self::unit($b);
        $shift = $aShift + $bShift;
        $cExponent = abs($c) < PHP_FLOAT_MIN ? -1023 : self::exponent($c);
        // |a * b| < 2^(2 - shift); below a quarter ulp of c it cannot move c.
        if (2 - $shift <= max($cExponent - 52, -1074) - 2) {
            return $c;
        }
        // A c this far below a * b only breaks ties; a stand-in of the same sign does too.
        $cScaled = $cExponent + $shift < -1000 ? ($c > 0.0 ? 1.0 : -1.0) * self::pow2(-600) : self::scale($c, $shift);
        $result = Fma::fma($aUnit, $bUnit, $cScaled);
        if ($result == 0.0) {
            return $result;
        }
        $exponent = self::exponent($result);
        if ($exponent - $shift >= -1022) {
            return self::scale($result, -$shift);
        }
        // A subnormal result is rounded twice: to 53 bits, then to multiples of
        // 2^-1074. That is wrong only when the first rounding lands exactly on
        // a midpoint of the coarser grid; then the exact remainder decides.
        if ($exponent + 1075 - $shift >= 0) {
            $halves = self::scale($result, 1075 - $shift); // in units of 2^-1075, exact
            if (abs(fmod($halves, 2.0)) == 1.0) {
                [$uh, $ul] = self::twoProduct($aUnit, $bUnit);
                [$th, $tl] = self::twoSum($cScaled, $uh);
                $direction = self::sumSign([-$result, $th, $tl, $ul]);
                if ($direction != 0.0) {
                    return self::scale($halves + $direction, -1075);
                }
            }
        }
        // Keep every step but the last exact.
        $result = self::scale($result, -($exponent + 1000));
        $rest = $shift - $exponent - 1000;
        if ($rest > 1074) {
            return $result * 0.0;
        }

        return $result * self::pow2(-$rest);
    }

    /** @return array{float, float} hi + lo == a * b exactly (Dekker; no overflow or underflow) */
    private static function twoProduct(float $a, float $b): array
    {
        $product = $a * $b;
        $t = 134217729.0 * $a;
        $aHi = $t - ($t - $a);
        $aLo = $a - $aHi;
        $t = 134217729.0 * $b;
        $bHi = $t - ($t - $b);
        $bLo = $b - $bHi;

        return [$product, (($aHi * $bHi - $product) + $aHi * $bLo + $aLo * $bHi) + $aLo * $bLo];
    }

    /** @return array{float, float} hi + lo == a + b exactly (Knuth) */
    private static function twoSum(float $a, float $b): array
    {
        $sum = $a + $b;
        $bp = $sum - $a;

        return [$sum, ($a - ($sum - $bp)) + ($b - $bp)];
    }

    /**
     * The sign (-1, 0 or 1) of an exact sum of floats, via Shewchuk's
     * grow-expansion: the largest nonzero component carries the sign.
     *
     * @param  list<float>  $terms
     */
    private static function sumSign(array $terms): float
    {
        $expansion = [];
        foreach ($terms as $q) {
            $next = [];
            foreach ($expansion as $component) {
                [$q, $error] = self::twoSum($q, $component);
                if ($error != 0.0) {
                    $next[] = $error;
                }
            }
            $next[] = $q;
            $expansion = $next;
        }
        for ($i = count($expansion) - 1; $i >= 0; $i--) {
            if ($expansion[$i] != 0.0) {
                return $expansion[$i] > 0.0 ? 1.0 : -1.0;
            }
        }

        return 0.0;
    }

    /** @return array{float, int} $x * 2^shift in [1, 2) and shift */
    private static function unit(float $x): array
    {
        $shift = 0;
        if (abs($x) < PHP_FLOAT_MIN) {
            $x *= self::pow2(600);
            $shift = 600;
        }
        $exponent = self::exponent($x);

        return [$x * self::pow2(-$exponent), $shift - $exponent];
    }

    /** Unbiased exponent of a normal float. */
    private static function exponent(float $x): int
    {
        return ((unpack('q', pack('e', $x))[1] >> 52) & 0x7FF) - 1023;
    }

    /** 2^n for -1074 <= n <= 1023, built from its bits. */
    private static function pow2(int $n): float
    {
        $bits = $n >= -1022 ? ($n + 1023) << 52 : 1 << ($n + 1074);

        return unpack('e', pack('q', $bits))[1];
    }

    /** $x * 2^n in steps that are exact while no step underflows. */
    private static function scale(float $x, int $n): float
    {
        while ($n > 1000) {
            $x *= self::pow2(1000);
            $n -= 1000;
        }
        while ($n < -1000) {
            $x *= self::pow2(-1000);
            $n += 1000;
        }

        return $x * self::pow2($n);
    }

    /**
     * nalgebra's `a * x` for a vector x: always the gemv path.
     *
     * @param  list<list<float>>  $a
     * @param  list<float>  $x
     * @return list<float>
     */
    public static function multiplyVector(array $a, array $x, ?int $rows = null): array
    {
        return self::multiply($a, [$x], $rows)[0];
    }

    /**
     * nalgebra's `a.tr_mul(b)` (and `ad_mul` for reals): entry (i, j) is the
     * dot product of column i of a with column j of b.
     *
     * @param  list<list<float>>  $a
     * @param  list<list<float>>  $b
     * @return list<list<float>> as columns
     */
    public static function transposeMultiply(array $a, array $b): array
    {
        $c = [];
        foreach ($b as $j => $bj) {
            $column = [];
            foreach ($a as $i => $ai) {
                $column[$i] = self::dot($ai, $bj);
            }
            $c[$j] = $column;
        }

        return $c;
    }

    /**
     * nalgebra's `a.tr_mul(&x)` for a vector x.
     *
     * @param  list<list<float>>  $a
     * @param  list<float>  $x
     * @return list<float>
     */
    public static function transposeMultiplyVector(array $a, array $x): array
    {
        $y = [];
        foreach ($a as $i => $ai) {
            $y[$i] = self::dot($ai, $x);
        }

        return $y;
    }

    /**
     * nalgebra's norm_squared: the sum over columns of each column's dot
     * product with itself.
     *
     * @param  list<list<float>>  $a
     */
    public static function normSquared(array $a): float
    {
        $result = 0.0;
        foreach ($a as $column) {
            $result += self::dot($column, $column);
        }

        return $result;
    }

    /** @param  list<list<float>>  $a */
    public static function norm(array $a): float
    {
        return sqrt(self::normSquared($a));
    }

    /** @param  list<float>  $x */
    public static function vectorNormSquared(array $x): float
    {
        return 0.0 + self::dot($x, $x);
    }

    /** @param  list<float>  $x */
    public static function vectorNorm(array $x): float
    {
        return sqrt(0.0 + self::dot($x, $x));
    }

    /**
     * @param  list<list<float>>  $a  as columns
     * @return list<list<float>> the transpose, as columns
     */
    public static function transpose(array $a, ?int $rows = null): array
    {
        $rows ??= count($a) > 0 ? count($a[0]) : 0;
        $t = [];
        for ($i = 0; $i < $rows; $i++) {
            $t[$i] = array_column($a, $i);
        }

        return $t;
    }

    /**
     * @param  list<list<float>>  $rows
     * @return list<list<float>> the same matrix as columns
     */
    public static function fromRows(array $rows, ?int $columns = null): array
    {
        return self::transpose($rows, $columns);
    }

    /**
     * Rust's f64 signum as num-traits defines it: NaN for NaN, -1 for any value
     * with the sign bit set (including -0.0), else 1.
     */
    public static function signum(float $x): float
    {
        if (is_nan($x)) {
            return NAN;
        }
        if ($x < 0.0 || ($x == 0.0 && fdiv(1.0, $x) < 0.0)) {
            return -1.0;
        }

        return 1.0;
    }

    /**
     * nalgebra's camax for reals: the first entry's magnitude, then simba's
     * simd_max (`a >= b ? a : b`), in column-major order.
     *
     * @param  list<list<float>>  $a
     */
    public static function camax(array $a): float
    {
        $max = null;
        foreach ($a as $column) {
            foreach ($column as $value) {
                $abs = abs($value);
                if ($max === null) {
                    $max = $abs;
                } elseif ($max >= $abs) {
                    continue;
                } else {
                    $max = $abs;
                }
            }
        }

        return $max ?? 0.0;
    }

    /**
     * nalgebra's householder::reflection_axis_mut on a whole vector. Returns the
     * reflected first component and whether a reflection was needed.
     *
     * @param  list<float>  $v
     * @return array{float, bool}
     */
    public static function reflectionAxis(array &$v): array
    {
        $squared = 0.0 + self::dot($v, $v);
        $norm = sqrt($squared);
        $first = $v[0];
        // simba's to_exp: (x, 1) when x >= 0 (so -0.0 keeps sign +1), else (-x, -1).
        if ($first >= 0.0) {
            $modulus = $first;
            $sign = 1.0;
        } else {
            $modulus = -$first;
            $sign = -1.0;
        }
        $signedNorm = $sign * $norm;
        $factor = ($squared + $modulus * $norm) * 2.0;
        $v[0] += $signedNorm;
        if ($factor == 0.0) {
            return [$signedNorm, false];
        }
        $root = sqrt($factor);
        $length = count($v);
        for ($i = 0; $i < $length; $i++) {
            $v[$i] = fdiv($v[$i], $root);
        }
        $norm = sqrt(0.0 + self::dot($v, $v));
        for ($i = 0; $i < $length; $i++) {
            $v[$i] = fdiv($v[$i], $norm);
        }

        return [-$signedNorm, true];
    }
}
