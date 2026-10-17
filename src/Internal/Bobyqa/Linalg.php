<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * The crate's float helpers and its linalg.rs subset, loop for loop. Matrices
 * are flat column-major lists: element (i, j) of an nr-row matrix is
 * $a[$i + $j * $nr], the layout of the crate's Mat.
 *
 * @internal
 */
final class Linalg
{
    public const float EPS = PHP_FLOAT_EPSILON;

    public const float REALMIN = PHP_FLOAT_MIN;

    public const float REALMAX = PHP_FLOAT_MAX;

    /** 2^1022, the norm rescue's upper scale. */
    private const float SCALMAX = 2.0 ** 1022;

    /**
     * Rust's f64::max: a NaN operand loses, and on aarch64 (fmaxnm) +0 beats -0.
     */
    public static function max(float $a, float $b): float
    {
        if (\is_nan($a)) {
            return $b;
        }
        if (\is_nan($b)) {
            return $a;
        }
        if ($a > $b) {
            return $a;
        }
        if ($b > $a) {
            return $b;
        }

        // Equal: only signed zeros can differ, and their sum is the larger.
        return $a == 0.0 ? $a + $b : $a;
    }

    /** Rust's f64::min: a NaN operand loses, and -0 beats +0. */
    public static function min(float $a, float $b): float
    {
        if (\is_nan($a)) {
            return $b;
        }
        if (\is_nan($b)) {
            return $a;
        }
        if ($a < $b) {
            return $a;
        }
        if ($b < $a) {
            return $b;
        }

        return $a == 0.0 ? -(-$a - $b) : $a;
    }

    /** |a| with the sign bit of b. */
    public static function copysign(float $a, float $b): float
    {
        $negative = (\ord(\pack('E', $b)[0]) & 0x80) !== 0;

        return $negative ? -\abs($a) : \abs($a);
    }

    /** Rust's f64::midpoint. */
    public static function midpoint(float $a, float $b): float
    {
        $lo = 2.0 * self::REALMIN;
        $hi = self::REALMAX / 2.0;
        $absA = \abs($a);
        $absB = \abs($b);
        if ($absA <= $hi && $absB <= $hi) {
            return ($a + $b) / 2.0;
        }
        if ($absA < $lo) {
            return $a + ($b / 2.0);
        }
        if ($absB < $lo) {
            return ($a / 2.0) + $b;
        }

        return ($a / 2.0) + ($b / 2.0);
    }

    /**
     * Rust's Iterator::sum for f64, which starts from -0.0.
     *
     * @param  list<float>  $x
     */
    public static function sum(array $x): float
    {
        $s = -0.0;
        foreach ($x as $v) {
            $s += $v;
        }

        return $s;
    }

    /**
     * Rust's sum of |x|, from -0.0.
     *
     * @param  list<float>  $x
     */
    public static function sumAbs(array $x): float
    {
        $s = -0.0;
        foreach ($x as $v) {
            $s += \abs($v);
        }

        return $s;
    }

    /** @param  list<float>  $x */
    public static function allFinite(array $x): bool
    {
        foreach ($x as $v) {
            if (! \is_finite($v)) {
                return false;
            }
        }

        return true;
    }

    /**
     * @param  list<float>  $x
     * @param  list<float>  $y
     */
    public static function inprod(array $x, array $y): float
    {
        $z = 0.0;
        $count = count($x);
        for ($i = 0; $i < $count; $i++) {
            $z += $x[$i] * $y[$i];
        }

        return $z;
    }

    /**
     * Row vector x times the nr x nc matrix y.
     *
     * @param  list<float>  $x
     * @param  list<float>  $y
     * @return list<float>
     */
    public static function matprod12(array $x, array $y, int $nr, int $nc): array
    {
        $z = [];
        for ($j = 0; $j < $nc; $j++) {
            $s = 0.0;
            $offset = $j * $nr;
            for ($i = 0; $i < $nr; $i++) {
                $s += $x[$i] * $y[$offset + $i];
            }
            $z[] = $s;
        }

        return $z;
    }

    /**
     * The nr x nc matrix x times column vector y, accumulated column by column.
     *
     * @param  list<float>  $x
     * @param  list<float>  $y
     * @return list<float>
     */
    public static function matprod21(array $x, int $nr, int $nc, array $y): array
    {
        if ($nc === 0) {
            return \array_fill(0, $nr, 0.0);
        }
        $z = [];
        $y0 = $y[0];
        for ($i = 0; $i < $nr; $i++) {
            // 0.0 + p, not p: the crate keeps the zero-start addition (it maps -0.0 to +0.0).
            $z[] = 0.0 + $x[$i] * $y0;
        }
        for ($j = 1; $j < $nc; $j++) {
            $yj = $y[$j];
            $offset = $j * $nr;
            for ($i = 0; $i < $nr; $i++) {
                $z[$i] += $x[$offset + $i] * $yj;
            }
        }

        return $z;
    }

    /**
     * The xr x xc matrix x times the xc x yc matrix y.
     *
     * @param  list<float>  $x
     * @param  list<float>  $y
     * @return list<float>
     */
    public static function matprod22(array $x, int $xr, int $xc, array $y, int $yc): array
    {
        $z = \array_fill(0, $xr * $yc, 0.0);
        for ($j = 0; $j < $yc; $j++) {
            $zOffset = $j * $xr;
            for ($i = 0; $i < $xc; $i++) {
                $yij = $y[$i + $j * $xc];
                $xOffset = $i * $xr;
                for ($r = 0; $r < $xr; $r++) {
                    $z[$zOffset + $r] += $x[$xOffset + $r] * $yij;
                }
            }
        }

        return $z;
    }

    /**
     * x * y^T as a count(x) x count(y) matrix.
     *
     * @param  list<float>  $x
     * @param  list<float>  $y
     * @return list<float>
     */
    public static function outprod(array $x, array $y): array
    {
        $z = [];
        foreach ($y as $yi) {
            foreach ($x as $xr) {
                $z[] = $xr * $yi;
            }
        }

        return $z;
    }

    /**
     * PRIMA's p_norm with p = 2, including its overflow/underflow rescue.
     *
     * @param  list<float>  $x
     */
    public static function norm(array $x): float
    {
        $maxabs = 0.0;
        foreach ($x as $v) {
            $maxabs = self::max($maxabs, \abs($v));
        }
        if ($x === []) {
            return 0.0;
        }
        if (! self::allFinite($x)) {
            $y = 0.0;
            foreach ($x as $v) {
                $y += \abs($v);
            }

            return $y;
        }
        if ($maxabs <= 0.0) {
            return 0.0;
        }
        $sumsq = 0.0;
        foreach ($x as $v) {
            $sumsq += $v * $v;
        }
        $y = \sqrt($sumsq);
        if ((\is_infinite($y) && $y > 0.0) || $y <= 0.0) {
            $scaling = self::min(self::SCALMAX, self::max(self::REALMIN, $maxabs));
            $scaled = 0.0;
            foreach ($x as $v) {
                $s = $v / $scaling;
                $scaled += $s * $s;
            }
            $y = $scaling * \sqrt($scaled);
        }

        return $y;
    }

    /**
     * PRIMA's planerot: the Givens rotation [[c, s], [-s, c]] zeroing x1, returned as [c, s].
     *
     * @return array{float, float}
     */
    public static function planerot(float $x0, float $x1): array
    {
        if (\is_nan($x0) || \is_nan($x1)) {
            return [1.0, 0.0];
        }
        if (\is_infinite($x0) && \is_infinite($x1)) {
            return [self::copysign(1.0 / \sqrt(2.0), $x0), self::copysign(1.0 / \sqrt(2.0), $x1)];
        }
        if (\abs($x0) <= 0.0 && \abs($x1) <= 0.0) {
            return [1.0, 0.0];
        }
        if (\abs($x1) <= self::EPS * \abs($x0)) {
            return [self::copysign(1.0, $x0), 0.0];
        }
        if (\abs($x0) <= self::EPS * \abs($x1)) {
            return [0.0, self::copysign(1.0, $x1)];
        }
        $lo = \sqrt(self::REALMIN);
        $hi = \sqrt(self::REALMAX / 2.1);
        if (\abs($x0) > $lo && \abs($x0) < $hi && \abs($x1) > $lo && \abs($x1) < $hi) {
            $r = self::norm([$x0, $x1]);

            return [$x0 / $r, $x1 / $r];
        }
        if (\abs($x0) > \abs($x1)) {
            $t = $x1 / $x0;
            $u = self::copysign(self::max(self::max(1.0, \abs($t)), \sqrt(1.0 + $t * $t)), $x0);

            return [1.0 / $u, $t / $u];
        }
        $t = $x0 / $x1;
        $u = self::copysign(self::max(self::max(1.0, \abs($t)), \sqrt(1.0 + $t * $t)), $x1);

        return [$t / $u, 1.0 / $u];
    }

    /**
     * A += alpha * x * x^T on the lower triangle, then symmetrize (n x n).
     *
     * @param  list<float>  $a
     * @param  list<float>  $x
     */
    public static function r1update(array &$a, int $n, float $alpha, array $x): void
    {
        for ($j = 0; $j < $n; $j++) {
            $xj = $x[$j];
            for ($i = $j; $i < $n; $i++) {
                $a[$i + $j * $n] += $alpha * $x[$i] * $xj;
            }
        }
        self::symmetrize($a, $n);
    }

    /**
     * A += alpha * (x * y^T + y * x^T) on the lower triangle, evaluated
     * (A + t1) + t2, then symmetrize.
     *
     * @param  list<float>  $a
     * @param  list<float>  $x
     * @param  list<float>  $y
     */
    public static function r2update(array &$a, int $n, float $alpha, array $x, array $y): void
    {
        for ($j = 0; $j < $n; $j++) {
            $xj = $x[$j];
            $yj = $y[$j];
            for ($i = $j; $i < $n; $i++) {
                $k = $i + $j * $n;
                $a[$k] = ($a[$k] + $alpha * $x[$i] * $yj) + $alpha * $y[$i] * $xj;
            }
        }
        self::symmetrize($a, $n);
    }

    /** @param  list<float>  $a */
    public static function symmetrize(array &$a, int $n): void
    {
        for ($j = 0; $j < $n; $j++) {
            for ($i = 0; $i < $j; $i++) {
                $a[$i + $j * $n] = $a[$j + $i * $n];
            }
        }
    }

    /**
     * Column j of an nr-row matrix.
     *
     * @param  list<float>  $a
     * @return list<float>
     */
    public static function col(array $a, int $nr, int $j): array
    {
        return \array_slice($a, $j * $nr, $nr);
    }
}
