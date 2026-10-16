<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/**
 * A correctly rounded fused multiply-add, which PHP lacks. The extension's
 * matrix products on aarch64 run matrixmultiply's NEON kernel, which
 * accumulates every product with a fused multiply-add.
 *
 * The emulation follows Boldo and Melquiond ("Emulation of FMA and correctly
 * rounded sums: proved algorithms using rounding to odd", 2008): the exact
 * product from Dekker's algorithm and the exact sum from Knuth's TwoSum are
 * combined with one addition rounded to odd. It is exact unless a product
 * or sum overflows or an error term falls into the subnormal range; those
 * operands fall back to the unfused expression.
 *
 * @internal
 */
final class Fma
{
    private const float SPLITTER = 134217729.0; // 2^27 + 1

    /** Products above this magnitude may overflow in Dekker's split. */
    private const float SPLIT_LIMIT = 1e300;

    /** Error terms below this magnitude may have lost bits to subnormals. */
    private const float TINY = 1e-290;

    public static function fma(float $a, float $b, float $c): float
    {
        if (is_infinite($c) && is_finite($a) && is_finite($b)) {
            return $c; // the exact product is finite even when a * b overflows
        }
        $product = $a * $b;
        if (is_nan($product) || is_infinite($product) || is_infinite($c) || is_nan($c)
            || abs($a) > self::SPLIT_LIMIT || abs($b) > self::SPLIT_LIMIT) {
            return $product + $c;
        }
        // Dekker: $product + $error == $a * $b exactly.
        $t = self::SPLITTER * $a;
        $aHi = $t - ($t - $a);
        $aLo = $a - $aHi;
        $t = self::SPLITTER * $b;
        $bHi = $t - ($t - $b);
        $bLo = $b - $bHi;
        $error = (($aHi * $bHi - $product) + $aHi * $bLo + $aLo * $bHi) + $aLo * $bLo;
        if ($error == 0.0) {
            return $product + $c;
        }
        if (abs($error) < self::TINY) {
            return $product + $c;
        }
        // TwoSum: $high + $low == $c + $product exactly.
        $high = $c + $product;
        $bp = $high - $c;
        $low = ($c - ($high - $bp)) + ($product - $bp);

        return $high + self::addOdd($low, $error);
    }

    /** $x + $y rounded to odd: exact sums are kept, inexact ones get an odd last bit. */
    private static function addOdd(float $x, float $y): float
    {
        $sum = $x + $y;
        $bp = $sum - $x;
        $error = ($x - ($sum - $bp)) + ($y - $bp);
        if ($error == 0.0) {
            return $sum;
        }
        $bits = unpack('q', pack('e', $sum))[1];
        if (($bits & 1) === 1) {
            return $sum;
        }
        // Step one unit in the last place toward the exact value.
        $towardLarger = ($error > 0.0) === ($sum > 0.0);
        $bits += $towardLarger ? 1 : -1;

        return unpack('e', pack('q', $bits))[1];
    }
}
