<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

use PhpStats\Polyfill\ErrorCode;

/**
 * Input/result policy shared by the scalar and batch distribution functions
 * (stats-core distributions/mod.rs and the extension's PHP glue).
 *
 * @internal
 */
final class Support
{
    public static function finite(string $name, float $value): float
    {
        if (! is_finite($value)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, "$name must be finite");
        }

        return $value;
    }

    public static function positive(string $name, float $value): float
    {
        if (self::finite($name, $value) <= 0.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "$name must be positive");
        }

        return $value;
    }

    public static function invalid(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    /**
     * The extension turns a Rust panic (an internal cdflib invariant, such as
     * the indeterminate incomplete-gamma regime) into a PHP \Error with code 0.
     */
    public static function panic(string $message): never
    {
        throw new RustPanic('Rust panic: '.$message);
    }

    /**
     * cdflib's TOMS 708/654 routines never return for some extreme arguments
     * and panic in the deep asymptotic regime of the incomplete gamma. As in
     * the extension (stats-core distributions/mod.rs), those regimes return
     * NaN up front and any remaining panic becomes NaN, which the distribution
     * functions report as a numerical failure (1006).
     */
    private const float CDFLIB_ARGUMENT_LIMIT = 1e300;

    /** @return array{0: float, 1: float} */
    public static function betaInc(float $a, float $b, float $x, float $y): array
    {
        $larger = is_nan($a) ? $b : (is_nan($b) ? $a : max($a, $b));
        $smaller = is_nan($x) ? $y : (is_nan($y) ? $x : min($x, $y));
        if (! is_finite($a + $b) || ($larger > self::CDFLIB_ARGUMENT_LIMIT && $smaller < 1.0 / self::CDFLIB_ARGUMENT_LIMIT)) {
            return [NAN, NAN];
        }
        try {
            return Cdflib::betaInc($a, $b, $x, $y);
        } catch (RustPanic) {
            return [NAN, NAN];
        }
    }

    /** @return array{0: float, 1: float} */
    public static function gammaInc(float $a, float $x): array
    {
        if (is_nan($a) || is_nan($x) || $a > self::CDFLIB_ARGUMENT_LIMIT || $x > self::CDFLIB_ARGUMENT_LIMIT) {
            return [NAN, NAN];
        }
        try {
            return Cdflib::gammaInc($a, $x);
        } catch (RustPanic) {
            return [NAN, NAN];
        }
    }

    /** Evaluate a density, CDF or survival function at x (±∞ allowed, NaN rejected). */
    public static function point(float $x, \Closure $f): float
    {
        if (is_nan($x)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'evaluation points must not be NaN');
        }

        return self::result($f($x));
    }

    /** Evaluate a quantile function at a probability in [0, 1]. */
    public static function probability(float $p, \Closure $f): float
    {
        if (is_nan($p)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'probabilities must not be NaN');
        }
        if (is_nan($p) || $p < 0.0 || $p > 1.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'probabilities must be in [0, 1]');
        }

        return self::result($f($p));
    }

    /** @return list<float> */
    public static function points(array $values, \Closure $f): array
    {
        $out = [];
        foreach ($values as $value) {
            $out[] = self::point(self::element($value), $f);
        }

        return $out;
    }

    /** @return list<float> */
    public static function probabilities(array $values, \Closure $f): array
    {
        $out = [];
        foreach ($values as $value) {
            $out[] = self::probability(self::element($value), $f);
        }

        return $out;
    }

    /** Array elements must be PHP ints or floats; there is no missing-value policy. */
    private static function element(mixed $value): float
    {
        if (is_float($value)) {
            return $value;
        }
        if (is_int($value)) {
            return (float) $value;
        }
        throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'array elements must be integers or floats');
    }

    private static function result(float $value): float
    {
        if (is_nan($value)) {
            throw ErrorCode::exception(ErrorCode::NONCONVERGENCE, 'numerical evaluation failed for these arguments');
        }

        return $value;
    }

    /** Rust f64::round: half away from zero, exact for every double. */
    public static function round(float $v): float
    {
        if (! is_finite($v)) {
            return $v;
        }
        $f = floor(abs($v));
        $r = abs($v) - $f >= 0.5 ? $f + 1.0 : $f;

        return $v < 0.0 ? -$r : $r;
    }
}
