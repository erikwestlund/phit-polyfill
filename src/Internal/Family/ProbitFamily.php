<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Family;

use PhpStats\Polyfill\ErrorCode;
use PhpStats\Polyfill\Internal\Statrs;

/** @internal */
final class ProbitFamily extends Binary
{
    protected function name(): string
    {
        return 'probit';
    }

    protected function solver(): string
    {
        return 'fisher_scoring';
    }

    protected function startThreshold(): float
    {
        return 8.0;
    }

    protected function link(float $mean): float
    {
        return Statrs::normalQuantile($mean);
    }

    /** ln Phi(x); a Mills-ratio continued fraction beyond -8 keeps tails finite. */
    private static function logProbability(float $value): float
    {
        if ($value < -8.0) {
            $x = -$value;
            $denominator = $x;
            for ($k = 100; $k >= 1; $k--) {
                $denominator = $x + $k / $denominator;
            }

            return -0.5 * $x * $x - 0.5 * log(2.0 * M_PI) - log($denominator);
        }

        return log(0.5 * Statrs::erfc(-$value / M_SQRT2));
    }

    private static function logDensity(float $value): float
    {
        return -0.5 * $value * $value - 0.5 * log(2.0 * M_PI);
    }

    private static function rootWeight(float $value): float
    {
        return exp(self::logDensity($value) - 0.5 * (self::logProbability($value) + self::logProbability(-$value)));
    }

    protected function likelihood(array $eta, array $y, array $weights): float
    {
        // Neumaier-compensated sum, as the extension.
        $sum = 0.0;
        $compensation = 0.0;
        foreach ($eta as $i => $e) {
            $v = $y[$i];
            $term = $weights[$i] * ($v * self::logProbability($e) + (1.0 - $v) * self::logProbability(-$e));
            $total = $sum + $term;
            $compensation += abs($sum) >= abs($term) ? ($sum - $total) + $term : ($term - $total) + $sum;
            $sum = $total;
        }

        return $sum + $compensation;
    }

    protected function rootWeights(array $eta): array
    {
        // Fisher scoring uses expected information, matching R glm(probit).
        $out = [];
        foreach ($eta as $i => $e) {
            $w = self::rootWeight($e);
            if (! is_finite($w)) {
                throw self::error(ErrorCode::NONFINITE, 'probit information exceeds floating-point range');
            }
            $out[$i] = $w;
        }

        return $out;
    }

    protected function finalRootWeights(array $eta): array
    {
        $out = [];
        foreach ($eta as $i => $e) {
            $w = self::rootWeight($e);
            if (! is_finite($w)) {
                throw self::error(ErrorCode::NONFINITE, 'final probit information exceeds floating-point range');
            }
            $out[$i] = $w;
        }

        return $out;
    }

    protected function scores(array $eta, array $y): array
    {
        $out = [];
        foreach ($eta as $i => $e) {
            $ld = self::logDensity($e);
            $out[$i] = $y[$i] * exp($ld - self::logProbability($e)) - (1.0 - $y[$i]) * exp($ld - self::logProbability(-$e));
        }

        return $out;
    }

    protected function tolerances(): array
    {
        return [1e-11, 1e-13];
    }

    protected function fittedValues(array $eta): array
    {
        $out = [];
        foreach ($eta as $i => $e) {
            $out[$i] = 0.5 * Statrs::erfc(-$e / M_SQRT2);
        }

        return $out;
    }
}
