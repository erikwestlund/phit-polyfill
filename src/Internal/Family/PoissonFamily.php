<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Family;

use PhpStats\Polyfill\ErrorCode;
use PhpStats\Polyfill\Internal\Glm;
use PhpStats\Polyfill\Internal\Statrs;

/** @internal */
final class PoissonFamily extends Glm
{
    protected function name(): string
    {
        return 'poisson';
    }

    protected function solver(): string
    {
        return 'newton_irls';
    }

    protected function prepareOutcomes(array $y, array $weights, array $offset, int $n): array
    {
        $allZero = true;
        foreach ($y as $value) {
            if ($value < 0.0 || fmod($value, 1.0) != 0.0) {
                throw self::error(ErrorCode::INVALID_INPUT, 'x must be finite and y must contain nonnegative integer counts');
            }
            $allZero = $allZero && $value == 0.0;
        }
        if ($allZero) {
            throw self::error(ErrorCode::INVALID_INPUT, 'poisson regression requires at least one positive count');
        }
        if (count($weights) !== $n || count($offset) !== $n) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'weights and offset must have one value per row');
        }
        self::checkWeights($weights, false);

        return [$y, $weights, $offset];
    }

    /** Log-sum-exp keeps exposures finite even for offsets beyond exp()'s range. */
    private static function logSumExp(array $values): float
    {
        $max = -INF;
        foreach ($values as $v) {
            $max = max($max, $v);
        }
        $sum = 0.0;
        foreach ($values as $v) {
            $sum += exp($v - $max);
        }

        return $max + log($sum);
    }

    private static function nullIntercept(array $y, array $weights, array $offset): float
    {
        $events = [];
        foreach ($y as $i => $v) {
            if ($v > 0.0) {
                $events[] = log($v) + log($weights[$i]);
            }
        }
        $exposure = [];
        foreach ($offset as $i => $o) {
            $exposure[] = $o + log($weights[$i]);
        }

        return self::logSumExp($events) - self::logSumExp($exposure);
    }

    protected function start(array $y, array $weights, array $offset, array $design, array $scales, int $n, int $p): ?array
    {
        $extreme = false;
        foreach ($offset as $o) {
            $extreme = $extreme || abs($o) > 20.0;
        }
        foreach ($y as $v) {
            $extreme = $extreme || $v > 1e8;
        }
        if (! $extreme) {
            return null;
        }
        foreach ($design as $j => $column) {
            foreach ($column as $value) {
                if ($value != 1.0) {
                    continue 2;
                }
            }
            $theta = array_fill(0, $p, 0.0);
            $theta[$j] = self::nullIntercept($y, $weights, $offset) * $scales[$j];

            return $theta;
        }

        return null;
    }

    protected function likelihood(array $eta, array $y, array $weights): float
    {
        $sum = 0.0;
        $half = 0.5;
        $ln2pi = log(2.0 * M_PI);
        foreach ($eta as $i => $e) {
            $v = $y[$i];
            if ($v >= 16.0) {
                // Deviance + Stirling error avoids subtracting O(y log y) terms.
                $d = $e - log($v);
                $inv = 1.0 / $v;
                $inv2 = $inv * $inv;
                $correction = $inv * (1.0 / 12.0 + $inv2 * (-1.0 / 360.0 + $inv2 * (1.0 / 1260.0 + $inv2 * (-1.0 / 1680.0 + $inv2 / 1188.0))));
                $term = -$v * (expm1($d) - $d) - $half * ($ln2pi + log($v)) - $correction;
            } else {
                $term = $v * $e - exp($e) - Statrs::lnGamma($v + 1.0);
            }
            $sum += $weights[$i] * $term;
        }

        return $sum;
    }

    private static function deviance(array $eta, array $y, array $weights): float
    {
        $sum = 0.0;
        foreach ($eta as $i => $e) {
            $v = $y[$i];
            if ($v == 0.0) {
                $term = exp($e);
            } else {
                $d = $e - log($v);
                $term = $v * (expm1($d) - $d);
            }
            $sum += $weights[$i] * $term;
        }

        return 2.0 * $sum;
    }

    protected function rootWeights(array $eta): array
    {
        $out = [];
        foreach ($eta as $i => $e) {
            $w = exp(0.5 * $e);
            if (! is_finite($w)) {
                throw self::error(ErrorCode::NONFINITE, 'Poisson information exceeds floating-point range');
            }
            $out[$i] = $w;
        }

        return $out;
    }

    protected function finalRootWeights(array $eta): array
    {
        $out = [];
        foreach ($eta as $i => $e) {
            $out[$i] = exp(0.5 * $e);
        }

        return $out;
    }

    protected function scores(array $eta, array $y): array
    {
        $out = [];
        foreach ($eta as $i => $e) {
            $out[$i] = $y[$i] - exp($e);
        }

        return $out;
    }

    protected function rankMessage(): string
    {
        return 'design/information is rank deficient; a boundary MLE may prevent a finite fit';
    }

    protected function overflowMessage(): string
    {
        return 'poisson calculation overflowed; a boundary MLE or extreme scaling may be present';
    }

    protected function nonconvergenceMessage(int $max): string
    {
        return "poisson regression did not converge in $max iterations; a boundary MLE or extreme scaling may be present";
    }

    protected function tolerances(): array
    {
        return [1e-8, 1e-10];
    }

    protected function fittedValues(array $eta): array
    {
        return array_map(exp(...), $eta);
    }

    protected function statistics(
        array $eta, array $y, array $weights, array $offset, float $likelihood,
        bool $intercept, bool $computeNull, int $n, int $p, int $maxIterations,
    ): array {
        $nullEta = $offset;
        if ($intercept && $computeNull) {
            $logMean = self::nullIntercept($y, $weights, $offset);
            foreach ($nullEta as $i => $o) {
                $nullEta[$i] = $o + $logMean;
            }
        }
        $deviance = self::deviance($eta, $y, $weights);
        $nullDeviance = self::deviance($nullEta, $y, $weights);
        $aic = 2.0 * $p - 2.0 * $likelihood;
        $finite = is_finite($deviance) && is_finite($nullDeviance) && is_finite($aic);
        foreach ($eta as $e) {
            $finite = $finite && is_finite(exp($e));
        }
        if (! $finite) {
            throw self::error(ErrorCode::NONFINITE, 'Poisson statistics exceed floating-point range');
        }

        return [
            'log_likelihood' => $likelihood,
            'deviance' => $deviance,
            'null_deviance' => $nullDeviance,
            'aic' => $aic,
        ];
    }
}
