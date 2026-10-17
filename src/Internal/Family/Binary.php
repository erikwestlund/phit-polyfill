<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Family;

use PhpStats\Polyfill\ErrorCode;
use PhpStats\Polyfill\Internal\Glm;
use PhpStats\Polyfill\Internal\Statrs;

/**
 * Shared Bernoulli / grouped-binomial pieces of the logistic and probit fits.
 *
 * @internal
 */
abstract class Binary extends Glm
{
    protected function prepareOutcomes(array $y, array $weights, array $offset, int $n): array
    {
        foreach ($y as $value) {
            if (is_nan($value) || $value < 0.0 || $value > 1.0) {
                throw self::error(ErrorCode::INVALID_INPUT, 'x must be finite and y must be finite proportions in [0,1], with integer successes under the supplied trial weights');
            }
        }
        if (count($weights) !== $n || count($offset) !== $n) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'weights and offset must have one value per row');
        }
        self::checkWeights($weights, true);
        foreach ($y as $i => $value) {
            $events = $weights[$i] * $value;
            if (abs($events - round($events)) > 1e-8) {
                throw self::error(ErrorCode::INVALID_INPUT, 'binary/grouped-binomial outcomes require integer trial counts and successes');
            }
        }
        // Count-rounding tolerance is representation tolerance, not fractional events.
        $allZero = $allOne = true;
        foreach ($y as $i => $value) {
            $value = round($value * $weights[$i]) / $weights[$i];
            $y[$i] = $value;
            $allZero = $allZero && $value == 0.0;
            $allOne = $allOne && $value == 1.0;
        }
        if ($allZero || $allOne) {
            throw self::error(ErrorCode::INVALID_INPUT, $this->name().' regression requires both outcome classes');
        }

        return [$y, $weights, $offset];
    }

    abstract protected function startThreshold(): float;

    abstract protected function link(float $mean): float;

    protected function start(array $y, array $weights, array $offset, array $design, array $scales, int $n, int $p): ?array
    {
        $extreme = false;
        foreach ($offset as $o) {
            if (abs($o) > $this->startThreshold()) {
                $extreme = true;
                break;
            }
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
            $num = 0.0;
            foreach ($y as $i => $value) {
                $num += $value * $weights[$i];
            }
            $mean = $num / array_sum($weights);
            $center = 0.0;
            foreach ($offset as $o) {
                $center += $o / $n;
            }
            $theta = array_fill(0, $p, 0.0);
            $theta[$j] = ($this->link($mean) - $center) * $scales[$j];

            return $theta;
        }

        return null;
    }

    /**
     * Whether eta itself certifies complete or quasi-complete separation: every
     * signed margin nonnegative and at least one positive (not an LP detector).
     */
    protected function separated(array $linear, array $y): bool
    {
        foreach ($linear as $i => $e) {
            if ($y[$i] > 0.0 && $y[$i] < 1.0 && $e != 0.0) {
                return false;
            }
        }
        $positive = false;
        foreach ($linear as $i => $e) {
            $margin = (2.0 * $y[$i] - 1.0) * $e;
            if ($margin < 0.0) {
                return false;
            }
            $positive = $positive || $margin > 0.0;
        }

        return $positive;
    }

    protected function rankMessage(): string
    {
        return 'design/information is rank deficient; separation may prevent a finite MLE';
    }

    protected function overflowMessage(): string
    {
        return $this->name().' calculation overflowed; separation or extreme scaling may be present';
    }

    protected function nonconvergenceMessage(int $max): string
    {
        return $this->name()." regression did not converge in $max iterations; separation may be present";
    }

    protected function statistics(
        array $eta, array $y, array $weights, array $offset, float $likelihood,
        bool $intercept, bool $computeNull, int $n, int $p, int $maxIterations,
    ): array {
        $constant = 0.0;
        foreach ($y as $i => $v) {
            if ($v == 0.0 || $v == 1.0) {
                continue;
            }
            $w = $weights[$i];
            $successes = $w * $v;
            $constant += Statrs::lnGamma($w + 1.0) - Statrs::lnGamma($successes + 1.0) - Statrs::lnGamma($w - $successes + 1.0);
        }
        if (! $computeNull) {
            $null = $likelihood;
        } elseif ($intercept) {
            // With offsets, the null intercept needs its own one-coefficient fit.
            $null = $this->fitData(array_fill(0, $n, [1.0]), $y, $weights, $offset, $maxIterations, false)['log_likelihood'] - $constant;
        } else {
            $null = $this->likelihood($offset, $y, $weights);
        }
        $saturated = 0.0;
        foreach ($y as $i => $v) {
            $logp = $v == 0.0 ? 0.0 : $v * log($v);
            $logq = $v == 1.0 ? 0.0 : (1.0 - $v) * log1p(-$v);
            $saturated += $weights[$i] * ($logp + $logq);
        }
        foreach ([$constant, $likelihood + $constant, $saturated, $null] as $value) {
            if (! is_finite($value)) {
                throw self::error(ErrorCode::NONFINITE, 'binomial likelihood statistics exceed floating-point range');
            }
        }
        $gap = 64.0 * PHP_FLOAT_EPSILON * (1.0 + abs($saturated) + abs($likelihood));
        if ($saturated - $likelihood < -$gap || $saturated - $null < -$gap) {
            throw self::error(ErrorCode::INTERNAL, 'binomial likelihood exceeds saturated likelihood beyond rounding precision');
        }
        $aic = 2.0 * $p - 2.0 * ($likelihood + $constant);
        foreach ([2.0 * ($saturated - $likelihood), 2.0 * ($saturated - $null), $aic] as $value) {
            if (! is_finite($value)) {
                throw self::error(ErrorCode::NONFINITE, 'binomial deviance or AIC exceeds floating-point range');
            }
        }

        return [
            'log_likelihood' => $likelihood + $constant,
            'deviance' => max(2.0 * ($saturated - $likelihood), 0.0),
            'null_deviance' => max(2.0 * ($saturated - $null), 0.0),
            'aic' => $aic,
        ];
    }
}
