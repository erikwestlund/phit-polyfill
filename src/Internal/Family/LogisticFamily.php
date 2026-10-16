<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Family;

/** @internal */
final class LogisticFamily extends Binary
{
    protected function name(): string
    {
        return 'logistic';
    }

    protected function solver(): string
    {
        return 'newton_irls';
    }

    protected function startThreshold(): float
    {
        return 20.0;
    }

    protected function link(float $mean): float
    {
        return log($mean) - log1p(-$mean);
    }

    /** Branching keeps exp() from overflowing for large logits. */
    private static function sigmoid(float $value): float
    {
        if ($value >= 0.0) {
            return 1.0 / (1.0 + exp(-$value));
        }
        $e = exp($value);

        return $e / (1.0 + $e);
    }

    protected function likelihood(array $eta, array $y, array $weights): float
    {
        // Neumaier-compensated, as the extension: near the optimum, rounding in a plain
        // sum can make an exact Newton step look like a decrease and stall the search.
        $sum = 0.0;
        $compensation = 0.0;
        foreach ($eta as $i => $e) {
            $v = $y[$i];
            $term = $weights[$i] * (-$v * self::softplus(-$e) - (1.0 - $v) * self::softplus($e));
            $total = $sum + $term;
            $compensation += abs($sum) >= abs($term) ? ($sum - $total) + $term : ($term - $total) + $sum;
            $sum = $total;
        }

        return $sum + $compensation;
    }

    protected function rootWeights(array $eta): array
    {
        // sqrt(p(1-p)) without subtracting a probability rounded to 1.
        $out = [];
        foreach ($eta as $i => $e) {
            $x = exp(-abs($e));
            $out[$i] = sqrt($x) / (1.0 + $x);
        }

        return $out;
    }

    protected function finalRootWeights(array $eta): array
    {
        return $this->rootWeights($eta);
    }

    protected function scores(array $eta, array $y): array
    {
        $out = [];
        foreach ($eta as $i => $e) {
            $out[$i] = $y[$i] * self::sigmoid(-$e) - (1.0 - $y[$i]) * self::sigmoid($e);
        }

        return $out;
    }

    protected function tolerances(): array
    {
        return [1e-8, 1e-10];
    }

    protected function fittedValues(array $eta): array
    {
        return array_map(self::sigmoid(...), $eta);
    }
}
