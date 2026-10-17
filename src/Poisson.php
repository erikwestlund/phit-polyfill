<?php

declare(strict_types=1);

namespace PhpStats\Polyfill;

use PhpStats\Polyfill\Internal\Family\PoissonFamily;

/**
 * Poisson count regression with exposure offsets by damped Newton IRLS (stats_poisson).
 */
final class Poisson
{
    public static function fit(
        array $x,
        array $y,
        string $missing = 'drop',
        ?array $weights = null,
        ?array $offset = null,
        bool $intercept = true,
    ): array {
        return (new PoissonFamily())->run($x, $y, $missing, $weights, $offset, $intercept, 'count', 'stats_poisson');
    }
}
