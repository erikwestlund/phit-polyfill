<?php

declare(strict_types=1);

namespace PhpStats\Polyfill;

use PhpStats\Polyfill\Internal\Family\ProbitFamily;

/**
 * Binary or grouped-binomial probit regression by Fisher scoring (stats_probit).
 */
final class Probit
{
    public static function fit(
        array $x,
        array $y,
        string $missing = 'drop',
        ?array $weights = null,
        ?array $offset = null,
        bool $intercept = true,
    ): array {
        return new ProbitFamily()->run($x, $y, $missing, $weights, $offset, $intercept, 'binary', 'stats_probit');
    }
}
