<?php

declare(strict_types=1);

use Phit\Polyfill\Polyfill;

// Extreme-parameter grid, after the extension's fuzzing of cdflib (phit
// bd98bf7): regimes where cdflib hangs or panics must fail fast with 1006 in
// both implementations, and everything else must agree exactly.

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

const EXTREME_VALUES = [1e-300, 1e-10, 0.5, 1.0, 3.0, 1e10, 1e30, 1e100, 1e299, 1e301, 1e308, PHP_FLOAT_MAX];

/** Value or exception code; Rust panics surface as \Error with code 0 in both. */
function extremeOutcome(callable $f, array $arguments): mixed
{
    try {
        return $f(...$arguments);
    } catch (\Throwable $e) {
        return 'E'.$e->getCode();
    }
}

it('agrees with the extension on extreme parameters', function (string $family, array $operations, int $parameters) {
    $checked = 0;
    foreach ($operations as $operation) {
        $fn = "Phit\\dist_{$family}_$operation";
        $points = $operation === 'quantile' ? [1e-300, 1e-10, 0.3, 0.5, 0.9, 1.0 - 1e-12] : [0.0, 1e-300, 0.5, 2.0, 1e30, 1e300];
        if ($family === 'beta' && $operation !== 'quantile') {
            $points = [1e-300, 1e-10, 0.5, 1.0 - 1e-12];
        }
        $grids = $parameters === 1
            ? array_map(fn ($a) => [$a], EXTREME_VALUES)
            : array_merge(...array_map(fn ($a) => array_map(fn ($b) => [$a, $b], EXTREME_VALUES), EXTREME_VALUES));
        foreach ($grids as $params) {
            foreach ($points as $x) {
                $arguments = [$x, ...$params];
                $polyfill = extremeOutcome(Polyfill::implementation($fn), $arguments);
                $native = extremeOutcome($fn, $arguments);
                parityLog('distributions-extreme', $fn, $polyfill, $native, is_string($polyfill) || is_string($native));
                $same = $polyfill === $native || (is_float($polyfill) && is_float($native) && is_nan($polyfill) && is_nan($native));
                expect($same)->toBeTrue(sprintf('%s(%s): polyfill %s, extension %s', $fn, implode(', ', array_map(fn ($v) => var_export($v, true), $arguments)), var_export($polyfill, true), var_export($native, true)));
                $checked++;
            }
        }
    }
    expect($checked)->toBeGreaterThan(0);
})->with([
    'F' => ['f', ['cdf', 'sf', 'logcdf', 'logsf', 'quantile'], 2],
    'beta' => ['beta', ['cdf', 'sf', 'logcdf', 'quantile'], 2],
    'gamma' => ['gamma', ['cdf', 'sf', 'logsf', 'quantile'], 2],
    'chi-square' => ['chisquare', ['cdf', 'sf', 'logcdf', 'quantile'], 1],
    't' => ['t', ['cdf', 'sf', 'logsf', 'quantile'], 1],
    'Poisson' => ['poisson', ['cdf', 'sf', 'logcdf', 'quantile'], 1],
    'negative binomial (mu)' => ['nbinom_mu', ['cdf', 'sf', 'logsf', 'quantile'], 2],
]);
