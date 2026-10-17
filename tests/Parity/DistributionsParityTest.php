<?php

declare(strict_types=1);

// Compares every polyfilled distribution function with the loaded extension on
// the fixture inputs and on seeded random inputs (tails, ±INF, p in {0, 1, tiny},
// extreme parameters), including exception categories. Skipped without it.
// PARITY_REPORT=1 prints the largest relative difference per family.

const DIST_PARITY_REL = 1e-13;

const DIST_PARITY_FAMILIES = [
    'normal', 't', 'chisquare', 'f', 'exponential', 'logistic', 'lognormal',
    'gamma', 'beta', 'binomial', 'poisson', 'nbinom', 'nbinom_mu',
];

if (! function_exists('distParityOutcome')) {
    /** [value or null, exception signature or null] */
    function distParityOutcome(callable $call): array
    {
        try {
            return [$call(), null];
        } catch (\Throwable $e) {
            return [null, get_class($e).'#'.$e->getCode()];
        }
    }

    /** Largest relative difference per family; pass a value to record it. */
    function distParityMax(?string $family = null, ?float $seen = null): array
    {
        static $max = [];
        if ($family !== null) {
            $max[$family] = max($max[$family] ?? 0.0, $seen ?? 0.0);
        }

        return $max;
    }

    function distParityDiff(float $a, float $b): float
    {
        if ($a === $b || (is_nan($a) && is_nan($b))) {
            return 0.0;
        }
        if (! is_finite($a) || ! is_finite($b)) {
            return INF;
        }
        $scale = max(abs($a), abs($b));

        return $scale < 1e-300 ? abs($a - $b) : abs($a - $b) / $scale;
    }

    function distParityCompare(string $family, string $fn, array $args): void
    {
        [$pf, $pfError] = distParityOutcome(fn () => pf($fn, ...$args));
        [$nat, $natError] = distParityOutcome(fn () => native($fn, ...$args));
        $label = fn () => $fn.' '.json_encode($args, JSON_PARTIAL_OUTPUT_ON_ERROR | JSON_PRESERVE_ZERO_FRACTION)
            .' pf='.var_export($pfError ?? $pf, true).' native='.var_export($natError ?? $nat, true);
        if ($pfError !== $natError) {
            throw new \RuntimeException('exception mismatch: '.$label());
        }
        expect($pfError)->toBe($natError);
        if ($pfError !== null) {
            return;
        }
        $pfValues = is_array($pf) ? $pf : [$pf];
        $natValues = is_array($nat) ? $nat : [$nat];
        if (array_keys($pfValues) !== array_keys($natValues)) {
            throw new \RuntimeException('shape mismatch: '.$label());
        }
        $worst = 0.0;
        foreach ($pfValues as $i => $value) {
            $worst = max($worst, distParityDiff($value, $natValues[$i]));
        }
        distParityMax($family, $worst);
        if ($worst > DIST_PARITY_REL) {
            throw new \RuntimeException(sprintf('value mismatch (%.3e): %s', $worst, $label()));
        }
    }

    function distParityUniform(float $lo, float $hi): float
    {
        return $lo + ($hi - $lo) * mt_rand() / mt_getrandmax();
    }

    function distParityLogUniform(float $lo, float $hi): float
    {
        return exp(distParityUniform(log($lo), log($hi)));
    }

    function distParityPick(array $values): mixed
    {
        return $values[mt_rand(0, count($values) - 1)];
    }

    /** A random parameter list for $family, extremes included. */
    function distParityParams(string $family): array
    {
        $scale = fn () => distParityPick([1.0, distParityLogUniform(1e-3, 1e3), distParityLogUniform(1e-300, 1e300)]);
        $df = fn (float $lo, float $hi) => distParityPick([1.0, 2.0, 0.5, distParityLogUniform($lo, $hi), distParityLogUniform($lo, $hi)]);

        return match ($family) {
            'normal' => [distParityPick([0.0, distParityUniform(-1e3, 1e3), 1e300]), $scale()],
            // df below ~0.05 makes the extension's quantile bracket walk take seconds per call;
            // that regime has its own test below.
            't' => [$df(5e-2, 1e12)],
            'chisquare' => [$df(1e-3, 1e8)],
            'f' => [$df(1e-2, 1e8), $df(1e-2, 1e8)],
            'exponential' => [distParityPick([1.0, distParityLogUniform(1e-100, 1e100)])],
            'logistic' => [distParityPick([0.0, distParityUniform(-1e3, 1e3)]), $scale()],
            'lognormal' => [distParityUniform(-50.0, 50.0), distParityLogUniform(1e-3, 1e2)],
            'gamma' => [$df(1e-3, 1e8), distParityPick([1.0, distParityLogUniform(1e-10, 1e10)])],
            'beta' => [$df(1e-3, 1e6), $df(1e-3, 1e6)],
            'binomial' => [
                distParityPick([0, 1, mt_rand(0, 100), mt_rand(0, 1000000), 2 ** 53, 1000000000000]),
                distParityPick([0.0, 1.0, distParityUniform(0.0, 1.0), 1e-10, 1.0 - 1e-10, distParityUniform(0.0, 1.0)]),
            ],
            'poisson' => [distParityPick([0.0, distParityLogUniform(1e-10, 1e9), distParityLogUniform(0.1, 100.0)])],
            'nbinom' => [$df(1e-3, 1e6), distParityPick([1.0, distParityLogUniform(1e-10, 1.0), distParityUniform(0.0, 1.0)])],
            'nbinom_mu' => [$df(1e-3, 1e6), distParityPick([0.0, distParityLogUniform(1e-6, 1e9), distParityUniform(0.0, 20.0)])],
        };
    }

    /** Evaluation points: random points around the bulk plus tails, zeros and infinities. */
    function distParityPoints(string $family, array $params): array
    {
        $points = [0.0, -0.0, INF, -INF, 1.0, -1.0, 1e-300, -1e-300, 5e-324, 1e10, -1e10, 1e300, -1e300];
        $discrete = in_array($family, ['binomial', 'poisson', 'nbinom', 'nbinom_mu'], true);
        if ($discrete) {
            [$mean, $sd] = match ($family) {
                'binomial' => [$params[0] * $params[1], sqrt($params[0] * $params[1] * (1 - $params[1]))],
                'poisson' => [$params[0], sqrt($params[0])],
                'nbinom' => [$params[0] * (1 - $params[1]) / $params[1], sqrt($params[0] * (1 - $params[1])) / $params[1]],
                'nbinom_mu' => [$params[1], sqrt($params[1] + $params[1] ** 2 / $params[0])],
            };
            for ($i = 0; $i < 8; $i++) {
                $points[] = max(0.0, floor($mean + distParityUniform(-8.0, 8.0) * $sd));
            }
            array_push($points, 2.5, 3, 7, floor($mean), floor($mean) + 0.5, 2 ** 53, 1e15);
        } else {
            for ($i = 0; $i < 6; $i++) {
                $points[] = distParityUniform(-5.0, 5.0);
                $points[] = distParityLogUniform(1e-300, 1e300);
                $points[] = -distParityLogUniform(1e-10, 1e10);
            }
            if ($family === 'normal' || $family === 'logistic') {
                for ($i = 0; $i < 6; $i++) {
                    $points[] = $params[0] + $params[1] * distParityUniform(-60.0, 60.0);
                }
            }
            if ($family === 'beta') {
                array_push($points, 0.5, 1e-10, 1 - 1e-10, distParityUniform(0.0, 1.0), distParityUniform(0.0, 1.0), distParityUniform(0.0, 1.0));
            }
        }

        return $points;
    }

    function distParityProbabilities(): array
    {
        $p = [0.0, 1.0, 0.5, 5e-324, 1e-300, 1e-100, 1e-10, 1e-3, 0.975, 1 - 1e-10, 1 - 2 ** -53, 1 - 2 ** -52];
        for ($i = 0; $i < 6; $i++) {
            $p[] = distParityUniform(0.0, 1.0);
        }
        $p[] = distParityLogUniform(1e-300, 1e-3);

        return $p;
    }

    function distParityOperations(string $family): array
    {
        $discrete = in_array($family, ['binomial', 'poisson', 'nbinom', 'nbinom_mu'], true);

        return [$discrete ? 'pmf' : 'pdf', $discrete ? 'logpmf' : 'logpdf', 'cdf', 'logcdf', 'sf', 'logsf', 'quantile'];
    }

    function distParityDecode(mixed $v): float
    {
        return is_string($v) ? ['Inf' => INF, '-Inf' => -INF, 'NaN' => NAN][$v] : (float) $v;
    }
}

beforeEach(function () {
    if (! extension_loaded('phit')) {
        $this->markTestSkipped('phit extension not loaded');
    }
});

dataset('parity families', DIST_PARITY_FAMILIES);

it('polyfills every distribution function', function (string $family) {
    foreach (distParityOperations($family) as $op) {
        foreach (['', '_array'] as $suffix) {
            expect(\Phit\Polyfill\Polyfill::supports("Phit\\dist_{$family}_{$op}{$suffix}"))->toBeTrue();
        }
    }
})->with('parity families');

it('agrees with the extension on the fixture inputs', function (string $family) {
    foreach (['distributions-r.json', 'distributions-scipy.json'] as $name) {
        foreach (loadFixture($name)['cases'] as $case) {
            if ($case['distribution'] !== $family) {
                continue;
            }
            $params = $case['params'];
            if ($family === 'binomial') {
                $params['size'] = (int) $params['size'];
            }
            $fn = "Phit\\dist_{$family}_{$case['operation']}";
            $input = array_map(distParityDecode(...), $case['input']);
            distParityCompare($family, $fn.'_array', [$input, ...$params]);
            foreach ($input as $x) {
                distParityCompare($family, $fn, [$x, ...$params]);
            }
        }
    }
})->with('parity families');

it('agrees with the extension on seeded random inputs', function (string $family) {
    mt_srand(crc32($family) + 20261015);
    for ($set = 0; $set < 25; $set++) {
        $params = distParityParams($family);
        $points = distParityPoints($family, $params);
        $probabilities = distParityProbabilities();
        foreach (distParityOperations($family) as $op) {
            $input = $op === 'quantile' ? $probabilities : $points;
            $fn = "Phit\\dist_{$family}_{$op}";
            distParityCompare($family, $fn.'_array', [$input, ...$params]);
            foreach ($input as $x) {
                distParityCompare($family, $fn, [$x, ...$params]);
            }
        }
    }
})->with('parity families');

it('raises the same exceptions on invalid input', function (string $family) {
    $valid = match ($family) {
        'normal', 'logistic', 'lognormal' => [0.5, 2.0],
        't', 'chisquare', 'exponential', 'poisson' => [3.0],
        'f', 'beta', 'nbinom_mu' => [2.0, 3.0],
        'gamma' => [2.0, 1.5],
        'binomial' => [10, 0.3],
        'nbinom' => [3.0, 0.4],
    };
    // Not PHP_FLOAT_MAX: the extension never returns from e.g. Phit\dist_f_cdf(2, DBL_MAX, 3)
    // (a NaN-driven infinite loop in cdflib's beta_frac); the polyfill reports 1006 there.
    $bad = [NAN, INF, -INF, 0.0, -0.0, -1.0, 5e-324, 1.5, 2.0, 1e300];
    $cases = [];
    foreach (array_keys($valid) as $i) {
        foreach ($bad as $value) {
            $params = $valid;
            if ($family === 'binomial' && $i === 0) {
                foreach ([-1, -PHP_INT_MAX, 0, 2 ** 53, 2 ** 53 + 1, PHP_INT_MAX] as $size) {
                    $cases[] = [$size, $valid[1]];
                }

                continue 2;
            }
            $params[$i] = $value;
            $cases[] = $params;
        }
    }
    foreach ($cases as $params) {
        foreach (distParityOperations($family) as $op) {
            $fn = "Phit\\dist_{$family}_{$op}";
            foreach ([0.5, 2.0, NAN, -INF] as $x) {
                distParityCompare($family, $fn, [$x, ...$params]);
            }
            distParityCompare($family, $fn.'_array', [[], ...$params]);
        }
    }
    foreach (distParityOperations($family) as $op) {
        $fn = "Phit\\dist_{$family}_{$op}";
        foreach ([-0.1, 1.1, NAN, INF, -INF] as $x) {
            distParityCompare($family, $fn, [$x, ...$valid]);
        }
        foreach ([[1.0, 'x'], [1.0, null], [true], [[1.0]], ['0.5'], [NAN, 'x'], ['x', NAN], [0.5, 1.5], ['k' => 0.25, 'j' => 0.75]] as $input) {
            distParityCompare($family, $fn.'_array', [$input, ...$valid]);
        }
    }
})->with('parity families');

it('agrees on Student t quantiles with tiny df (long upstream bracket walks)', function () {
    foreach ([[0.3, 0.01], [0.1, 0.01], [0.45, 0.02], [0.2, 0.05], [0.01, 0.05], [0.7, 0.01],
        [0.3937733622238847, 0.0038418888682897044]] as [$p, $df]) {
        distParityCompare('t', 'Phit\dist_t_quantile', [$p, $df]);
    }
    distParityCompare('t', 'Phit\dist_t_quantile_array', [[1e-300, 0.01, 0.25, 0.5, 0.75, 0.99], 0.02]);
});

it('agrees in the indeterminate incomplete-gamma regime (a Rust panic upstream)', function () {
    distParityCompare('gamma', 'Phit\dist_gamma_cdf', [1e30, 1e30]);
    distParityCompare('chisquare', 'Phit\dist_chisquare_sf', [2e30, 2e30]);
    distParityCompare('poisson', 'Phit\dist_poisson_cdf', [1e30, 1e30]);
});

it('reports the largest relative difference per family', function () {
    $max = distParityMax();
    if (getenv('PARITY_REPORT')) {
        foreach ($max as $family => $value) {
            file_put_contents('php://stderr', sprintf("%-12s max relative difference: %.3e\n", $family, $value));
        }
    }
    expect(count($max))->toBe(count(DIST_PARITY_FAMILIES));
    foreach ($max as $value) {
        expect($value)->toBeLessThanOrEqual(DIST_PARITY_REL);
    }
});
