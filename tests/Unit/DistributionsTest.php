<?php

declare(strict_types=1);

// Port of the extension's tests/distributions.php: R references (adjudicated
// against mpmath where R itself is off), scalar/batch consistency, identities,
// boundary conventions and exception categories. Every call goes through pf(),
// so the polyfill is tested whether or not the extension is loaded.

const DIST_FAMILIES = [
    'normal', 't', 'chisquare', 'f', 'exponential', 'logistic', 'lognormal',
    'gamma', 'beta', 'binomial', 'poisson', 'nbinom', 'nbinom_mu',
];

if (! function_exists('distDecode')) {
    function distDecode(mixed $v): float
    {
        return is_string($v) ? ['Inf' => INF, '-Inf' => -INF, 'NaN' => NAN][$v] : (float) $v;
    }

    /** The original suite's tolerance by operation. */
    function distClose(float $actual, float $expected, string $operation): bool
    {
        if ($actual === $expected || (is_nan($actual) && is_nan($expected))) {
            return true;
        }
        if (! is_finite($actual) || ! is_finite($expected)) {
            return false;
        }
        // Probabilities/densities: relative 1e-12 with a floor below the normal range.
        // Log scale: 1e-12 relative to max(1, |value|). Quantiles: 1e-12 relative plus 1e-15.
        $tolerance = match (true) {
            $operation === 'quantile' => 1e-15 + 1e-12 * abs($expected),
            str_starts_with($operation, 'log') => 1e-12 * max(1.0, abs($expected)),
            default => 1e-300 + 1e-12 * abs($expected),
        };

        return abs($actual - $expected) <= $tolerance;
    }

    /** Documented differences from R (adjudicated against mpmath upstream); skipped for R only. */
    function distKnownDifference(array $case, float $x, float $expected, float $actual): bool
    {
        $known = ['f|pdf|Inf' => true, 'f|logpdf|Inf' => true];
        $key = $case['distribution'].'|'.$case['operation'].'|'.(is_infinite($x) ? ($x > 0 ? 'Inf' : '-Inf') : '');
        if (isset($known[$key])) {
            return true;
        }
        if ($case['operation'] !== 'quantile') {
            return false;
        }
        // R's qf underflows to 0 for tiny p with df1 <= 2 (the true quantile is positive).
        if ($case['distribution'] === 'f' && $expected === 0.0 && $x > 0 && $actual > 0 && $actual < 1e-15) {
            return true;
        }
        // R adds 64 machine epsilons to p before discrete searches.
        $discrete = in_array($case['distribution'], ['binomial', 'poisson', 'nbinom', 'nbinom_mu'], true);
        if ($discrete && $x === 1 - 2 ** -52 && $actual > $expected) {
            return true;
        }
        // A point mass at 0 (prob 0): R returns size at p = 1; the definition gives 0.
        return $case['distribution'] === 'binomial' && $case['params']['prob'] == 0 && $x === 1.0;
    }

    /** Running totals across the fixture tests: [compared, adjudicated, skipped]. */
    function distTally(int $compared = 0, int $adjudicated = 0, int $skipped = 0): array
    {
        static $totals = [0, 0, 0];
        $totals[0] += $compared;
        $totals[1] += $adjudicated;
        $totals[2] += $skipped;

        return $totals;
    }

    function distParams(array $case): array
    {
        $params = $case['params'];
        if ($case['distribution'] === 'binomial') {
            $params['size'] = (int) $params['size'];
        }

        return $params;
    }

    function distCode(callable $f): ?int
    {
        try {
            $f();
        } catch (\Throwable $e) {
            return $e->getCode();
        }

        return null;
    }
}

dataset('distribution families', DIST_FAMILIES);

it('matches the R references (mpmath where R is adjudicated wrong)', function (string $family) {
    $r = loadFixture('distributions-r.json');
    $independent = loadFixture('distributions-mpmath.json');
    $compared = $adjudicated = $skipped = 0;
    foreach ($r['cases'] as $caseIndex => $case) {
        if ($case['distribution'] !== $family) {
            continue;
        }
        $params = distParams($case);
        $scalar = 'Phit\dist_'.$case['distribution'].'_'.$case['operation'];
        $input = array_map(distDecode(...), $case['input']);
        $batch = pf($scalar.'_array', $input, ...$params);
        expect(array_is_list($batch) && count($batch) === count($input))->toBeTrue("$scalar batch shape");
        foreach ($input as $i => $x) {
            $one = pf($scalar, $x, ...$params);
            expect($one === $batch[$i] || (is_nan($one) && is_nan($batch[$i])))
                ->toBeTrue("$scalar scalar/batch differ at $x");
            $expected = distDecode($case['expected'][$i]);
            if (distKnownDifference($case, $x, $expected, $one)) {
                $skipped++;

                continue;
            }
            $reference = $independent['cases'][$caseIndex];
            expect($reference['distribution'] === $case['distribution']
                && $reference['operation'] === $case['operation'] && $reference['params'] == $params
                && $reference['input'] == $case['input'])->toBeTrue('independent reference alignment');
            if ($reference['expected'][$i] !== null
                && ! distClose($expected, distDecode($reference['expected'][$i]), $case['operation'])) {
                $mp = distDecode($reference['expected'][$i]);
                expect(distClose($one, $mp, $case['operation']))
                    ->toBeTrue("$scalar(".json_encode($params).") at $x: $one != independent mpmath reference $mp");
                $adjudicated++;
            } else {
                expect(distClose($one, $expected, $case['operation']))
                    ->toBeTrue("$scalar(".json_encode($params).") at $x: $one != R $expected");
            }
            $compared++;
        }
    }
    distTally($compared, $adjudicated, $skipped);
    expect($compared)->toBeGreaterThan(0);
})->with('distribution families');

it('compares more than 5000 reference values', function () {
    [$compared] = distTally();
    expect($compared)->toBeGreaterThan(5000);
});

it('matches SciPy wherever SciPy agrees with the mpmath reference', function (string $family) {
    $scipy = loadFixture('distributions-scipy.json');
    $independent = loadFixture('distributions-mpmath.json');
    $checked = 0;
    foreach ($scipy['cases'] as $caseIndex => $case) {
        if ($case['distribution'] !== $family) {
            continue;
        }
        $params = distParams($case);
        $scalar = 'Phit\dist_'.$case['distribution'].'_'.$case['operation'];
        $reference = $independent['cases'][$caseIndex];
        foreach ($case['input'] as $i => $raw) {
            $x = distDecode($raw);
            $expected = distDecode($case['expected'][$i]);
            if ($reference['expected'][$i] === null
                || ! distClose($expected, distDecode($reference['expected'][$i]), $case['operation'])) {
                continue;
            }
            $one = pf($scalar, $x, ...$params);
            expect(distClose($one, $expected, $case['operation']))
                ->toBeTrue("$scalar(".json_encode($params).") at $x: $one != SciPy $expected");
            $checked++;
        }
    }
    expect($checked)->toBeGreaterThan(0);
})->with('distribution families');

it('satisfies continuous identities: F(Q(p)) = p, F + S = 1, exp(log f) = f, monotone CDF', function (string $name, array $params) {
    foreach ([1e-12, 0.001, 0.025, 0.3, 0.5, 0.8, 0.975, 0.999] as $p) {
        $q = pf("Phit\\dist_{$name}_quantile", $p, ...$params);
        $back = pf("Phit\\dist_{$name}_cdf", $q, ...$params);
        expect(abs($back - $p) <= 1e-12 * $p + 1e-15)->toBeTrue("$name F(Q($p)) = $back");
        $sum = $back + pf("Phit\\dist_{$name}_sf", $q, ...$params);
        expect(abs($sum - 1.0) <= 4 * PHP_FLOAT_EPSILON)->toBeTrue("$name F + S at $q");
        $f = pf("Phit\\dist_{$name}_pdf", $q, ...$params);
        expect(abs(exp(pf("Phit\\dist_{$name}_logpdf", $q, ...$params)) - $f) <= 1e-13 * $f)->toBeTrue("$name exp(logpdf)");
    }
    $grid = array_map(fn ($p) => pf("Phit\\dist_{$name}_quantile", $p, ...$params), range(0.01, 0.99, 0.01));
    $cdf = pf("Phit\\dist_{$name}_cdf_array", $grid, ...$params);
    for ($i = 1; $i < count($cdf); $i++) {
        expect($cdf[$i] >= $cdf[$i - 1])->toBeTrue("$name CDF monotone");
    }
})->with([
    ['normal', [0.0, 1.0]], ['normal', [3.0, 2.0]], ['t', [3.0]], ['t', [30.0]], ['chisquare', [4.0]],
    ['f', [4.0, 40.0]], ['beta', [2.0, 5.0]], ['gamma', [2.0, 3.0]], ['exponential', [0.5]],
    ['lognormal', [0.0, 1.0]], ['logistic', [0.0, 1.0]],
]);

it('satisfies discrete identities and step-function conventions', function (string $name, array $params) {
    $total = 0.0;
    foreach (range(0, 200) as $k) {
        $total += pf("Phit\\dist_{$name}_pmf", $k, ...$params);
    }
    expect(abs($total - 1.0) < 1e-12)->toBeTrue("$name pmf sums to one");
    foreach ([0.01, 0.2, 0.5, 0.9, 0.999] as $p) {
        $q = pf("Phit\\dist_{$name}_quantile", $p, ...$params);
        expect(pf("Phit\\dist_{$name}_cdf", $q, ...$params) >= $p)->toBeTrue("$name F(Q(p)) >= p");
        expect($q == 0 || pf("Phit\\dist_{$name}_cdf", $q - 1, ...$params) < $p)->toBeTrue("$name Q(p) is the smallest such k");
    }
    // Step-function conventions: cdf floors x; P(X > k) excludes k; pmf is 0 off the integers.
    expect(pf("Phit\\dist_{$name}_cdf", 2.7, ...$params))->toBe(pf("Phit\\dist_{$name}_cdf", 2, ...$params));
    $sf = pf("Phit\\dist_{$name}_sf", 2, ...$params);
    expect(abs($sf - (1 - pf("Phit\\dist_{$name}_cdf", 2, ...$params))) < 1e-15)->toBeTrue("$name sf is P(X > k)");
    expect(abs(pf("Phit\\dist_{$name}_sf", 1, ...$params) - ($sf + pf("Phit\\dist_{$name}_pmf", 2, ...$params))) < 1e-15)
        ->toBeTrue("$name P(X >= k) = sf(k - 1)");
    expect(pf("Phit\\dist_{$name}_pmf", 2.5, ...$params))->toBe(0.0)
        ->and(pf("Phit\\dist_{$name}_logpmf", 2.5, ...$params))->toBe(-INF)
        ->and(pf("Phit\\dist_{$name}_pmf", -1, ...$params))->toBe(0.0)
        ->and(pf("Phit\\dist_{$name}_sf", -1, ...$params))->toBe(1.0);
})->with([
    ['binomial', [10, 0.3]], ['poisson', [4.5]], ['nbinom', [3.0, 0.4]], ['nbinom_mu', [2.0, 4.5]],
]);

it('keeps extreme tails finite on the log scale', function () {
    $sf = pf('Phit\dist_normal_sf', 40.0);
    expect($sf === 0.0 || $sf < 1e-300)->toBeTrue('normal tail underflows')
        ->and(abs(pf('Phit\dist_normal_logsf', 40.0) - -804.6084420137538) < 1e-11)->toBeTrue('normal log tail')
        ->and(abs(pf('Phit\dist_normal_logcdf', -1e10) - -5e19) < 1e5)->toBeTrue('normal log tail at -1e10')
        ->and(is_finite(pf('Phit\dist_t_logsf', 1e5, 100.0)))->toBeTrue()
        ->and(is_finite(pf('Phit\dist_chisquare_logsf', 1e4, 0.5)))->toBeTrue()
        ->and(is_finite(pf('Phit\dist_poisson_logsf', 1e7, 3.2)))->toBeTrue()
        ->and(is_finite(pf('Phit\dist_binomial_logcdf', 0, 1000, 0.999)))->toBeTrue()
        ->and(pf('Phit\dist_t_quantile', 0.5, 3))->toBe(0.0)
        ->and(pf('Phit\dist_normal_quantile', 0.5))->toBe(0.0);
});

it('handles boundaries: infinite points and endpoint quantiles', function () {
    expect(pf('Phit\dist_normal_cdf', INF))->toBe(1.0)
        ->and(pf('Phit\dist_normal_cdf', -INF))->toBe(0.0)
        ->and(pf('Phit\dist_normal_pdf', INF))->toBe(0.0)
        ->and(pf('Phit\dist_normal_quantile', 0.0))->toBe(-INF)
        ->and(pf('Phit\dist_normal_quantile', 1.0))->toBe(INF)
        ->and(pf('Phit\dist_gamma_quantile', 0.0, 2.0))->toBe(0.0)
        ->and(pf('Phit\dist_beta_quantile', 1.0, 2.0, 3.0))->toBe(1.0)
        ->and(pf('Phit\dist_binomial_quantile', 1.0, 7, 0.0))->toBe(0.0)
        ->and(pf('Phit\dist_poisson_quantile', 1.0, 3.0))->toBe(INF)
        ->and(pf('Phit\dist_binomial_quantile', 1.0, 7, 0.4))->toBe(7.0)
        ->and(pf('Phit\dist_gamma_pdf', 0.0, 0.5))->toBe(INF)
        ->and(pf('Phit\dist_gamma_pdf', 0.0, 1.0, 2.0))->toBe(0.5)
        ->and(pf('Phit\dist_gamma_pdf', -1.0, 2.0))->toBe(0.0)
        ->and(pf('Phit\dist_poisson_pmf', 0, 0.0))->toBe(1.0)
        ->and(pf('Phit\dist_binomial_pmf', 3, 3, 1.0))->toBe(1.0)
        ->and(pf('Phit\dist_nbinom_mu_pmf', 0, 1e308, 1e308))->toBe(0.0)
        ->and(pf('Phit\dist_nbinom_mu_logpmf', 0, 1e308, 1e308))->toBe(-1e308 * log(2.0))
        ->and(pf('Phit\dist_nbinom_mu_cdf', 0, 1e308, 1e308))->toBe(0.0)
        ->and(pf('Phit\dist_nbinom_mu_sf', 0, 1e308, 1e308))->toBe(1.0)
        ->and(pf('Phit\dist_nbinom_mu_pmf', 0, 3.0, 0.0))->toBe(1.0);
});

it('supports defaults, integer arguments and named arguments', function () {
    expect(pf('Phit\dist_normal_cdf', 1.96))->toBe(pf('Phit\dist_normal_cdf', x: 1.96, mean: 0, sd: 1))
        ->and(pf('Phit\dist_t_quantile', 0.975, 30))->toBe(pf('Phit\dist_t_quantile', p: 0.975, df: 30.0))
        ->and(pf('Phit\dist_gamma_cdf', 2.0, 2.0))->toBe(pf('Phit\dist_gamma_cdf', 2.0, shape: 2.0, scale: 1.0))
        ->and(pf('Phit\dist_exponential_cdf', 1.0))->toBe(pf('Phit\dist_exponential_cdf', 1.0, rate: 1.0));
});

it('validates batch input: empty arrays, order, keys and element types', function () {
    expect(pf('Phit\dist_normal_cdf_array', []))->toBe([])
        ->and(pf('Phit\dist_normal_cdf_array', ['b' => 1.0, 'a' => -1.0]))
        ->toBe([pf('Phit\dist_normal_cdf', 1.0), pf('Phit\dist_normal_cdf', -1.0)])
        ->and(pf('Phit\dist_poisson_logpmf_array', counts: [3, 4], lambda: 4.5))
        ->toBe([pf('Phit\dist_poisson_logpmf', 3, 4.5), pf('Phit\dist_poisson_logpmf', 4, 4.5)])
        ->and(distCode(fn () => pf('Phit\dist_normal_cdf_array', [], sd: -1.0)))->toBe(1001);
    foreach ([[1.0, 'x'], [1.0, null], [true], [[1.0]], ['1.5']] as $bad) {
        expect(distCode(fn () => pf('Phit\dist_normal_cdf_array', $bad)))->toBe(1001);
    }
    expect(distCode(fn () => pf('Phit\dist_normal_cdf_array', [1.0, NAN])))->toBe(1003)
        ->and(distCode(fn () => pf('Phit\dist_normal_quantile_array', [0.5, 1.5])))->toBe(1001)
        ->and(distCode(fn () => pf('Phit\dist_normal_quantile_array', [0.5, NAN])))->toBe(1003)
        // The first bad element decides the error.
        ->and(distCode(fn () => pf('Phit\dist_normal_cdf_array', [NAN, 'x'])))->toBe(1003)
        ->and(distCode(fn () => pf('Phit\dist_normal_cdf_array', ['x', NAN])))->toBe(1001);
});

it('raises 1001 for invalid parameters and probabilities', function (Closure $f) {
    expect(distCode($f))->toBe(1001);
})->with([
    'normal sd 0' => fn () => pf('Phit\dist_normal_cdf', 0.0, 0.0, 0.0),
    'normal sd < 0' => fn () => pf('Phit\dist_normal_cdf', 0.0, 0.0, -1.0),
    't df 0' => fn () => pf('Phit\dist_t_cdf', 0.0, 0.0),
    'chisquare df < 0' => fn () => pf('Phit\dist_chisquare_cdf', 1.0, -2.0),
    'f df1 0' => fn () => pf('Phit\dist_f_cdf', 1.0, 0.0, 3.0),
    'beta alpha 0' => fn () => pf('Phit\dist_beta_cdf', 0.5, 0.0, 1.0),
    'gamma scale 0' => fn () => pf('Phit\dist_gamma_cdf', 1.0, 2.0, 0.0),
    'exponential rate < 0' => fn () => pf('Phit\dist_exponential_cdf', 1.0, -1.0),
    'lognormal sdlog 0' => fn () => pf('Phit\dist_lognormal_cdf', 1.0, 0.0, 0.0),
    'logistic scale < 0' => fn () => pf('Phit\dist_logistic_cdf', 1.0, 0.0, -2.0),
    'binomial size < 0' => fn () => pf('Phit\dist_binomial_pmf', 1, -1, 0.5),
    'binomial prob > 1' => fn () => pf('Phit\dist_binomial_pmf', 1, 3, 1.5),
    'poisson lambda < 0' => fn () => pf('Phit\dist_poisson_pmf', 1, -0.5),
    'nbinom size 0' => fn () => pf('Phit\dist_nbinom_pmf', 1, 0.0, 0.5),
    'nbinom prob 0' => fn () => pf('Phit\dist_nbinom_pmf', 1, 2.0, 0.0),
    'nbinom_mu mu < 0' => fn () => pf('Phit\dist_nbinom_mu_pmf', 1, 2.0, -1.0),
    'quantile p < 0' => fn () => pf('Phit\dist_normal_quantile', -0.1),
    'quantile p > 1' => fn () => pf('Phit\dist_poisson_quantile', 1.1, 2.0),
    'chisquare df/2 underflows' => fn () => pf('Phit\dist_chisquare_cdf', 1.0, 5e-324),
]);

it('raises 1003 for nonfinite parameters and NaN inputs', function (Closure $f) {
    expect(distCode($f))->toBe(1003);
})->with([
    'NaN point' => fn () => pf('Phit\dist_normal_cdf', NAN),
    'infinite mean' => fn () => pf('Phit\dist_normal_cdf', 0.0, INF),
    'NaN probability' => fn () => pf('Phit\dist_t_quantile', NAN, 3.0),
    'NaN count' => fn () => pf('Phit\dist_poisson_pmf', NAN, 2.0),
    'NaN shape' => fn () => pf('Phit\dist_gamma_pdf', 1.0, NAN),
    'infinite df' => fn () => pf('Phit\dist_t_cdf', 1.0, INF),
    'NaN binomial prob beats size limit' => fn () => pf('Phit\dist_binomial_pmf', 1, 2 ** 54, NAN),
]);

it('raises 1007 for binomial sizes beyond 2^53', function () {
    expect(distCode(fn () => pf('Phit\dist_binomial_pmf', 1, 2 ** 54, 0.5)))->toBe(1007)
        ->and(distCode(fn () => pf('Phit\dist_binomial_cdf_array', [], 2 ** 53 + 1, 0.5)))->toBe(1007)
        ->and(pf('Phit\dist_binomial_cdf', 0, 2 ** 53, 0.0))->toBe(1.0);
});

it('throws plain exceptions carrying the error code', function () {
    try {
        pf('Phit\dist_normal_cdf', 0.0, 0.0, -1.0);
        $this->fail('expected an exception');
    } catch (\Exception $e) {
        expect(get_class($e))->toBe(\Exception::class)->and($e->getCode())->toBe(1001);
    }
});

it('returns lists with the same order for every batch function', function (string $family) {
    $values = ['z' => 0.25, 'y' => 3, 7 => 0.75];
    $params = match ($family) {
        'normal', 'lognormal', 'logistic', 'exponential' => [],
        't', 'chisquare' => [3.0],
        'f', 'beta' => [3.0, 4.0],
        'gamma' => [2.0],
        'binomial' => [10, 0.4],
        'poisson' => [2.5],
        'nbinom' => [2.0, 0.4],
        'nbinom_mu' => [2.0, 3.0],
    };
    $discrete = in_array($family, ['binomial', 'poisson', 'nbinom', 'nbinom_mu'], true);
    foreach ([$discrete ? 'pmf' : 'pdf', $discrete ? 'logpmf' : 'logpdf', 'cdf', 'logcdf', 'sf', 'logsf', 'quantile'] as $op) {
        $input = $op === 'quantile' ? ['a' => 0.1, 'b' => 0.9, 3 => 0.5] : $values;
        $batch = pf("Phit\\dist_{$family}_{$op}_array", $input, ...$params);
        $scalars = array_map(fn ($v) => pf("Phit\\dist_{$family}_{$op}", $v, ...$params), array_values($input));
        expect($batch)->toBe($scalars);
    }
})->with('distribution families');

it('finds the same quantile bracket as the literal upstream step-by-2 walk', function () {
    $bracket = (new ReflectionClass(\Phit\Polyfill\Distribution\Quantile::class))->getMethod('bracket');
    $expand = (new ReflectionClass(\Phit\Polyfill\Distribution\Quantile::class))->getMethod('expand');
    // Literal port of the extension's loop (quantile.rs), for comparison.
    $walk = function (float $x, float $direction, Closure $stop) use ($expand): ?array {
        $before = $x;
        $next = $x;
        while (true) {
            $next = $expand->invoke(null, $next, $direction);
            if ($next == $direction * INF) {
                return null;
            }
            if ($stop($next)) {
                return [$next, $before];
            }
            $before = $next;
        }
    };
    $cases = [
        // [start, direction, root]: stop once the point passes the root.
        [1000.5, -1.0, 37.25], [1000.5, -1.0, 998.7], [1000.5, -1.0, 999.0], [1000.5, -1.0, -5.0],
        [1000.5, -1.0, -1e40], [3.0, -1.0, 0.3], [2.0, -1.0, 1.0], [123456789.123, -1.0, 54321.5],
        [-1000.5, 1.0, -37.25], [-1000.5, 1.0, -998.5], [-1000.5, 1.0, 0.7], [-1000.5, 1.0, 1e40],
        [-0.25, 1.0, 1e300], [-7.0, 1.0, -5.0], [-987654.321, 1.0, -2.5], [-4.0, 1.0, -2.0],
    ];
    foreach ($cases as [$start, $direction, $root]) {
        $stop = $direction < 0.0
            ? static fn (float $v): bool => $v <= $root
            : static fn (float $v): bool => $v >= $root;
        expect($bracket->invoke(null, $start, $direction, $stop))->toBe($walk($start, $direction, $stop));
    }
});

// phit bd98bf7: a NaN tail stops the quantile bracket search and fails with
// 1006 instead of walking to ±INF. The stop test checks NaN explicitly, because
// with opcache enabled a negated comparison did not treat NaN as stopping.
it('fails a quantile whose tail evaluation fails during bracketing', function () {
    expect(thrownCode(fn () => pf('Phit\dist_t_quantile', 1.0 - 0.999999999999, 1e299)))->toBe(1006)
        ->and(thrownCode(fn () => pf('Phit\dist_t_quantile', 0.999999999999, 1e299)))->toBe(1006)
        ->and(thrownCode(fn () => pf('Phit\dist_f_quantile', 0.999999999999, 1e-10, 1e299)))->toBe(1006);
});

it('fails fast on the cdflib regimes that hang or panic upstream', function (string $fn, array $arguments) {
    expect(thrownCode(fn () => pf($fn, ...$arguments)))->toBe(1006);
})->with([
    'F cdf, df2 = MAX' => ['Phit\dist_f_cdf', [2.0, 3.0, PHP_FLOAT_MAX]],
    'NB(mu) cdf, size = MAX' => ['Phit\dist_nbinom_mu_cdf', [2.0, PHP_FLOAT_MAX, 3.0]],
    'beta quantile, beta = MAX' => ['Phit\dist_beta_quantile', [1e-10, 1.5, PHP_FLOAT_MAX]],
    'gamma cdf, shape = x = 1e30' => ['Phit\dist_gamma_cdf', [1e30, 1e30]],
    'chi-square cdf, df = x = 2e30' => ['Phit\dist_chisquare_cdf', [2e30, 2e30]],
    'Poisson cdf, lambda = k = 1e30' => ['Phit\dist_poisson_cdf', [1e30, 1e30]],
    'beta cdf, alpha = beta = 1e308' => ['Phit\dist_beta_cdf', [0.5, 1e308, 1e308]],
]);
