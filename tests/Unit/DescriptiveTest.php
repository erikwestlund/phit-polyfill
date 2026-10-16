<?php

declare(strict_types=1);

// Port of the extension's tests/descriptive.php. Relative accuracy 1e-12; exact zeros stay exact.
function dclose(mixed $actual, float $expected, string $name = ''): void
{
    expect(is_finite($actual))->toBeTrue($name);
    expect($actual)->toBeClose($expected, 0.0, 1e-12);
}

function dcode(callable $f, int $expected): void
{
    expect(thrownCode($f))->toBe($expected);
}

$fixture = loadFixture('descriptive-r.json');

it('matches R for single-sample statistics', function () use ($fixture) {
    foreach ($fixture['descriptive'] as $c) {
        foreach (['mean' => 'mean', 'variance' => 'var', 'sd' => 'sd', 'median' => 'median'] as $key => $fn) {
            dclose(pf('stats_'.$fn, $c['values']), $c[$key], "{$c['name']} $key");
        }
        dclose(pf('stats_var', $c['values'], type: 'population'), $c['population_variance']);
        dclose(pf('stats_sd', $c['values'], type: 'population'), sqrt($c['population_variance']));
        $batch = pf('stats_quantile_array', $c['values'], $c['quantile_probabilities']);
        expect(array_is_list($batch))->toBeTrue()->and($batch)->toHaveCount(count($c['quantiles']));
        foreach ($batch as $j => $v) {
            dclose($v, $c['quantiles'][$j], "{$c['name']} quantile $j");
            expect($v)->toBe(pf('stats_quantile', $c['values'], $c['quantile_probabilities'][$j]));
        }
    }
});

it('matches R for covariance and correlation', function () use ($fixture) {
    foreach ($fixture['pairs'] as $c) {
        dclose(pf('stats_cov', $c['x'], $c['y']), $c['covariance']);
        dclose(pf('stats_cov', $c['x'], $c['y'], type: 'population'), $c['population_covariance']);
        dclose(pf('stats_cor', $c['x'], $c['y']), $c['correlation']);
    }
});

it('matches R for weighted statistics', function () use ($fixture) {
    foreach ($fixture['weighted'] as $c) {
        dclose(pf('stats_weighted_mean', $c['values'], $c['weights']), $c['mean']);
        foreach (['frequency', 'reliability', 'population'] as $type) {
            dclose(pf('stats_weighted_var', $c['values'], $c['weights'], type: $type), $c[$type]);
        }
        if (isset($c['expanded_variance'])) {
            dclose(pf('stats_weighted_var', $c['values'], $c['weights']), $c['expanded_variance']);
        }
        if (isset($c['cov_wt_unbiased'])) {
            dclose(pf('stats_weighted_var', $c['values'], $c['weights'], 'reliability'), $c['cov_wt_unbiased']);
        }
    }
});

it('rejects missing values by default and drops them on request without warnings', function () {
    foreach ([null, NAN] as $missing) {
        foreach (['stats_mean', 'stats_var', 'stats_sd', 'stats_median'] as $fn) {
            dcode(fn () => pf($fn, [1, $missing, 3]), 1001);
            expect(warnings(fn () => pf($fn, [1, $missing, 3], missing: 'drop')))->toBe([]);
            dclose(pf($fn, [1, $missing, 3], missing: 'drop'), pf($fn, [1, 3]));
        }
        dcode(fn () => pf('stats_quantile', [1, $missing, 3], .5), 1001);
        dclose(pf('stats_quantile', [1, $missing, 3], .5, missing: 'drop'), 2);
        expect(pf('stats_quantile_array', [1, $missing, 3], [0, .5, 1], 'drop'))->toBe([1.0, 2.0, 3.0]);
        foreach (['stats_cov', 'stats_cor'] as $fn) {
            dcode(fn () => pf($fn, [1, $missing, 3, 4], [2, 4, $missing, 8]), 1001);
            dclose(pf($fn, [1, $missing, 3, 4], [2, 4, $missing, 8], missing: 'drop'), pf($fn, [1, 4], [2, 8]));
        }
        foreach (['stats_weighted_mean', 'stats_weighted_var'] as $fn) {
            dcode(fn () => pf($fn, [1, $missing, 3, 4], [1, 2, $missing, 3]), 1001);
            dclose(pf($fn, [1, $missing, 3, 4], [1, 2, $missing, 3], missing: 'drop'), pf($fn, [1, 4], [1, 3]));
        }
    }
});

it('rejects non-numeric elements', function () {
    foreach ([true, '1', [], new stdClass()] as $bad) {
        dcode(fn () => pf('stats_mean', [1, $bad]), 1001);
        dcode(fn () => pf('stats_cov', [1, null], [2, $bad], missing: 'drop'), 1001);
        dcode(fn () => pf('stats_weighted_mean', [1, null], [2, $bad], missing: 'drop'), 1001);
        dcode(fn () => pf('stats_quantile_array', [1, 2], [$bad]), 1001);
    }
});

it('rejects infinite elements', function () {
    foreach ([INF, -INF] as $bad) {
        dcode(fn () => pf('stats_mean', [1, $bad]), 1003);
        dcode(fn () => pf('stats_cov', [1, null], [2, $bad], missing: 'drop'), 1003);
        dcode(fn () => pf('stats_weighted_mean', [1, null], [2, $bad], missing: 'drop'), 1003);
    }
});

it('validates sizes, options and probabilities', function () {
    dcode(fn () => pf('stats_mean', []), 1001);
    dcode(fn () => pf('stats_mean', [null, NAN], 'drop'), 1001);
    dcode(fn () => pf('stats_var', [1]), 1001);
    dcode(fn () => pf('stats_sd', [1]), 1001);
    dcode(fn () => pf('stats_cov', [1], [2], 'population'), 1001);
    dcode(fn () => pf('stats_cov', [1, 2], [3]), 1002);
    dcode(fn () => pf('stats_weighted_mean', [1, 2], [1]), 1002);
    dcode(fn () => pf('stats_mean', [1], missing: 'ignore'), 1001);
    dcode(fn () => pf('stats_var', [1, 2], type: 'unbiased'), 1001);
    dcode(fn () => pf('stats_cov', [1, 2], [3, 4], type: 'unbiased'), 1001);
    dcode(fn () => pf('stats_weighted_var', [1, 2], [1, 2], type: 'sample'), 1001);
    foreach ([-.1, 1.1, INF, -INF] as $bad) {
        dcode(fn () => pf('stats_quantile', [1, 2], $bad), 1001);
        dcode(fn () => pf('stats_quantile_array', [1, 2], [$bad]), 1001);
    }
    dcode(fn () => pf('stats_quantile', [1, 2], NAN), 1003);
    dcode(fn () => pf('stats_quantile_array', [1, 2], [NAN]), 1003);
});

it('validates weights and overflow', function () {
    dcode(fn () => pf('stats_weighted_mean', [1, 2], [0, 0]), 1001);
    dcode(fn () => pf('stats_weighted_mean', [1, 2], [-1, 2]), 1001);
    dcode(fn () => pf('stats_weighted_mean', [1, null], [1, -2], missing: 'drop'), 1001);
    dcode(fn () => pf('stats_weighted_var', [1, null], [1, -2], missing: 'drop'), 1001);
    dcode(fn () => pf('stats_var', [-1e200, 1e200]), 1003);
    dcode(fn () => pf('stats_cov', [-1e200, 1e200], [1, 2]), 1003);
    dcode(fn () => pf('stats_weighted_var', [-1e200, 1e200], [1, 1]), 1003);
    dcode(fn () => pf('stats_weighted_mean', [1, 2], [PHP_FLOAT_MAX, PHP_FLOAT_MAX]), 1003);
    dcode(fn () => pf('stats_mean', [PHP_FLOAT_MAX, PHP_FLOAT_MAX]), 1003);
    dcode(fn () => pf('stats_weighted_var', [1, 2], [.25, .5]), 1001);
    dcode(fn () => pf('stats_weighted_var', [1, 2], [0, 2], 'reliability'), 1001);
});

it('handles degenerate inputs', function () {
    expect(pf('stats_cor', [2, 2], [1, 3]))->toBeNan();
    expect(pf('stats_var', [7], 'population'))->toBe(0.0);
    expect(pf('stats_quantile_array', [1, 2], []))->toBe([]);
    expect(pf('stats_mean', ['later' => 3, 'earlier' => 1]))->toBe(2.0);
    expect(pf('stats_quantile_array', [3, 1, 2], ['high' => 1, 'low' => 0]))->toBe([3.0, 1.0]);
});

it('does not lose variability to large offsets', function () {
    $values = [1e9 + 4.0, 1e9 + 7.0, 1e9 + 13.0, 1e9 + 16.0];
    expect(pf('stats_var', $values))->toBe(30.0)
        ->and(pf('stats_mean', $values))->toBe(1e9 + 10.0);
    expect(pf('stats_cor', [1, 2, 3], [1, 1, 1]))->toBeNan();
});

it('exposes the global functions with matching signatures', function () {
    foreach (['stats_mean', 'stats_var', 'stats_sd', 'stats_median', 'stats_quantile', 'stats_quantile_array',
        'stats_cov', 'stats_cor', 'stats_weighted_mean', 'stats_weighted_var'] as $fn) {
        expect(function_exists($fn))->toBeTrue($fn);
        expect(\PhpStats\Polyfill\Polyfill::supports($fn))->toBeTrue($fn);
    }
});
