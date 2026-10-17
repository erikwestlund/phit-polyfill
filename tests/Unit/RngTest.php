<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\ChaCha20;

// Port of the extension's tests/rng.php against the R/NumPy ChaCha20 reference
// vectors in rng-python-r.json. Continuous draws within 1e-12 relative; counts exact.

$fixture = loadFixture('rng-python-r.json');

it('reproduces the reference sample vectors', function () use ($fixture) {
    $compared = 0;
    foreach ($fixture['cases'] as $c) {
        if ($c['distribution'] === 'uniform') {
            continue; // stats-core engine reference only.
        }
        $fn = 'Phit\dist_'.$c['distribution'].'_rng';
        $draws = pf($fn, ...['n' => $c['n'], 'seed' => $c['seed'], ...$c['params']]);
        expect(array_is_list($draws))->toBeTrue($fn)->and($draws)->toHaveCount($c['n']);
        expect(pf($fn, $c['n'], $c['seed'], ...$c['params']))->toBe($draws, "$fn repeatability");
        expect(pf($fn, 7, $c['seed'], ...$c['params']))->toBe(array_slice($draws, 0, 7), "$fn prefix stability");
        expect(pf($fn, 0, $c['seed'], ...$c['params']))->toBe([]);
        $discrete = in_array($c['distribution'], ['binomial', 'poisson', 'nbinom', 'nbinom_mu'], true);
        foreach ($draws as $i => $v) {
            $expected = $c['expected'][$i];
            expect(is_float($v) && is_finite($v))->toBeTrue("$fn numeric output");
            if ($discrete) {
                expect($v == $expected && $v >= 0 && floor($v) === $v)->toBeTrue("$fn count $i");
            } else {
                expect(abs($v - $expected) <= 1e-15 + 1e-12 * abs($expected))->toBeTrue("$fn: $v != $expected");
            }
            $compared++;
        }
        expect(thrownCode(fn () => pf($fn, -1, $c['seed'], ...$c['params'])))->toBe(1001);
        expect(thrownCode(fn () => pf($fn, 1, -1, ...$c['params'])))->toBe(1001);
        expect(thrownCode(fn () => pf($fn, 1000001, $c['seed'], ...$c['params'])))->toBe(1007);
    }
    expect($compared)->toBeGreaterThan(1000);
});

it('keeps every call on its own stream', function () {
    $first = pf('Phit\dist_normal_rng', 128, 42);
    pf('Phit\dist_gamma_rng', 128, 99, shape: 3);
    expect(pf('Phit\dist_normal_rng', 128, 42))->toBe($first)
        ->and(pf('Phit\dist_normal_rng', 128, 43))->not->toBe($first)
        ->and(pf('Phit\dist_normal_rng', 128, 42, mean: 0, sd: 1))->toBe($first)
        ->and(pf('Phit\dist_lognormal_rng', 5, 42))->toBe(pf('Phit\dist_lognormal_rng', 5, 42, meanlog: 0, sdlog: 1))
        ->and(pf('Phit\dist_gamma_rng', 5, 42, shape: 2))->toBe(pf('Phit\dist_gamma_rng', 5, 42, shape: 2, scale: 1))
        ->and(pf('Phit\dist_exponential_rng', 5, 42))->toBe(pf('Phit\dist_exponential_rng', 5, 42, rate: 1))
        ->and(pf('Phit\dist_logistic_rng', 5, 42))->toBe(pf('Phit\dist_logistic_rng', 5, 42, location: 0, scale: 1));
});

it('validates parameters even for an empty batch', function () {
    foreach ([
        fn () => pf('Phit\dist_normal_rng', 0, 42, sd: 0), fn () => pf('Phit\dist_t_rng', 0, 42, df: 0),
        fn () => pf('Phit\dist_beta_rng', 0, 42, alpha: -1, beta: 2), fn () => pf('Phit\dist_gamma_rng', 0, 42, shape: 0),
        fn () => pf('Phit\dist_binomial_rng', 0, 42, size: -1, prob: .5), fn () => pf('Phit\dist_poisson_rng', 0, 42, lambda: -1),
        fn () => pf('Phit\dist_nbinom_rng', 0, 42, size: 2, prob: 0), fn () => pf('Phit\dist_nbinom_mu_rng', 0, 42, size: 0, mu: 3),
    ] as $f) {
        expect(thrownCode($f))->toBe(1001);
    }
    foreach ([NAN, INF, -INF] as $v) {
        expect(thrownCode(fn () => pf('Phit\dist_normal_rng', 0, 42, mean: $v)))->toBe(1003)
            ->and(thrownCode(fn () => pf('Phit\dist_poisson_rng', 0, 42, lambda: $v)))->toBe(1003);
    }
    expect(thrownCode(fn () => pf('Phit\dist_binomial_rng', 0, 42, size: 2 ** 54, prob: .5)))->toBe(1007)
        ->and(thrownCode(fn () => pf('Phit\dist_normal_rng', 1, 42, mean: PHP_FLOAT_MAX, sd: PHP_FLOAT_MAX)))->toBe(1003);
});

it('has the expected moments', function () {
    $normal = pf('Phit\dist_normal_rng', 50000, 173, mean: 3, sd: 2);
    $poisson = pf('Phit\dist_poisson_rng', 50000, 174, lambda: 12);
    expect(abs(pf('Phit\stat_mean', $normal) - 3))->toBeLessThan(.06)
        ->and(abs(pf('Phit\stat_var', $normal, 'population') - 4))->toBeLessThan(.15)
        ->and(abs(pf('Phit\stat_mean', $poisson) - 12))->toBeLessThan(.1)
        ->and(abs(pf('Phit\stat_var', $poisson, 'population') - 12))->toBeLessThan(.4);
});

it('produces the same keystream with and without OpenSSL', function () {
    foreach ([0, 1, 42, PHP_INT_MAX, 2 ** 40 + 7] as $seed) {
        $fast = new ChaCha20($seed);
        $plain = ChaCha20::withoutOpenssl($seed);
        for ($i = 0; $i < 300; $i++) {
            expect($fast->open01())->toBe($plain->open01());
        }
    }
});
