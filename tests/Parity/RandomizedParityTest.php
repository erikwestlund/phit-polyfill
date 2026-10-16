<?php

declare(strict_types=1);

// Randomized differential tests: seeded random inputs (including invalid and
// non-finite ones) go to the polyfill and to the loaded extension, and the
// results must be identical bit for bit, or both must throw the same exception
// class, code, and message. Skipped without the extension.
//
// PARITY_SEED picks the generator seed (default 20261016) and
// PARITY_ITERATIONS the rounds per function group (default 60).

use PhpStats\Polyfill\Polyfill;

if (! function_exists('randomParityOutcome')) {
    function randomParityOutcome(callable $call): array
    {
        try {
            return ['value', $call()];
        } catch (\Throwable $e) {
            return ['error', get_class($e), $e->getCode(), $e->getMessage()];
        }
    }

    function randomParityCompare(string $fn, array $arguments): void
    {
        $native = randomParityOutcome(fn () => $fn(...$arguments));
        $polyfill = randomParityOutcome(fn () => Polyfill::implementation($fn)(...$arguments));
        if (serialize($native) !== serialize($polyfill)) {
            $show = fn (mixed $v) => substr(json_encode($v, JSON_PARTIAL_OUTPUT_ON_ERROR | JSON_PRESERVE_ZERO_FRACTION), 0, 400);
            throw new \RuntimeException("$fn differs for ".$show($arguments)."\nnative:   ".$show($native)."\npolyfill: ".$show($polyfill));
        }
    }

    function randomUniform(float $lo, float $hi): float
    {
        return $lo + ($hi - $lo) * mt_rand() / mt_getrandmax();
    }

    function randomParityIterations(): int
    {
        return (int) (getenv('PARITY_ITERATIONS') ?: 60);
    }

    function randomParitySeed(string $group): void
    {
        mt_srand(crc32($group) ^ (int) (getenv('PARITY_SEED') ?: 20261016));
    }

    /** A covariance matrix A A' + D with a random positive diagonal. */
    function randomCovariance(int $p): array
    {
        $a = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $a[$i][$j] = randomUniform(-1, 1);
            }
        }
        $c = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $s = 0.0;
                for ($k = 0; $k < $p; $k++) {
                    $s += $a[$i][$k] * $a[$j][$k];
                }
                $c[$i][$j] = $s + ($i === $j ? randomUniform(0, 0.5) : 0.0);
            }
        }

        return $c;
    }
}

beforeEach(function () {
    if (! extension_loaded('php_stats')) {
        $this->markTestSkipped('php_stats extension not loaded');
    }
});

it('samples exactly like the extension', function () {
    randomParitySeed('rng');
    $families = [
        'normal' => fn () => [randomUniform(-5, 5), randomUniform(0.01, 10)],
        't' => fn () => [mt_rand(0, 3) > 0 ? randomUniform(0.1, 50) : randomUniform(1e-3, 1e6)],
        'chisquare' => fn () => [randomUniform(0.05, 100)],
        'f' => fn () => [randomUniform(0.1, 50), randomUniform(0.1, 50)],
        'exponential' => fn () => [randomUniform(0.01, 10)],
        'logistic' => fn () => [randomUniform(-5, 5), randomUniform(0.01, 10)],
        'lognormal' => fn () => [randomUniform(-3, 3), randomUniform(0.01, 3)],
        'gamma' => fn () => [randomUniform(0.01, 50), randomUniform(0.01, 10)],
        'beta' => fn () => [randomUniform(0.01, 20), randomUniform(0.01, 20)],
        'binomial' => fn () => [mt_rand(0, 500), randomUniform(0, 1)],
        'poisson' => fn () => [randomUniform(0, 300)],
        'nbinom' => fn () => [randomUniform(0.01, 50), randomUniform(0.001, 1)],
        'nbinom_mu' => fn () => [randomUniform(0.01, 50), randomUniform(0, 100)],
    ];
    $calls = 0;
    for ($round = 0; $round < randomParityIterations(); $round++) {
        foreach ($families as $family => $parameters) {
            $params = $parameters();
            if (mt_rand(0, 20) === 0) {
                $k = array_rand($params);
                // Integer parameters stay integers: the extension rejects floats
                // for them where PHP would coerce (a documented difference).
                $params[$k] = is_int($params[$k])
                    ? [-1, 0, PHP_INT_MAX][mt_rand(0, 2)]
                    : [-1.0, 0.0, NAN, INF, 1e308][mt_rand(0, 4)];
            }
            $n = [0, 1, 7, 64, 333][mt_rand(0, 4)];
            if (mt_rand(0, 40) === 0) {
                $n = [-1, 1000001][mt_rand(0, 1)];
            }
            $seed = mt_rand(0, 3) > 0 ? mt_rand(0, PHP_INT_MAX) : [0, -1, 1][mt_rand(0, 2)];
            randomParityCompare("stats_{$family}_rng", [$n, $seed, ...$params]);
            $calls++;
        }
    }
    expect($calls)->toBeGreaterThan(0);
});

it('computes draw kernels exactly like the extension', function () {
    randomParitySeed('draws');
    for ($round = 0; $round < randomParityIterations() * 3; $round++) {
        $n = mt_rand(1, 12);
        $p = mt_rand(1, 9);
        $scale = 10 ** mt_rand(-3, 6);
        $offset = mt_rand(0, 3) > 0 ? 0 : 10 ** mt_rand(2, 9);
        $draws = [];
        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $draws[$i][$j] = $offset + $scale * randomUniform(-1, 1);
            }
        }
        $weights = mt_rand(0, 2) > 0 ? null : array_map(fn () => mt_rand(0, 4) > 0 ? randomUniform(0, 5) : 0.0, range(1, $n));
        randomParityCompare('stats_draw_moments', [$draws, $weights, ['population', 'unbiased'][mt_rand(0, 1)]]);

        randomParityCompare('stats_draw_covariance_map', [randomCovariance($p), randomCovariance(mt_rand(0, 8) > 0 ? $p : $p + 1)]);

        $matrix = [];
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $matrix[$i][$j] = randomUniform(-2, 2) * 10 ** mt_rand(-2, 2);
            }
        }
        if ($p > 1 && mt_rand(0, 10) === 0) {
            $matrix[1] = $matrix[0];
        }
        $shift = array_map(fn () => randomUniform(-10, 10), range(1, $p));
        // More than five rows and columns takes the extension's fused kernel.
        $rows = mt_rand(0, 1) === 1 ? mt_rand(6, 40) : mt_rand(1, 5);
        $input = [];
        for ($i = 0; $i < $rows; $i++) {
            for ($j = 0; $j < $p; $j++) {
                $input[$i][$j] = randomUniform(-3, 3);
            }
        }
        randomParityCompare('stats_draw_affine', [$input, $matrix, $shift]);
    }
    randomParityCompare('stats_hello', []);
    expect(true)->toBeTrue();
});
