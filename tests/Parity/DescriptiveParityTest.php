<?php

declare(strict_types=1);

// Compares the polyfill with the loaded extension; skipped without it.
const PARITY_REL = 1e-13;

/** Largest relative difference seen so far (pass a reset to read it). */
function parityMax(?float $seen = null): float
{
    static $max = 0.0;
    if ($seen !== null) {
        $max = max($max, $seen);
    }

    return $max;
}

function parityCompare(string $fn, array $args): void
{
    $pf = $nat = $pfCode = $natCode = null;
    try {
        $pf = pf($fn, ...$args);
    } catch (\Throwable $e) {
        $pfCode = $e instanceof \Exception ? $e->getCode() : get_class($e);
    }
    try {
        $nat = native($fn, ...$args);
    } catch (\Throwable $e) {
        $natCode = $e instanceof \Exception ? $e->getCode() : get_class($e);
    }
    $label = fn () => $fn.' '.json_encode($args, JSON_PARTIAL_OUTPUT_ON_ERROR).' pf='
        .json_encode($pf, JSON_PARTIAL_OUTPUT_ON_ERROR).' native='.json_encode($nat, JSON_PARTIAL_OUTPUT_ON_ERROR);
    if ($pfCode !== $natCode) {
        throw new \RuntimeException("exception mismatch ($pfCode vs $natCode): ".$label());
    }
    expect($pfCode)->toBe($natCode);
    if ($pfCode !== null) {
        return;
    }
    if (! close($pf, $nat, 0.0, PARITY_REL)) {
        throw new \RuntimeException('value mismatch: '.$label());
    }
    foreach (is_array($pf) ? $pf : [$pf] as $k => $v) {
        $n = is_array($nat) ? $nat[$k] : $nat;
        if (is_finite($v) && is_finite($n) && $n != 0.0) {
            parityMax(abs($v - $n) / abs($n));
        }
    }
}

function parityVector(int $n, string $kind): array
{
    $v = [];
    for ($i = 0; $i < $n; $i++) {
        $x = match ($kind) {
            'offset' => 1e9 + mt_rand(0, 1000000) / 1000.0,
            'huge' => 1e12 + mt_rand(0, 100) / 7.0,
            'int' => mt_rand(-50, 50),
            'tiny' => (mt_rand(0, 1000000) - 500000) * 1e-150,
            default => (mt_rand(0, 2000000) - 1000000) / 1000.0,
        };
        $v[] = $x;
    }

    return $v;
}

beforeEach(function () {
    if (! extension_loaded('php_stats')) {
        $this->markTestSkipped('php_stats extension not loaded');
    }
});

it('agrees on fixture inputs', function () {
    $fixture = loadFixture('descriptive-r.json');
    foreach ($fixture['descriptive'] as $c) {
        foreach (['mean', 'median'] as $f) {
            parityCompare("stats_$f", [$c['values']]);
        }
        foreach (['var', 'sd'] as $f) {
            parityCompare("stats_$f", [$c['values']]);
            parityCompare("stats_$f", [$c['values'], 'population']);
        }
        parityCompare('stats_quantile_array', [$c['values'], $c['quantile_probabilities']]);
        foreach ($c['quantile_probabilities'] as $p) {
            parityCompare('stats_quantile', [$c['values'], $p]);
        }
    }
    foreach ($fixture['pairs'] as $c) {
        parityCompare('stats_cov', [$c['x'], $c['y']]);
        parityCompare('stats_cov', [$c['x'], $c['y'], 'population']);
        parityCompare('stats_cor', [$c['x'], $c['y']]);
    }
    foreach ($fixture['weighted'] as $c) {
        parityCompare('stats_weighted_mean', [$c['values'], $c['weights']]);
        foreach (['frequency', 'reliability', 'population'] as $t) {
            parityCompare('stats_weighted_var', [$c['values'], $c['weights'], $t]);
        }
    }
});

it('agrees on seeded random inputs', function () {
    mt_srand(20261015);
    foreach (['plain', 'offset', 'huge', 'int', 'tiny'] as $kind) {
        foreach ([1, 2, 3, 7, 50, 500] as $n) {
            foreach (range(1, 4) as $rep) {
                $x = parityVector($n, $kind);
                $y = parityVector($n, $kind === 'int' ? 'plain' : $kind);
                $w = array_map(fn () => mt_rand(0, 40) / 8.0 + ($rep === 1 ? 1 : 0), range(1, $n));
                $p = array_map(fn () => mt_rand(0, 1000) / 1000.0, range(1, 5));
                foreach (['mean', 'median'] as $f) {
                    parityCompare("stats_$f", [$x]);
                }
                foreach (['var', 'sd'] as $f) {
                    parityCompare("stats_$f", [$x]);
                    parityCompare("stats_$f", [$x, 'population']);
                }
                parityCompare('stats_quantile', [$x, $p[0]]);
                parityCompare('stats_quantile_array', [$x, [...$p, 0, 1, 0.5]]);
                parityCompare('stats_cov', [$x, $y]);
                parityCompare('stats_cov', [$x, $y, 'population']);
                parityCompare('stats_cor', [$x, $y]);
                parityCompare('stats_weighted_mean', [$x, $w]);
                foreach (['frequency', 'reliability', 'population'] as $t) {
                    parityCompare('stats_weighted_var', [$x, $w, $t]);
                }
                // Missing values, dropped and raised.
                $xm = $x;
                $ym = $y;
                $wm = $w;
                $xm[mt_rand(0, $n - 1)] = null;
                $ym[mt_rand(0, $n - 1)] = NAN;
                $wm[mt_rand(0, $n - 1)] = mt_rand(0, 1) ? null : NAN;
                foreach (['drop', 'raise'] as $m) {
                    parityCompare('stats_mean', [$xm, $m]);
                    parityCompare('stats_median', [$xm, $m]);
                    parityCompare('stats_var', [$xm, 'sample', $m]);
                    parityCompare('stats_sd', [$xm, 'population', $m]);
                    parityCompare('stats_quantile', [$xm, 0.3, $m]);
                    parityCompare('stats_quantile_array', [$xm, [0.1, 0.9], $m]);
                    parityCompare('stats_cov', [$xm, $ym, 'sample', $m]);
                    parityCompare('stats_cor', [$xm, $ym, $m]);
                    parityCompare('stats_weighted_mean', [$xm, $wm, $m]);
                    parityCompare('stats_weighted_var', [$xm, $wm, 'reliability', $m]);
                }
            }
        }
    }
});

it('agrees on invalid inputs and options', function () {
    $inputs = [[], [1], [1, 2], [1, 2, 3, 4], [null], [NAN, NAN], [1, null, 3], [1, INF], [-INF, null], [INF, null],
        [null, INF], [1, 'a'], [1, true], [1, [2]], ['1', 2], [1, new stdClass()], [PHP_FLOAT_MAX, PHP_FLOAT_MAX],
        [-1e200, 1e200], [1e308, -1e308, 1e308], ['a' => 3, 'b' => 1], [5 => 1.5, 2 => 2]];
    $options = ['raise', 'drop', 'ignore', ''];
    $types = ['sample', 'population', 'frequency', 'reliability', 'bogus'];
    foreach ($inputs as $a) {
        foreach ($options as $m) {
            parityCompare('stats_mean', [$a, $m]);
            parityCompare('stats_median', [$a, $m]);
            parityCompare('stats_quantile', [$a, 0.5, $m]);
            parityCompare('stats_quantile_array', [$a, [0.5], $m]);
            foreach ($types as $t) {
                parityCompare('stats_var', [$a, $t, $m]);
                parityCompare('stats_sd', [$a, $t, $m]);
            }
        }
        foreach ($inputs as $b) {
            foreach (['raise', 'drop', 'ignore'] as $m) {
                parityCompare('stats_cor', [$a, $b, $m]);
                parityCompare('stats_weighted_mean', [$a, $b, $m]);
                foreach ($types as $t) {
                    parityCompare('stats_cov', [$a, $b, $t, $m]);
                    parityCompare('stats_weighted_var', [$a, $b, $t, $m]);
                }
            }
        }
        foreach ([[1], [null], [INF], ['x'], [NAN], [-0.5], [1.5], [], [0, 1, 0.5], [true], ['a' => 1.0, 'b' => 0]] as $p) {
            parityCompare('stats_quantile_array', [$a, $p]);
        }
        foreach ([0, 1, 0.5, -0.1, 1.1, NAN, INF, -INF] as $p) {
            parityCompare('stats_quantile', [$a, $p]);
        }
    }
    foreach ([[1, 2], [0, 0], [-1, 2], [1, null], [1, NAN], [1, INF], [PHP_FLOAT_MAX, PHP_FLOAT_MAX], [0.25, 0.5], [0, 2]] as $w) {
        foreach ([[1, 2], [1, null], [-1e200, 1e200]] as $v) {
            foreach (['frequency', 'reliability', 'population', 'x'] as $t) {
                foreach (['raise', 'drop'] as $m) {
                    parityCompare('stats_weighted_var', [$v, $w, $t, $m]);
                    parityCompare('stats_weighted_mean', [$v, $w, $m]);
                }
            }
        }
    }
    parityCompare('stats_cor', [[2, 2], [1, 3]]);
    parityCompare('stats_cor', [[1, 2, 3], [1, 2]]);
});

it('reports the largest relative difference', function () {
    if (getenv('PARITY_REPORT')) {
        file_put_contents('php://stderr', sprintf("max relative difference: %.3e\n", parityMax()));
    }
    expect(parityMax())->toBeLessThanOrEqual(PARITY_REL);
});
