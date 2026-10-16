<?php

// Polyfill vs native php_stats timings for stats_ols and stats_logistic.
//
// Run with the extension and JIT enabled, for example:
//   php85 -d extension=.ext/libphp_stats.dylib -d opcache.enable_cli=1 \
//         -d opcache.jit=tracing -d opcache.jit_buffer_size=64M benchmarks/run.php
// Options: --sizes=10000,100000 --p=10 --repeat=3 --seed=1

declare(strict_types=1);

require __DIR__.'/../vendor/autoload.php';

use PhpStats\Polyfill\Logistic;
use PhpStats\Polyfill\Ols;

$options = getopt('', ['sizes::', 'p::', 'repeat::', 'seed::']);
$sizes = array_map(intval(...), explode(',', $options['sizes'] ?? '10000,100000'));
$p = (int) ($options['p'] ?? 10);
$repeat = (int) ($options['repeat'] ?? 3);
$seed = (int) ($options['seed'] ?? 1);
$native = extension_loaded('php_stats');
$jit = function_exists('opcache_get_status') && (opcache_get_status(false)['jit']['on'] ?? false);

function gaussian(): float
{
    $u = (mt_rand() + 1) / (mt_getrandmax() + 2.0);
    $v = mt_rand() / (mt_getrandmax() + 1.0);

    return sqrt(-2.0 * log($u)) * cos(2.0 * M_PI * $v);
}

/** Median wall time in milliseconds. */
function time_ms(callable $f, int $repeat): float
{
    $times = [];
    for ($i = 0; $i < $repeat; $i++) {
        $start = hrtime(true);
        $f();
        $times[] = (hrtime(true) - $start) / 1e6;
    }
    sort($times);

    return $times[intdiv(count($times), 2)];
}

printf("PHP %s, JIT %s, extension %s, p = %d predictors + intercept, median of %d\n\n",
    PHP_VERSION, $jit ? 'on' : 'off', $native ? 'loaded' : 'not loaded', $p, $repeat);
printf("%-10s %9s %14s %14s %8s\n", 'model', 'n', 'polyfill ms', 'extension ms', 'ratio');

foreach ($sizes as $n) {
    mt_srand($seed + $n);
    $x = [];
    $continuous = [];
    $binary = [];
    for ($i = 0; $i < $n; $i++) {
        $row = [];
        $eta = 0.3;
        for ($j = 0; $j < $p; $j++) {
            $value = gaussian();
            $row[] = $value;
            $eta += $value * (($j % 3) - 1) * 0.2;
        }
        $x[] = $row;
        $continuous[] = $eta + gaussian();
        $binary[] = (mt_rand() / mt_getrandmax()) < 1 / (1 + exp(-$eta)) ? 1 : 0;
    }
    foreach ([['ols', Ols::fit(...), 'stats_ols', $continuous], ['logistic', Logistic::fit(...), 'stats_logistic', $binary]] as [$name, $polyfill, $fn, $y]) {
        try {
            $ours = time_ms(fn () => $polyfill($x, $y), $repeat);
            $theirs = $native ? time_ms(fn () => $fn($x, $y), $repeat) : NAN;
        } catch (Exception $e) {
            printf("%-10s %9d   failed with code %d: %s\n", $name, $n, $e->getCode(), $e->getMessage());

            continue;
        }
        printf("%-10s %9d %14.1f %14.1f %8s\n", $name, $n, $ours, $theirs, $native ? sprintf('%.1fx', $ours / $theirs) : '-');
    }
}
