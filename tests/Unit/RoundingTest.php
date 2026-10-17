<?php

declare(strict_types=1);

use PhpStats\Polyfill\Distribution\Support;

// PHP 8.3's round() pre-rounds its argument to 15 significant digits, so
// round(0.49999999999999994) is 1.0 there and 0.0 on PHP 8.4 and Rust. Source
// code rounds through Support::round, which is exact on every PHP version.

it('does not call round() in the source', function () {
    $offenders = [];
    $files = new RecursiveIteratorIterator(new RecursiveDirectoryIterator(dirname(__DIR__, 2).'/src'));
    foreach ($files as $file) {
        if ($file->getExtension() !== 'php') {
            continue;
        }
        foreach (file($file->getPathname()) as $number => $line) {
            if (preg_match('/^\s*(\/\/|\*)/', $line)) {
                continue;
            }
            if (preg_match('/(?<![\w:>$])\\\\?round\s*\(/', $line) && ! str_contains($line, 'function round(')) {
                $offenders[] = $file->getFilename().':'.($number + 1).': '.trim($line);
            }
        }
    }
    expect($offenders)->toBe([]);
});

it('rounds like f64::round', function () {
    $bits = static fn (float $v): string => bin2hex(pack('E', $v));
    $cases = [
        [0.49999999999999994, 0.0],
        [-0.49999999999999994, -0.0],
        [1.4999999999999998, 1.0],
        [2.4999999999999996, 2.0],
        [0.5, 1.0],
        [-0.5, -1.0],
        [2.5, 3.0],
        [-2.5, -3.0],
        [4503599627370495.5, 4503599627370496.0],
        [4503599627370497.0, 4503599627370497.0],
        [-0.0, -0.0],
        [INF, INF],
        [-INF, -INF],
    ];
    foreach ($cases as [$in, $out]) {
        expect($bits(Support::round($in)))->toBe($bits($out), (string) $in);
    }
    expect(is_nan(Support::round(NAN)))->toBeTrue();
});
