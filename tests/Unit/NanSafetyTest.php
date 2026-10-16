<?php

declare(strict_types=1);

// OPcache's CFG pass (optimization pass 5) rewrites a negated float comparison
// such as !($x > 0.0) into the inverted comparison ($x <= 0.0). The two differ
// for NaN, so with opcache enabled the negation evaluates to false where PHP's
// semantics (and Rust's) give true. See docs/notes/opcache-nan.md. Source code
// must therefore test NaN explicitly instead of relying on a negated comparison.

it('does not negate comparisons in the source', function () {
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
            // A `!` applied to a parenthesised expression containing <, >, <=, >= or ==.
            if (preg_match('/!\s*\((?:[^()]|\([^()]*\))*(?<![-=!])(?:<|>|<=|>=|==)(?![>=])/', $line)) {
                $offenders[] = $file->getFilename().':'.($number + 1).': '.trim($line);
            }
        }
    }
    expect($offenders)->toBe([]);
});

it('keeps NaN out of range checks under opcache', function () {
    expect(thrownCode(fn () => pf('stats_binomial_cdf', 1.0, 10, NAN)))->not->toBeNull()
        ->and(thrownCode(fn () => pf('stats_nbinom_cdf', 1.0, 2.0, NAN)))->not->toBeNull()
        ->and(pf('stats_beta_logpdf', 1.5, 2.0, 3.0))->toBe(-INF);
});
