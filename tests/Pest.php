<?php

declare(strict_types=1);

use PhpStats\Polyfill\Polyfill;

/*
 * The ported suites call the polyfill implementation directly through pf(),
 * so they test the polyfill whether or not the extension is loaded. The
 * Parity suite compares those implementations with the loaded extension.
 */

/** Decoded JSON fixture from tests/fixtures. */
function loadFixture(string $name): array
{
    static $cache = [];

    return $cache[$name] ??= json_decode(
        file_get_contents(__DIR__.'/fixtures/'.$name),
        true,
        flags: JSON_THROW_ON_ERROR,
    );
}

/** Call the polyfill implementation of a stats_* function (named arguments allowed). */
function pf(string $fn, mixed ...$arguments): mixed
{
    return Polyfill::implementation($fn)(...$arguments);
}

/** Call the native extension's function; only valid when it is loaded. */
function native(string $fn, mixed ...$arguments): mixed
{
    return $fn(...$arguments);
}

/** |actual - expected| <= abs + rel * |expected|, elementwise, with matching NaN/INF. */
function close(mixed $actual, mixed $expected, float $abs = 1e-9, float $rel = 1e-7, string $label = ''): bool
{
    if (is_array($expected)) {
        if (! is_array($actual) || count($actual) !== count($expected)) {
            return false;
        }
        foreach ($expected as $key => $value) {
            if (! array_key_exists($key, $actual) || ! close($actual[$key], $value, $abs, $rel)) {
                return false;
            }
        }

        return true;
    }
    if (! is_int($actual) && ! is_float($actual)) {
        return $actual === $expected;
    }
    if (is_string($expected)) {
        // JSON fixtures may spell non-finite values as strings.
        $expected = match (strtolower($expected)) {
            'nan' => NAN, 'inf', 'infinity' => INF, '-inf', '-infinity' => -INF,
            default => (float) $expected,
        };
    }
    if ($expected === null) {
        return is_nan((float) $actual);
    }
    $actual = (float) $actual;
    $expected = (float) $expected;
    if (is_nan($expected) || is_nan($actual)) {
        return is_nan($expected) && is_nan($actual);
    }
    if (is_infinite($expected) || is_infinite($actual)) {
        return $actual === $expected;
    }

    return abs($actual - $expected) <= $abs + $rel * abs($expected);
}

/** Collect PHP warnings raised while $callback runs. */
function warnings(callable $callback, mixed &$result = null): array
{
    $warnings = [];
    set_error_handler(function (int $severity, string $message) use (&$warnings): bool {
        $warnings[] = [$severity, $message];

        return true;
    });
    try {
        $result = $callback();
    } finally {
        restore_error_handler();
    }

    return $warnings;
}

/** The exception code thrown by $callback, or null. */
function thrownCode(callable $callback): ?int
{
    try {
        warnings($callback);
    } catch (\Throwable $e) {
        return $e->getCode();
    }

    return null;
}

expect()->extend('toBeClose', function (mixed $expected, float $abs = 1e-9, float $rel = 1e-7) {
    $ok = close($this->value, $expected, $abs, $rel);
    expect($ok)->toBeTrue(sprintf(
        "Values differ beyond abs %.1e + rel %.1e:\nactual:   %s\nexpected: %s",
        $abs,
        $rel,
        json_encode($this->value, JSON_PARTIAL_OUTPUT_ON_ERROR | JSON_PRESERVE_ZERO_FRACTION),
        json_encode($expected, JSON_PARTIAL_OUTPUT_ON_ERROR | JSON_PRESERVE_ZERO_FRACTION),
    ));

    return $this;
});
