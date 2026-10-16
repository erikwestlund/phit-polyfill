<?php

declare(strict_types=1);

use PhpStats\Polyfill\Polyfill;

// The global functions the polyfill defines must have the extension's exact
// signatures (tests/fixtures/php-api.json), so code runs unchanged on either.

$api = loadFixture('php-api.json');

it('declares every polyfilled function with the extension signature', function () use ($api) {
    foreach (Polyfill::functions() as $fn) {
        expect($api)->toHaveKey($fn);
        $spec = $api[$fn];
        $reflection = extension_loaded('php_stats')
            ? new ReflectionMethod(...Polyfill::implementation($fn))
            : new ReflectionFunction($fn);
        $parameters = $reflection->getParameters();
        expect($parameters)->toHaveCount(count($spec['parameters']), $fn);
        foreach ($spec['parameters'] as $i => $expected) {
            $actual = $parameters[$i];
            expect($actual->getName())->toBe($expected['name'], "$fn parameter $i")
                ->and((string) $actual->getType())->toBe($expected['type'], "$fn \${$expected['name']} type")
                ->and($actual->isOptional())->toBe($expected['optional'], "$fn \${$expected['name']} optional");
            if ($expected['optional']) {
                expect($actual->getDefaultValue())->toEqual($expected['default'], "$fn \${$expected['name']} default");
            }
        }
        expect((string) $reflection->getReturnType())->toBe($spec['return_type'], "$fn return type");
    }
});

it('defines the global functions only when the extension is absent', function () {
    foreach (Polyfill::functions() as $fn) {
        expect(function_exists($fn))->toBeTrue($fn);
        $function = new ReflectionFunction($fn);
        expect($function->isInternal())->toBe(extension_loaded('php_stats'), $fn);
    }
});

it('leaves extension-only functions undefined', function () use ($api) {
    $extensionOnly = array_diff(array_keys($api), Polyfill::functions());
    foreach ($extensionOnly as $fn) {
        expect(Polyfill::supports($fn))->toBeFalse($fn);
        if (! extension_loaded('php_stats')) {
            expect(function_exists($fn))->toBeFalse($fn);
        }
    }
    foreach (['stats_lmm', 'stats_lmm_fit', 'stats_logistic_mixed', 'stats_hello', 'stats_normal_rng', 'stats_draw_affine'] as $fn) {
        expect(Polyfill::supports($fn))->toBeFalse($fn);
    }
});

it('routes the global functions to the polyfill implementation', function () {
    if (extension_loaded('php_stats')) {
        $this->markTestSkipped('the extension provides the global functions');
    }
    $x = [[1.0], [2.0], [3.0], [4.0], [5.0]];
    $y = [2.1, 4.2, 5.8, 8.1, 10.2];
    expect(stats_ols($x, $y))->toBe(Polyfill::implementation('stats_ols')($x, $y))
        ->and(stats_mean([1, 2, 3]))->toBe(2.0);
});
