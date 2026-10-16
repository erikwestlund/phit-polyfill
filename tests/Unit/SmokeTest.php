<?php

use PhpStats\Polyfill\Polyfill;

it('reports unsupported functions', function () {
    expect(Polyfill::supports('stats_hello'))->toBeFalse()
        ->and(Polyfill::supports('stats_lmm'))->toBeFalse();
});
