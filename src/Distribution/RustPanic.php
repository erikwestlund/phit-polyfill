<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/**
 * A panic in the extension's Rust code, which PHP sees as an \Error with
 * code 0. Throwing a subclass lets the polyfill catch exactly the panics the
 * extension turns into NaN, and leave the others as they are.
 *
 * @internal
 */
final class RustPanic extends \Error {}
