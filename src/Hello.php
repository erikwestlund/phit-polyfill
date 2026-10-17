<?php

declare(strict_types=1);

namespace Phit\Polyfill;

/** Phit\hello, the extension's smoke-test function. */
final class Hello
{
    /** Phit\hello */
    public static function hello(): string
    {
        return 'Hello from Rust!';
    }
}
