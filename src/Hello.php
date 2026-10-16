<?php

declare(strict_types=1);

namespace PhpStats\Polyfill;

/** stats_hello, the extension's smoke-test function. */
final class Hello
{
    /** stats_hello */
    public static function hello(): string
    {
        return 'Hello from Rust!';
    }
}
