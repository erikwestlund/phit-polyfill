<?php

declare(strict_types=1);

namespace PhpStats\Polyfill;

/**
 * Exception codes shared with the php_stats extension. Both throw a plain
 * \Exception whose code identifies the category; messages are not part of
 * the contract.
 */
final class ErrorCode
{
    public const int INVALID_INPUT = 1001;
    public const int DIMENSION_MISMATCH = 1002;
    public const int NONFINITE = 1003;
    public const int RANK_DEFICIENCY = 1004;
    public const int FACTORIZATION_FAILURE = 1005;
    public const int NONCONVERGENCE = 1006;
    public const int UNSUPPORTED = 1007;
    public const int INTERNAL = 1008;

    public static function exception(int $code, string $message): \Exception
    {
        return new \Exception($message, $code);
    }
}
