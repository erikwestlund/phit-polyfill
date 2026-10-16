<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/** @internal */
final class Limits
{
    /** Largest binomial trial weight; larger-count likelihood arithmetic is not verified. */
    public const float MAX_TRIALS = 1000000.0;
}
