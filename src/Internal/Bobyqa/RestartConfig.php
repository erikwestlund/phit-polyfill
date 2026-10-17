<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal\Bobyqa;

/**
 * The crate's RestartConfig. The defaults are its recommended schedule: one
 * restart, triggered once a cycle has spent an eighth of the budget that
 * remained when it started, settling below a 1e-6 relative improvement.
 *
 * @internal
 */
final class RestartConfig
{
    public int $maxRestarts = 1;

    public float $improveRelTol = 1e-6;

    public int $stallReductions = 0;

    public float $cycleBudgetFrac = 0.125;
}
