<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal\Bobyqa;

/**
 * The crate's Config: PRIMA's defaults for an n-dimensional problem, each field
 * public to override before constructing the solver (which validates them).
 *
 * @internal
 */
final class Config
{
    /** Interpolation points, in n + 2 ..= (n + 1)(n + 2) / 2. */
    public int $npt;

    /** Initial trust-region radius, positive and finite. */
    public float $rhoBegin = 1.0;

    /** Final trust-region radius, in (0, rhoBegin]. */
    public float $rhoEnd = 1.0e-6;

    /** Objective-evaluation budget; must exceed npt. */
    public int $maxFun;

    /** Stop once an evaluation reaches f <= fTarget; -INF disables it. */
    public float $fTarget = -INF;

    public ?RestartConfig $restart = null;

    public ?float $ftolRel = null;

    public ?float $ftolAbs = null;

    /** false applies the crate's documented deviation on a fully determined model. */
    public bool $primaParity = true;

    public function __construct(int $n)
    {
        $this->npt = 2 * $n + 1;
        $this->maxFun = 500 * $n;
    }
}
