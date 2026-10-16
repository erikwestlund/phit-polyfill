<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * The result of one Bobyqa::minimize call.
 *
 * @internal
 */
final readonly class Outcome
{
    public function __construct(
        /** Best objective value found; NAN when nothing was evaluated. */
        public float $f,
        public int $nEval,
        public Status $status,
    ) {}

    /** Whether the solver ever evaluated a genuinely finite objective value (f < FUNCMAX). */
    public function foundFinite(): bool
    {
        return $this->f < Util::FUNCMAX;
    }
}
