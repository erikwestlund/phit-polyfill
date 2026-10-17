<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

/**
 * The final IterState of a Nelder-Mead run that returned Ok.
 *
 * @internal
 */
final readonly class NelderMeadResult
{
    /**
     * @param  list<float>|null  $bestParam  get_best_param(); null when no cost ever beat +INF (all NaN)
     * @param  float  $bestCost  get_best_cost()
     * @param  int  $iterations  get_iter()
     */
    public function __construct(
        public NelderMeadTermination $termination,
        public ?array $bestParam,
        public float $bestCost,
        public int $iterations,
    ) {}
}
