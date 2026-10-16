<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * The crate's RestartState: bookkeeping threaded through the solve when
 * Config::$restart is set.
 *
 * @internal
 */
final class RestartState
{
    public ?float $lastFopt = null;

    public int $stallCount = 0;

    public ?float $stallFopt = null;

    public int $nfCycleStart = 0;

    /** @var list<int> cumulative evaluation count at each restart */
    public array $cycleBoundaries = [];

    public function __construct(public readonly RestartConfig $config) {}

    public function shouldRestart(float $fopt, int $nf, int $maxfun): bool
    {
        if (count($this->cycleBoundaries) >= $this->config->maxRestarts) {
            return false;
        }
        if ($nf >= $maxfun) {
            return false;
        }
        if ($this->lastFopt === null) {
            return true;
        }
        $improved = $this->lastFopt - $fopt;

        return $improved >= $this->config->improveRelTol * Linalg::max(abs($fopt), 1.0);
    }

    public function noteReduction(float $fopt): bool
    {
        if ($this->config->stallReductions === 0) {
            return false;
        }
        if ($this->stallFopt === null) {
            $this->stallCount = 0;
        } else {
            $improved = $this->stallFopt - $fopt;
            if ($improved >= $this->config->improveRelTol * Linalg::max(abs($fopt), 1.0)) {
                $this->stallCount = 0;
            } else {
                $this->stallCount++;
            }
        }
        $this->stallFopt = $fopt;

        return $this->stallCount >= $this->config->stallReductions;
    }

    public function overCycleBudget(int $nf, int $maxfun): bool
    {
        $frac = $this->config->cycleBudgetFrac;

        return $frac > 0.0
            && (float) ($nf - $this->nfCycleStart) >= $frac * (float) ($maxfun - $this->nfCycleStart);
    }

    public function onRestart(float $fopt, int $nf): void
    {
        $this->lastFopt = $fopt;
        $this->stallCount = 0;
        $this->stallFopt = $fopt;
        $this->nfCycleStart = $nf;
        $this->cycleBoundaries[] = $nf;
    }
}
