<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

/**
 * The argmin TerminationReason variants a Nelder-Mead run can end with; the
 * value is the variant's Display text. Interrupt, Timeout and SolverExit
 * cannot occur without ctrlc, a timeout or a solver that exits on its own.
 *
 * @internal
 */
enum NelderMeadTermination: string
{
    case MaxItersReached = 'Maximum number of iterations reached';
    case TargetCostReached = 'Target cost value reached';
    case SolverConverged = 'Solver converged';
}
