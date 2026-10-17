<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal\Bobyqa;

/**
 * Why the solver stopped (or why construction failed), the crate's Status.
 *
 * @internal
 */
enum Status
{
    case Converged;
    case TargetReached;
    case FtolReached;
    case MaxFunReached;
    case ModelDegenerate;
    case InvalidArgs;
    case AllocationFailed;

    /** The crate's Display text. */
    public function label(): string
    {
        return match ($this) {
            self::Converged => 'the trust-region radius reached rho_end, with any restart schedule settled',
            self::TargetReached => 'an evaluation reached f_target',
            self::FtolReached => 'the best f improved by less than the ftol tolerance over one rho stage',
            self::MaxFunReached => 'the max_fun evaluation budget was exhausted',
            self::ModelDegenerate => 'the interpolation model degenerated beyond rescue',
            self::InvalidArgs => 'invalid arguments: bad bounds, npt, sizes, or config',
            self::AllocationFailed => 'the solver workspace could not be allocated',
        };
    }

    /** The PRIMA exit code to Status mapping (status_from_info). */
    public static function fromInfo(int $info): self
    {
        return match ($info) {
            Util::SMALL_TR_RADIUS => self::Converged,
            Util::FTARGET_ACHIEVED => self::TargetReached,
            Util::FTOL_REACHED => self::FtolReached,
            Util::MAXFUN_REACHED, Util::MAXTR_REACHED => self::MaxFunReached,
            default => self::ModelDegenerate,
        };
    }
}
