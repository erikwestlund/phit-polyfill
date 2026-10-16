<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/**
 * An argmin `Error` (anyhow::Error): the message is its Display text.
 *
 * A cost callable throws it where the Rust CostFunction returns `Err(...)`,
 * e.g. `new NelderMeadError('Nelder-Mead evaluation budget exhausted')` for
 * `Error::msg(...)`. NelderMead also throws it for the crate's own errors
 * (invalid parameters, the "unreachable point" PotentialBug).
 *
 * @internal
 */
final class NelderMeadError extends \RuntimeException {}
