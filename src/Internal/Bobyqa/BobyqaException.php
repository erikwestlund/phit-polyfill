<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * A rejected solver configuration; the message is the crate's Display text of the status.
 *
 * @internal
 */
final class BobyqaException extends \InvalidArgumentException
{
    public function __construct(public readonly Status $status)
    {
        parent::__construct($status->label());
    }
}
