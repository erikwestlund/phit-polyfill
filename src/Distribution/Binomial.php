<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

use PhpStats\Polyfill\ErrorCode;

/** Binomial distribution: successes in `size` trials with success probability `prob`. */
final class Binomial
{
    /** Largest supported size: integers are exact in a double up to 2^53. */
    public const int MAX_SIZE = 1 << 53;

    private function __construct(
        private readonly float $n,
        private readonly float $p,
        private readonly float $q,
        private readonly float $lnP,
        private readonly float $lnQ,
    ) {}

    /** @internal */
    public static function create(int $size, float $prob): self
    {
        if ($size < 0) {
            throw Support::invalid('size must be nonnegative');
        }
        Support::finite('prob', $prob);
        if (! ($prob >= 0.0 && $prob <= 1.0)) {
            throw Support::invalid('prob must be in [0, 1]');
        }
        if ($size > self::MAX_SIZE) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'size must not exceed 2^53 (exact integer range of a double)');
        }

        return new self((float) $size, $prob, 1.0 - $prob, log($prob), log1p(-$prob));
    }

    /** stats_binomial_pmf */
    public static function pmf(float $k, int $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->pmfAt(...));
    }

    /** stats_binomial_pmf_array */
    public static function pmfArray(array $counts, int $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->pmfAt(...));
    }

    /** stats_binomial_logpmf */
    public static function logpmf(float $k, int $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->logpmfAt(...));
    }

    /** stats_binomial_logpmf_array */
    public static function logpmfArray(array $counts, int $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->logpmfAt(...));
    }

    /** stats_binomial_cdf */
    public static function cdf(float $k, int $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->cdfAt(...));
    }

    /** stats_binomial_cdf_array */
    public static function cdfArray(array $counts, int $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->cdfAt(...));
    }

    /** stats_binomial_logcdf */
    public static function logcdf(float $k, int $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->logcdfAt(...));
    }

    /** stats_binomial_logcdf_array */
    public static function logcdfArray(array $counts, int $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->logcdfAt(...));
    }

    /** stats_binomial_sf */
    public static function sf(float $k, int $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->sfAt(...));
    }

    /** stats_binomial_sf_array */
    public static function sfArray(array $counts, int $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->sfAt(...));
    }

    /** stats_binomial_logsf */
    public static function logsf(float $k, int $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->logsfAt(...));
    }

    /** stats_binomial_logsf_array */
    public static function logsfArray(array $counts, int $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->logsfAt(...));
    }

    /** stats_binomial_quantile */
    public static function quantile(float $p, int $size, float $prob): float
    {
        return Support::probability($p, self::create($size, $prob)->quantileAt(...));
    }

    /** stats_binomial_quantile_array */
    public static function quantileArray(array $probabilities, int $size, float $prob): array
    {
        return Support::probabilities($probabilities, self::create($size, $prob)->quantileAt(...));
    }

    private function exact(float $k): ?array
    {
        if ($k < 0.0) {
            return [0.0, 1.0];
        }
        if ($k >= $this->n || $this->p == 0.0) {
            return [1.0, 0.0];
        }
        if ($this->q == 0.0) {
            return [0.0, 1.0];
        }

        return null;
    }

    private function tails(float $k): array
    {
        if (($exact = $this->exact($k)) !== null) {
            return $exact;
        }
        $k = floor($k);
        [$lower, $upper] = Cdflib::betaInc($this->n - $k, $k + 1.0, $this->q, $this->p);

        return [
            Special::refine($lower, fn (): float => $this->lnLower($k)),
            Special::refine($upper, fn (): float => $this->lnUpper($k)),
        ];
    }

    private function lnLower(float $k): float
    {
        return Special::lnBetaLowerFraction($this->n - $k, $k + 1.0, $this->q, $this->lnQ, $this->lnP);
    }

    private function lnUpper(float $k): float
    {
        return Special::lnBetaLowerFraction($k + 1.0, $this->n - $k, $this->p, $this->lnP, $this->lnQ);
    }

    /** @internal */
    public function pmfAt(float $k): float
    {
        return exp($this->logpmfAt($k));
    }

    /** @internal */
    public function logpmfAt(float $k): float
    {
        if (! Quantile::isCount($k) || $k > $this->n) {
            return -INF;
        }

        return Special::lnBinomialRaw($k, $this->n, $this->p, $this->q);
    }

    /** @internal */
    public function cdfAt(float $k): float
    {
        return $this->tails($k)[0];
    }

    /** @internal */
    public function sfAt(float $k): float
    {
        return $this->tails($k)[1];
    }

    /** @internal */
    public function logcdfAt(float $k): float
    {
        if (($exact = $this->exact($k)) !== null) {
            return log($exact[0]);
        }
        [$lower, $upper] = $this->tails($k);

        return Special::logTail($lower, $upper, fn (): float => $this->lnLower(floor($k)));
    }

    /** @internal */
    public function logsfAt(float $k): float
    {
        if (($exact = $this->exact($k)) !== null) {
            return log($exact[1]);
        }
        [$lower, $upper] = $this->tails($k);

        return Special::logTail($upper, $lower, fn (): float => $this->lnUpper(floor($k)));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        if ($this->p == 0.0 || $this->n == 0.0) {
            return 0.0;
        }
        if ($this->q == 0.0) {
            return $p == 0.0 ? 0.0 : $this->n;
        }
        $sd = sqrt($this->n * $this->p * $this->q);
        $first = Quantile::cornishFisher($p, $this->n * $this->p, $sd, fdiv($this->q - $this->p, $sd));

        return Quantile::discrete($p, $first, $this->n, $this->cdfAt(...), $this->sfAt(...), $this->logcdfAt(...));
    }
}
