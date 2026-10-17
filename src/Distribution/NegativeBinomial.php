<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/**
 * Negative binomial distribution: failures before the `size`-th success,
 * P(X = k) = Γ(k + size)/(Γ(size) k!) prob^size (1 - prob)^k.
 * NegativeBinomialMu gives the same distribution parameterized by its mean.
 */
final class NegativeBinomial
{
    private function __construct(
        private readonly float $size,
        private readonly float $p,
        private readonly float $q,
        private readonly float $lnP,
        private readonly float $lnQ,
    ) {}

    /** @internal */
    public static function create(float $size, float $prob): self
    {
        Support::positive('size', $size);
        Support::finite('prob', $prob);
        if (\is_nan($prob) || $prob <= 0.0 || $prob > 1.0) {
            throw Support::invalid('prob must be in (0, 1]');
        }

        return new self($size, $prob, 1.0 - $prob, \log($prob), \log1p(-$prob));
    }

    /** @internal Mean parameterization: Var = mu + mu²/size. */
    public static function withMean(float $size, float $mu): self
    {
        Support::positive('size', $size);
        Support::finite('mu', $mu);
        if ($mu < 0.0) {
            throw Support::invalid('mu must be nonnegative');
        }
        // p = size/(size + mu) and q = mu/(size + mu), each without cancellation;
        // an overflowing sum must not become p = q = 0 (a false point mass).
        $total = $size + $mu;
        if (\is_finite($total)) {
            $p = $size / $total;
            $q = $mu / $total;
        } elseif ($size >= $mu) {
            $ratio = $mu / $size;
            $p = 1.0 / (1.0 + $ratio);
            $q = $ratio / (1.0 + $ratio);
        } else {
            $ratio = $size / $mu;
            $p = $ratio / (1.0 + $ratio);
            $q = 1.0 / (1.0 + $ratio);
        }

        return new self($size, $p, $q, -\log1p($mu / $size), $mu == 0.0 ? -INF : -\log1p($size / $mu));
    }

    /** stats_nbinom_pmf */
    public static function pmf(float $k, float $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->pmfAt(...));
    }

    /** stats_nbinom_pmf_array */
    public static function pmfArray(array $counts, float $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->pmfAt(...));
    }

    /** stats_nbinom_logpmf */
    public static function logpmf(float $k, float $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->logpmfAt(...));
    }

    /** stats_nbinom_logpmf_array */
    public static function logpmfArray(array $counts, float $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->logpmfAt(...));
    }

    /** stats_nbinom_cdf */
    public static function cdf(float $k, float $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->cdfAt(...));
    }

    /** stats_nbinom_cdf_array */
    public static function cdfArray(array $counts, float $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->cdfAt(...));
    }

    /** stats_nbinom_logcdf */
    public static function logcdf(float $k, float $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->logcdfAt(...));
    }

    /** stats_nbinom_logcdf_array */
    public static function logcdfArray(array $counts, float $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->logcdfAt(...));
    }

    /** stats_nbinom_sf */
    public static function sf(float $k, float $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->sfAt(...));
    }

    /** stats_nbinom_sf_array */
    public static function sfArray(array $counts, float $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->sfAt(...));
    }

    /** stats_nbinom_logsf */
    public static function logsf(float $k, float $size, float $prob): float
    {
        return Support::point($k, self::create($size, $prob)->logsfAt(...));
    }

    /** stats_nbinom_logsf_array */
    public static function logsfArray(array $counts, float $size, float $prob): array
    {
        return Support::points($counts, self::create($size, $prob)->logsfAt(...));
    }

    /** stats_nbinom_quantile */
    public static function quantile(float $p, float $size, float $prob): float
    {
        return Support::probability($p, self::create($size, $prob)->quantileAt(...));
    }

    /** stats_nbinom_quantile_array */
    public static function quantileArray(array $probabilities, float $size, float $prob): array
    {
        return Support::probabilities($probabilities, self::create($size, $prob)->quantileAt(...));
    }

    /** stats_nbinom_rng */
    public static function rng(int $n, int $seed, float $size, float $prob): array
    {
        return Support::sample($n, $seed, self::create($size, $prob)->quantileAt(...));
    }

    private function exact(float $k): ?array
    {
        if ($k < 0.0) {
            return [0.0, 1.0];
        }
        if ($this->q == 0.0 || \is_infinite($k)) {
            return [1.0, 0.0];
        }

        return null;
    }

    private function tails(float $k): array
    {
        if (($exact = $this->exact($k)) !== null) {
            return $exact;
        }
        $k = \floor($k);
        [$lower, $upper] = Support::betaInc($this->size, $k + 1.0, $this->p, $this->q);

        return [
            Special::refine($lower, fn (): float => $this->lnLower($k)),
            Special::refine($upper, fn (): float => $this->lnUpper($k)),
        ];
    }

    private function lnLower(float $k): float
    {
        return Special::lnBetaLowerFraction($this->size, $k + 1.0, $this->p, $this->lnP, $this->lnQ);
    }

    private function lnUpper(float $k): float
    {
        return Special::lnBetaLowerFraction($k + 1.0, $this->size, $this->q, $this->lnQ, $this->lnP);
    }

    /** @internal */
    public function pmfAt(float $k): float
    {
        return \exp($this->logpmfAt($k));
    }

    /** @internal */
    public function logpmfAt(float $k): float
    {
        if (! Quantile::isCount($k)) {
            return -INF;
        }
        if ($this->q == 0.0) {
            return $k == 0.0 ? 0.0 : -INF;
        }
        if ($k == 0.0) {
            return $this->size * $this->lnP;
        }

        // size/(size + k) times the binomial density of size successes in size + k trials.
        return \log($this->size / ($this->size + $k)) + Special::lnBinomialRawSplit($this->size, $k, $this->p, $this->q);
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
            return \log($exact[0]);
        }
        [$lower, $upper] = $this->tails($k);

        return Special::logTail($lower, $upper, fn (): float => $this->lnLower(\floor($k)));
    }

    /** @internal */
    public function logsfAt(float $k): float
    {
        if (($exact = $this->exact($k)) !== null) {
            return \log($exact[1]);
        }
        [$lower, $upper] = $this->tails($k);

        return Special::logTail($upper, $lower, fn (): float => $this->lnUpper(\floor($k)));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        if ($this->q == 0.0) {
            return 0.0;
        }
        $mean = \fdiv($this->size * $this->q, $this->p);
        $sd = \sqrt(\fdiv($mean, $this->p));
        $first = Quantile::cornishFisher($p, $mean, $sd, \fdiv(1.0 + $this->q, \sqrt($this->size * $this->q)));

        return Quantile::discrete($p, $first, INF, $this->cdfAt(...), $this->sfAt(...), $this->logcdfAt(...));
    }
}
