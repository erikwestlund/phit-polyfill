<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** Poisson distribution with mean `lambda`. */
final class Poisson
{
    private function __construct(private readonly float $lambda) {}

    /** @internal */
    public static function create(float $lambda): self
    {
        Support::finite('lambda', $lambda);
        if ($lambda < 0.0) {
            throw Support::invalid('lambda must be nonnegative');
        }

        return new self($lambda);
    }

    /** stats_poisson_pmf */
    public static function pmf(float $k, float $lambda): float
    {
        return Support::point($k, self::create($lambda)->pmfAt(...));
    }

    /** stats_poisson_pmf_array */
    public static function pmfArray(array $counts, float $lambda): array
    {
        return Support::points($counts, self::create($lambda)->pmfAt(...));
    }

    /** stats_poisson_logpmf */
    public static function logpmf(float $k, float $lambda): float
    {
        return Support::point($k, self::create($lambda)->logpmfAt(...));
    }

    /** stats_poisson_logpmf_array */
    public static function logpmfArray(array $counts, float $lambda): array
    {
        return Support::points($counts, self::create($lambda)->logpmfAt(...));
    }

    /** stats_poisson_cdf */
    public static function cdf(float $k, float $lambda): float
    {
        return Support::point($k, self::create($lambda)->cdfAt(...));
    }

    /** stats_poisson_cdf_array */
    public static function cdfArray(array $counts, float $lambda): array
    {
        return Support::points($counts, self::create($lambda)->cdfAt(...));
    }

    /** stats_poisson_logcdf */
    public static function logcdf(float $k, float $lambda): float
    {
        return Support::point($k, self::create($lambda)->logcdfAt(...));
    }

    /** stats_poisson_logcdf_array */
    public static function logcdfArray(array $counts, float $lambda): array
    {
        return Support::points($counts, self::create($lambda)->logcdfAt(...));
    }

    /** stats_poisson_sf */
    public static function sf(float $k, float $lambda): float
    {
        return Support::point($k, self::create($lambda)->sfAt(...));
    }

    /** stats_poisson_sf_array */
    public static function sfArray(array $counts, float $lambda): array
    {
        return Support::points($counts, self::create($lambda)->sfAt(...));
    }

    /** stats_poisson_logsf */
    public static function logsf(float $k, float $lambda): float
    {
        return Support::point($k, self::create($lambda)->logsfAt(...));
    }

    /** stats_poisson_logsf_array */
    public static function logsfArray(array $counts, float $lambda): array
    {
        return Support::points($counts, self::create($lambda)->logsfAt(...));
    }

    /** stats_poisson_quantile */
    public static function quantile(float $p, float $lambda): float
    {
        return Support::probability($p, self::create($lambda)->quantileAt(...));
    }

    /** stats_poisson_quantile_array */
    public static function quantileArray(array $probabilities, float $lambda): array
    {
        return Support::probabilities($probabilities, self::create($lambda)->quantileAt(...));
    }

    /** stats_poisson_rng */
    public static function rng(int $n, int $seed, float $lambda): array
    {
        return Support::sample($n, $seed, self::create($lambda)->quantileAt(...));
    }

    private function exact(float $k): ?array
    {
        if ($k < 0.0) {
            return [0.0, 1.0];
        }
        if ($this->lambda == 0.0 || \is_infinite($k)) {
            return [1.0, 0.0];
        }

        return null;
    }

    private function tails(float $k): array
    {
        if (($exact = $this->exact($k)) !== null) {
            return $exact;
        }
        $a = \floor($k) + 1.0;
        [$lowerGamma, $upperGamma] = Support::gammaInc($a, $this->lambda);

        return [
            Special::refine($upperGamma, fn (): float => Special::lnGammaUpperFraction($a, $this->lambda)),
            Special::refine($lowerGamma, fn (): float => Special::lnGammaLowerSeries($a, $this->lambda)),
        ];
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

        return Special::lnPoissonRaw($k, $this->lambda);
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

        return Special::logTail($lower, $upper, fn (): float => Special::lnGammaUpperFraction(\floor($k) + 1.0, $this->lambda));
    }

    /** @internal */
    public function logsfAt(float $k): float
    {
        if (($exact = $this->exact($k)) !== null) {
            return \log($exact[1]);
        }
        [$lower, $upper] = $this->tails($k);

        return Special::logTail($upper, $lower, fn (): float => Special::lnGammaLowerSeries(\floor($k) + 1.0, $this->lambda));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        if ($this->lambda == 0.0) {
            return 0.0;
        }
        $sd = \sqrt($this->lambda);
        $first = Quantile::cornishFisher($p, $this->lambda, $sd, 1.0 / $sd);

        return Quantile::discrete($p, $first, INF, $this->cdfAt(...), $this->sfAt(...), $this->logcdfAt(...));
    }
}
