<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/** Negative binomial distribution with size and mean mu (Var = mu + mu²/size). */
final class NegativeBinomialMu
{
    private function __construct() {}

    /** @internal */
    public static function create(float $size, float $mu): NegativeBinomial
    {
        return NegativeBinomial::withMean($size, $mu);
    }

    /** Phit\dist_nbinom_mu_pmf */
    public static function pmf(float $k, float $size, float $mu): float
    {
        return Support::point($k, self::create($size, $mu)->pmfAt(...));
    }

    /** Phit\dist_nbinom_mu_pmf_array */
    public static function pmfArray(array $counts, float $size, float $mu): array
    {
        return Support::points($counts, self::create($size, $mu)->pmfAt(...));
    }

    /** Phit\dist_nbinom_mu_logpmf */
    public static function logpmf(float $k, float $size, float $mu): float
    {
        return Support::point($k, self::create($size, $mu)->logpmfAt(...));
    }

    /** Phit\dist_nbinom_mu_logpmf_array */
    public static function logpmfArray(array $counts, float $size, float $mu): array
    {
        return Support::points($counts, self::create($size, $mu)->logpmfAt(...));
    }

    /** Phit\dist_nbinom_mu_cdf */
    public static function cdf(float $k, float $size, float $mu): float
    {
        return Support::point($k, self::create($size, $mu)->cdfAt(...));
    }

    /** Phit\dist_nbinom_mu_cdf_array */
    public static function cdfArray(array $counts, float $size, float $mu): array
    {
        return Support::points($counts, self::create($size, $mu)->cdfAt(...));
    }

    /** Phit\dist_nbinom_mu_logcdf */
    public static function logcdf(float $k, float $size, float $mu): float
    {
        return Support::point($k, self::create($size, $mu)->logcdfAt(...));
    }

    /** Phit\dist_nbinom_mu_logcdf_array */
    public static function logcdfArray(array $counts, float $size, float $mu): array
    {
        return Support::points($counts, self::create($size, $mu)->logcdfAt(...));
    }

    /** Phit\dist_nbinom_mu_sf */
    public static function sf(float $k, float $size, float $mu): float
    {
        return Support::point($k, self::create($size, $mu)->sfAt(...));
    }

    /** Phit\dist_nbinom_mu_sf_array */
    public static function sfArray(array $counts, float $size, float $mu): array
    {
        return Support::points($counts, self::create($size, $mu)->sfAt(...));
    }

    /** Phit\dist_nbinom_mu_logsf */
    public static function logsf(float $k, float $size, float $mu): float
    {
        return Support::point($k, self::create($size, $mu)->logsfAt(...));
    }

    /** Phit\dist_nbinom_mu_logsf_array */
    public static function logsfArray(array $counts, float $size, float $mu): array
    {
        return Support::points($counts, self::create($size, $mu)->logsfAt(...));
    }

    /** Phit\dist_nbinom_mu_quantile */
    public static function quantile(float $p, float $size, float $mu): float
    {
        return Support::probability($p, self::create($size, $mu)->quantileAt(...));
    }

    /** Phit\dist_nbinom_mu_quantile_array */
    public static function quantileArray(array $probabilities, float $size, float $mu): array
    {
        return Support::probabilities($probabilities, self::create($size, $mu)->quantileAt(...));
    }

    /** Phit\dist_nbinom_mu_rng */
    public static function rng(int $n, int $seed, float $size, float $mu): array
    {
        return Support::sample($n, $seed, self::create($size, $mu)->quantileAt(...));
    }
}
