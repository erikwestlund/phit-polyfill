<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** Log-normal distribution: ln X ~ Normal(meanlog, sdlog). */
final class LogNormal
{
    private function __construct(private readonly float $meanlog, private readonly float $sdlog) {}

    /** @internal */
    public static function create(float $meanlog, float $sdlog): self
    {
        return new self(Support::finite('meanlog', $meanlog), Support::positive('sdlog', $sdlog));
    }

    /** stats_lognormal_pdf */
    public static function pdf(float $x, float $meanlog = 0.0, float $sdlog = 1.0): float
    {
        return Support::point($x, self::create($meanlog, $sdlog)->pdfAt(...));
    }

    /** stats_lognormal_pdf_array */
    public static function pdfArray(array $values, float $meanlog = 0.0, float $sdlog = 1.0): array
    {
        return Support::points($values, self::create($meanlog, $sdlog)->pdfAt(...));
    }

    /** stats_lognormal_logpdf */
    public static function logpdf(float $x, float $meanlog = 0.0, float $sdlog = 1.0): float
    {
        return Support::point($x, self::create($meanlog, $sdlog)->logpdfAt(...));
    }

    /** stats_lognormal_logpdf_array */
    public static function logpdfArray(array $values, float $meanlog = 0.0, float $sdlog = 1.0): array
    {
        return Support::points($values, self::create($meanlog, $sdlog)->logpdfAt(...));
    }

    /** stats_lognormal_cdf */
    public static function cdf(float $x, float $meanlog = 0.0, float $sdlog = 1.0): float
    {
        return Support::point($x, self::create($meanlog, $sdlog)->cdfAt(...));
    }

    /** stats_lognormal_cdf_array */
    public static function cdfArray(array $values, float $meanlog = 0.0, float $sdlog = 1.0): array
    {
        return Support::points($values, self::create($meanlog, $sdlog)->cdfAt(...));
    }

    /** stats_lognormal_logcdf */
    public static function logcdf(float $x, float $meanlog = 0.0, float $sdlog = 1.0): float
    {
        return Support::point($x, self::create($meanlog, $sdlog)->logcdfAt(...));
    }

    /** stats_lognormal_logcdf_array */
    public static function logcdfArray(array $values, float $meanlog = 0.0, float $sdlog = 1.0): array
    {
        return Support::points($values, self::create($meanlog, $sdlog)->logcdfAt(...));
    }

    /** stats_lognormal_sf */
    public static function sf(float $x, float $meanlog = 0.0, float $sdlog = 1.0): float
    {
        return Support::point($x, self::create($meanlog, $sdlog)->sfAt(...));
    }

    /** stats_lognormal_sf_array */
    public static function sfArray(array $values, float $meanlog = 0.0, float $sdlog = 1.0): array
    {
        return Support::points($values, self::create($meanlog, $sdlog)->sfAt(...));
    }

    /** stats_lognormal_logsf */
    public static function logsf(float $x, float $meanlog = 0.0, float $sdlog = 1.0): float
    {
        return Support::point($x, self::create($meanlog, $sdlog)->logsfAt(...));
    }

    /** stats_lognormal_logsf_array */
    public static function logsfArray(array $values, float $meanlog = 0.0, float $sdlog = 1.0): array
    {
        return Support::points($values, self::create($meanlog, $sdlog)->logsfAt(...));
    }

    /** stats_lognormal_quantile */
    public static function quantile(float $p, float $meanlog = 0.0, float $sdlog = 1.0): float
    {
        return Support::probability($p, self::create($meanlog, $sdlog)->quantileAt(...));
    }

    /** stats_lognormal_quantile_array */
    public static function quantileArray(array $probabilities, float $meanlog = 0.0, float $sdlog = 1.0): array
    {
        return Support::probabilities($probabilities, self::create($meanlog, $sdlog)->quantileAt(...));
    }

    private function z(float $x): float
    {
        return $x <= 0.0 ? -INF : (log($x) - $this->meanlog) / $this->sdlog;
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        if ($x <= 0.0 || is_infinite($x)) {
            return 0.0;
        }

        return exp($this->logpdfAt($x));
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        if ($x <= 0.0 || is_infinite($x)) {
            return -INF;
        }

        return Normal::stdLogpdf($this->z($x)) - log($this->sdlog) - log($x);
    }

    /** @internal */
    public function cdfAt(float $x): float
    {
        return Normal::stdTails($this->z($x))[0];
    }

    /** @internal */
    public function sfAt(float $x): float
    {
        return Normal::stdTails($this->z($x))[1];
    }

    /** @internal */
    public function logcdfAt(float $x): float
    {
        return Normal::stdLogCdf($this->z($x));
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        return Normal::stdLogCdf(-$this->z($x));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        return exp($this->meanlog + $this->sdlog * Normal::stdQuantile($p));
    }
}
