<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/** Exponential distribution with rate λ (mean 1/λ). */
final class Exponential
{
    private function __construct(private readonly float $rate) {}

    /** @internal */
    public static function create(float $rate): self
    {
        return new self(Support::positive('rate', $rate));
    }

    /** Phit\dist_exponential_pdf */
    public static function pdf(float $x, float $rate = 1.0): float
    {
        return Support::point($x, self::create($rate)->pdfAt(...));
    }

    /** Phit\dist_exponential_pdf_array */
    public static function pdfArray(array $values, float $rate = 1.0): array
    {
        return Support::points($values, self::create($rate)->pdfAt(...));
    }

    /** Phit\dist_exponential_logpdf */
    public static function logpdf(float $x, float $rate = 1.0): float
    {
        return Support::point($x, self::create($rate)->logpdfAt(...));
    }

    /** Phit\dist_exponential_logpdf_array */
    public static function logpdfArray(array $values, float $rate = 1.0): array
    {
        return Support::points($values, self::create($rate)->logpdfAt(...));
    }

    /** Phit\dist_exponential_cdf */
    public static function cdf(float $x, float $rate = 1.0): float
    {
        return Support::point($x, self::create($rate)->cdfAt(...));
    }

    /** Phit\dist_exponential_cdf_array */
    public static function cdfArray(array $values, float $rate = 1.0): array
    {
        return Support::points($values, self::create($rate)->cdfAt(...));
    }

    /** Phit\dist_exponential_logcdf */
    public static function logcdf(float $x, float $rate = 1.0): float
    {
        return Support::point($x, self::create($rate)->logcdfAt(...));
    }

    /** Phit\dist_exponential_logcdf_array */
    public static function logcdfArray(array $values, float $rate = 1.0): array
    {
        return Support::points($values, self::create($rate)->logcdfAt(...));
    }

    /** Phit\dist_exponential_sf */
    public static function sf(float $x, float $rate = 1.0): float
    {
        return Support::point($x, self::create($rate)->sfAt(...));
    }

    /** Phit\dist_exponential_sf_array */
    public static function sfArray(array $values, float $rate = 1.0): array
    {
        return Support::points($values, self::create($rate)->sfAt(...));
    }

    /** Phit\dist_exponential_logsf */
    public static function logsf(float $x, float $rate = 1.0): float
    {
        return Support::point($x, self::create($rate)->logsfAt(...));
    }

    /** Phit\dist_exponential_logsf_array */
    public static function logsfArray(array $values, float $rate = 1.0): array
    {
        return Support::points($values, self::create($rate)->logsfAt(...));
    }

    /** Phit\dist_exponential_quantile */
    public static function quantile(float $p, float $rate = 1.0): float
    {
        return Support::probability($p, self::create($rate)->quantileAt(...));
    }

    /** Phit\dist_exponential_quantile_array */
    public static function quantileArray(array $probabilities, float $rate = 1.0): array
    {
        return Support::probabilities($probabilities, self::create($rate)->quantileAt(...));
    }

    /** Phit\dist_exponential_rng */
    public static function rng(int $n, int $seed, float $rate = 1.0): array
    {
        return Support::sample($n, $seed, self::create($rate)->quantileAt(...));
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        return $x < 0.0 ? 0.0 : $this->rate * exp(-$this->rate * $x);
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        return $x < 0.0 ? -INF : log($this->rate) - $this->rate * $x;
    }

    /** @internal */
    public function cdfAt(float $x): float
    {
        return $x <= 0.0 ? 0.0 : -expm1(-$this->rate * $x);
    }

    /** @internal */
    public function sfAt(float $x): float
    {
        return $x <= 0.0 ? 1.0 : exp(-$this->rate * $x);
    }

    /** @internal */
    public function logcdfAt(float $x): float
    {
        return $x <= 0.0 ? -INF : Special::lnOneMinusExpNeg($this->rate * $x);
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        return $x <= 0.0 ? 0.0 : -$this->rate * $x;
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        return $p <= 0.5 ? -log1p(-$p) / $this->rate : -log(1.0 - $p) / $this->rate;
    }
}
