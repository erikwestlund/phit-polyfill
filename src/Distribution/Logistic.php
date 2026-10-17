<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** Logistic distribution with location and scale. */
final class Logistic
{
    private function __construct(private readonly float $location, private readonly float $scale) {}

    /** @internal */
    public static function create(float $location, float $scale): self
    {
        return new self(Support::finite('location', $location), Support::positive('scale', $scale));
    }

    /** stats_logistic_pdf */
    public static function pdf(float $x, float $location = 0.0, float $scale = 1.0): float
    {
        return Support::point($x, self::create($location, $scale)->pdfAt(...));
    }

    /** stats_logistic_pdf_array */
    public static function pdfArray(array $values, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::points($values, self::create($location, $scale)->pdfAt(...));
    }

    /** stats_logistic_logpdf */
    public static function logpdf(float $x, float $location = 0.0, float $scale = 1.0): float
    {
        return Support::point($x, self::create($location, $scale)->logpdfAt(...));
    }

    /** stats_logistic_logpdf_array */
    public static function logpdfArray(array $values, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::points($values, self::create($location, $scale)->logpdfAt(...));
    }

    /** stats_logistic_cdf */
    public static function cdf(float $x, float $location = 0.0, float $scale = 1.0): float
    {
        return Support::point($x, self::create($location, $scale)->cdfAt(...));
    }

    /** stats_logistic_cdf_array */
    public static function cdfArray(array $values, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::points($values, self::create($location, $scale)->cdfAt(...));
    }

    /** stats_logistic_logcdf */
    public static function logcdf(float $x, float $location = 0.0, float $scale = 1.0): float
    {
        return Support::point($x, self::create($location, $scale)->logcdfAt(...));
    }

    /** stats_logistic_logcdf_array */
    public static function logcdfArray(array $values, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::points($values, self::create($location, $scale)->logcdfAt(...));
    }

    /** stats_logistic_sf */
    public static function sf(float $x, float $location = 0.0, float $scale = 1.0): float
    {
        return Support::point($x, self::create($location, $scale)->sfAt(...));
    }

    /** stats_logistic_sf_array */
    public static function sfArray(array $values, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::points($values, self::create($location, $scale)->sfAt(...));
    }

    /** stats_logistic_logsf */
    public static function logsf(float $x, float $location = 0.0, float $scale = 1.0): float
    {
        return Support::point($x, self::create($location, $scale)->logsfAt(...));
    }

    /** stats_logistic_logsf_array */
    public static function logsfArray(array $values, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::points($values, self::create($location, $scale)->logsfAt(...));
    }

    /** stats_logistic_quantile */
    public static function quantile(float $p, float $location = 0.0, float $scale = 1.0): float
    {
        return Support::probability($p, self::create($location, $scale)->quantileAt(...));
    }

    /** stats_logistic_quantile_array */
    public static function quantileArray(array $probabilities, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::probabilities($probabilities, self::create($location, $scale)->quantileAt(...));
    }

    /** stats_logistic_rng */
    public static function rng(int $n, int $seed, float $location = 0.0, float $scale = 1.0): array
    {
        return Support::sample($n, $seed, self::create($location, $scale)->quantileAt(...));
    }

    private function z(float $x): float
    {
        return ($x - $this->location) / $this->scale;
    }

    private static function standardCdf(float $z): float
    {
        if ($z >= 0.0) {
            return 1.0 / (1.0 + \exp(-$z));
        }
        $e = \exp($z);

        return $e / (1.0 + $e);
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        $e = \exp(-\abs($this->z($x)));

        return $e / ($this->scale * (1.0 + $e) * (1.0 + $e));
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        $a = \abs($this->z($x));

        return -$a - 2.0 * \log1p(\exp(-$a)) - \log($this->scale);
    }

    /** @internal */
    public function cdfAt(float $x): float
    {
        return self::standardCdf($this->z($x));
    }

    /** @internal */
    public function sfAt(float $x): float
    {
        return self::standardCdf(-$this->z($x));
    }

    /** @internal */
    public function logcdfAt(float $x): float
    {
        return -Special::softplus(-$this->z($x));
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        return -Special::softplus($this->z($x));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        if ($p <= 0.5) {
            $logit = \log($p) - \log1p(-$p);
        } else {
            $q = 1.0 - $p;
            $logit = \log1p(-$q) - \log($q);
        }

        return $this->location + $this->scale * $logit;
    }
}
