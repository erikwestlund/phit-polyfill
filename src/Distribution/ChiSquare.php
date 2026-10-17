<?php

declare(strict_types=1);

namespace Phit\Polyfill\Distribution;

/** Chi-square distribution with df degrees of freedom: Gamma(df/2, scale 2). */
final class ChiSquare
{
    private function __construct(private readonly Gamma $gamma) {}

    /** @internal */
    public static function create(float $df): self
    {
        return new self(Gamma::create(0.5 * Support::positive('df', $df), 2.0));
    }

    /** Phit\dist_chisquare_pdf */
    public static function pdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->pdfAt(...));
    }

    /** Phit\dist_chisquare_pdf_array */
    public static function pdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->pdfAt(...));
    }

    /** Phit\dist_chisquare_logpdf */
    public static function logpdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->logpdfAt(...));
    }

    /** Phit\dist_chisquare_logpdf_array */
    public static function logpdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->logpdfAt(...));
    }

    /** Phit\dist_chisquare_cdf */
    public static function cdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->cdfAt(...));
    }

    /** Phit\dist_chisquare_cdf_array */
    public static function cdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->cdfAt(...));
    }

    /** Phit\dist_chisquare_logcdf */
    public static function logcdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->logcdfAt(...));
    }

    /** Phit\dist_chisquare_logcdf_array */
    public static function logcdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->logcdfAt(...));
    }

    /** Phit\dist_chisquare_sf */
    public static function sf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->sfAt(...));
    }

    /** Phit\dist_chisquare_sf_array */
    public static function sfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->sfAt(...));
    }

    /** Phit\dist_chisquare_logsf */
    public static function logsf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->logsfAt(...));
    }

    /** Phit\dist_chisquare_logsf_array */
    public static function logsfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->logsfAt(...));
    }

    /** Phit\dist_chisquare_quantile */
    public static function quantile(float $p, float $df): float
    {
        return Support::probability($p, self::create($df)->quantileAt(...));
    }

    /** Phit\dist_chisquare_quantile_array */
    public static function quantileArray(array $probabilities, float $df): array
    {
        return Support::probabilities($probabilities, self::create($df)->quantileAt(...));
    }

    /** Phit\dist_chisquare_rng */
    public static function rng(int $n, int $seed, float $df): array
    {
        return Support::sample($n, $seed, self::create($df)->quantileAt(...));
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        return $this->gamma->pdfAt($x);
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        return $this->gamma->logpdfAt($x);
    }

    /** @internal */
    public function cdfAt(float $x): float
    {
        return $this->gamma->cdfAt($x);
    }

    /** @internal */
    public function sfAt(float $x): float
    {
        return $this->gamma->sfAt($x);
    }

    /** @internal */
    public function logcdfAt(float $x): float
    {
        return $this->gamma->logcdfAt($x);
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        return $this->gamma->logsfAt($x);
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        return $this->gamma->quantileAt($p);
    }
}
