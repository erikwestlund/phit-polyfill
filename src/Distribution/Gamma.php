<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** Gamma distribution with shape k and scale θ (mean kθ). */
final class Gamma
{
    private function __construct(private readonly float $shape, private readonly float $scale) {}

    /** @internal */
    public static function create(float $shape, float $scale): self
    {
        return new self(Support::positive('shape', $shape), Support::positive('scale', $scale));
    }

    /** stats_gamma_pdf */
    public static function pdf(float $x, float $shape, float $scale = 1.0): float
    {
        return Support::point($x, self::create($shape, $scale)->pdfAt(...));
    }

    /** stats_gamma_pdf_array */
    public static function pdfArray(array $values, float $shape, float $scale = 1.0): array
    {
        return Support::points($values, self::create($shape, $scale)->pdfAt(...));
    }

    /** stats_gamma_logpdf */
    public static function logpdf(float $x, float $shape, float $scale = 1.0): float
    {
        return Support::point($x, self::create($shape, $scale)->logpdfAt(...));
    }

    /** stats_gamma_logpdf_array */
    public static function logpdfArray(array $values, float $shape, float $scale = 1.0): array
    {
        return Support::points($values, self::create($shape, $scale)->logpdfAt(...));
    }

    /** stats_gamma_cdf */
    public static function cdf(float $x, float $shape, float $scale = 1.0): float
    {
        return Support::point($x, self::create($shape, $scale)->cdfAt(...));
    }

    /** stats_gamma_cdf_array */
    public static function cdfArray(array $values, float $shape, float $scale = 1.0): array
    {
        return Support::points($values, self::create($shape, $scale)->cdfAt(...));
    }

    /** stats_gamma_logcdf */
    public static function logcdf(float $x, float $shape, float $scale = 1.0): float
    {
        return Support::point($x, self::create($shape, $scale)->logcdfAt(...));
    }

    /** stats_gamma_logcdf_array */
    public static function logcdfArray(array $values, float $shape, float $scale = 1.0): array
    {
        return Support::points($values, self::create($shape, $scale)->logcdfAt(...));
    }

    /** stats_gamma_sf */
    public static function sf(float $x, float $shape, float $scale = 1.0): float
    {
        return Support::point($x, self::create($shape, $scale)->sfAt(...));
    }

    /** stats_gamma_sf_array */
    public static function sfArray(array $values, float $shape, float $scale = 1.0): array
    {
        return Support::points($values, self::create($shape, $scale)->sfAt(...));
    }

    /** stats_gamma_logsf */
    public static function logsf(float $x, float $shape, float $scale = 1.0): float
    {
        return Support::point($x, self::create($shape, $scale)->logsfAt(...));
    }

    /** stats_gamma_logsf_array */
    public static function logsfArray(array $values, float $shape, float $scale = 1.0): array
    {
        return Support::points($values, self::create($shape, $scale)->logsfAt(...));
    }

    /** stats_gamma_quantile */
    public static function quantile(float $p, float $shape, float $scale = 1.0): float
    {
        return Support::probability($p, self::create($shape, $scale)->quantileAt(...));
    }

    /** stats_gamma_quantile_array */
    public static function quantileArray(array $probabilities, float $shape, float $scale = 1.0): array
    {
        return Support::probabilities($probabilities, self::create($shape, $scale)->quantileAt(...));
    }

    // Standard (scale 1) pieces in y = x / scale.

    private function stdLogpdf(float $y): float
    {
        $a = $this->shape;
        if ($y < 0.0 || is_infinite($y)) {
            return -INF;
        }
        if ($y == 0.0) {
            return $a < 1.0 ? INF : ($a == 1.0 ? 0.0 : -INF);
        }
        if ($y < 1e-100) {
            return ($a - 1.0) * log($y) - $y - Cdflib::gammaLog($a);
        }
        if ($a < 1.0) {
            return Special::lnPoissonRaw($a, $y) + log($a) - log($y);
        }

        return Special::lnPoissonRaw($a - 1.0, $y);
    }

    private static function stdExact(float $y): ?array
    {
        if ($y <= 0.0) {
            return [0.0, 1.0];
        }
        if (is_infinite($y)) {
            return [1.0, 0.0];
        }

        return null;
    }

    private function stdTails(float $y): array
    {
        if (($exact = self::stdExact($y)) !== null) {
            return $exact;
        }
        [$p, $q] = Support::gammaInc($this->shape, $y);

        return [
            Special::refine($p, fn (): float => Special::lnGammaLowerSeries($this->shape, $y)),
            Special::refine($q, fn (): float => Special::lnGammaUpperFraction($this->shape, $y)),
        ];
    }

    private function stdLogcdf(float $y): float
    {
        if (($exact = self::stdExact($y)) !== null) {
            return log($exact[0]);
        }
        [$p, $q] = $this->stdTails($y);

        return Special::logTail($p, $q, fn (): float => Special::lnGammaLowerSeries($this->shape, $y));
    }

    private function stdLogsf(float $y): float
    {
        if (($exact = self::stdExact($y)) !== null) {
            return log($exact[1]);
        }
        [$p, $q] = $this->stdTails($y);

        return Special::logTail($q, $p, fn (): float => Special::lnGammaUpperFraction($this->shape, $y));
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        return exp($this->logpdfAt($x));
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        return $this->stdLogpdf($x / $this->scale) - log($this->scale);
    }

    /** @internal */
    public function cdfAt(float $x): float
    {
        return $this->stdTails($x / $this->scale)[0];
    }

    /** @internal */
    public function sfAt(float $x): float
    {
        return $this->stdTails($x / $this->scale)[1];
    }

    /** @internal */
    public function logcdfAt(float $x): float
    {
        return $this->stdLogcdf($x / $this->scale);
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        return $this->stdLogsf($x / $this->scale);
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        $a = $this->shape;
        // Wilson–Hilferty start, or the small-x power law when it is not positive.
        $c = 1.0 / (9.0 * $a);
        $w = 1.0 - $c + Normal::stdQuantile($p) * sqrt($c);
        $start = $w > 0.0
            ? $a * $w * $w * $w
            : exp((log($p) + Cdflib::gammaLog($a + 1.0)) / $a);

        return $this->scale * Quantile::continuous(
            $p,
            $start,
            0.0,
            INF,
            $this->stdLogcdf(...),
            $this->stdLogsf(...),
            $this->stdLogpdf(...),
        );
    }
}
