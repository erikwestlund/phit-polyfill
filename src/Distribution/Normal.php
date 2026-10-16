<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** Normal distribution with mean and standard deviation. */
final class Normal
{
    private const float INV_SQRT_2PI = 0.398_942_280_401_432_7;

    private function __construct(private readonly float $mean, private readonly float $sd) {}

    /** @internal */
    public static function create(float $mean, float $sd): self
    {
        return new self(Support::finite('mean', $mean), Support::positive('sd', $sd));
    }

    /** stats_normal_pdf */
    public static function pdf(float $x, float $mean = 0.0, float $sd = 1.0): float
    {
        return Support::point($x, self::create($mean, $sd)->pdfAt(...));
    }

    /** stats_normal_pdf_array */
    public static function pdfArray(array $values, float $mean = 0.0, float $sd = 1.0): array
    {
        return Support::points($values, self::create($mean, $sd)->pdfAt(...));
    }

    /** stats_normal_logpdf */
    public static function logpdf(float $x, float $mean = 0.0, float $sd = 1.0): float
    {
        return Support::point($x, self::create($mean, $sd)->logpdfAt(...));
    }

    /** stats_normal_logpdf_array */
    public static function logpdfArray(array $values, float $mean = 0.0, float $sd = 1.0): array
    {
        return Support::points($values, self::create($mean, $sd)->logpdfAt(...));
    }

    /** stats_normal_cdf */
    public static function cdf(float $x, float $mean = 0.0, float $sd = 1.0): float
    {
        return Support::point($x, self::create($mean, $sd)->cdfAt(...));
    }

    /** stats_normal_cdf_array */
    public static function cdfArray(array $values, float $mean = 0.0, float $sd = 1.0): array
    {
        return Support::points($values, self::create($mean, $sd)->cdfAt(...));
    }

    /** stats_normal_logcdf */
    public static function logcdf(float $x, float $mean = 0.0, float $sd = 1.0): float
    {
        return Support::point($x, self::create($mean, $sd)->logcdfAt(...));
    }

    /** stats_normal_logcdf_array */
    public static function logcdfArray(array $values, float $mean = 0.0, float $sd = 1.0): array
    {
        return Support::points($values, self::create($mean, $sd)->logcdfAt(...));
    }

    /** stats_normal_sf */
    public static function sf(float $x, float $mean = 0.0, float $sd = 1.0): float
    {
        return Support::point($x, self::create($mean, $sd)->sfAt(...));
    }

    /** stats_normal_sf_array */
    public static function sfArray(array $values, float $mean = 0.0, float $sd = 1.0): array
    {
        return Support::points($values, self::create($mean, $sd)->sfAt(...));
    }

    /** stats_normal_logsf */
    public static function logsf(float $x, float $mean = 0.0, float $sd = 1.0): float
    {
        return Support::point($x, self::create($mean, $sd)->logsfAt(...));
    }

    /** stats_normal_logsf_array */
    public static function logsfArray(array $values, float $mean = 0.0, float $sd = 1.0): array
    {
        return Support::points($values, self::create($mean, $sd)->logsfAt(...));
    }

    /** stats_normal_quantile */
    public static function quantile(float $p, float $mean = 0.0, float $sd = 1.0): float
    {
        return Support::probability($p, self::create($mean, $sd)->quantileAt(...));
    }

    /** stats_normal_quantile_array */
    public static function quantileArray(array $probabilities, float $mean = 0.0, float $sd = 1.0): array
    {
        return Support::probabilities($probabilities, self::create($mean, $sd)->quantileAt(...));
    }

    // ---------------------------------------------------------------- standard normal pieces (shared)

    /** @internal */
    public static function stdPdf(float $z): float
    {
        if (is_infinite($z)) {
            return 0.0;
        }
        if (abs($z) <= 5.0) {
            return self::INV_SQRT_2PI * exp(-0.5 * $z * $z);
        }
        // Split z so that z² is formed exactly in its leading part (as R's dnorm).
        $z1 = Support::round($z * 65536.0) / 65536.0;
        $z2 = $z - $z1;

        return self::INV_SQRT_2PI * exp(-0.5 * $z1 * $z1) * exp((-0.5 * $z2 - $z1) * $z2);
    }

    /** @internal */
    public static function stdLogpdf(float $z): float
    {
        return -0.5 * $z * $z - Special::LN_SQRT_2PI;
    }

    /**
     * @internal
     *
     * @return array{0: float, 1: float}
     */
    public static function stdTails(float $z): array
    {
        if ($z == INF) {
            return [1.0, 0.0];
        }
        if ($z == -INF) {
            return [0.0, 1.0];
        }

        return Cdflib::cumnor($z);
    }

    /** @internal */
    public static function stdLogCdf(float $z): float
    {
        [$lower, $upper] = self::stdTails($z);

        return Special::logTail($lower, $upper, static fn (): float => Special::normalLogLowerTail($z));
    }

    /** Standard normal quantile, polished on the log scale from CDFLIB's dinvnr. @internal */
    public static function stdQuantile(float $p): float
    {
        if ($p == 0.5) {
            return 0.0;
        }
        if ($p > 0.5) {
            return -self::stdQuantile(1.0 - $p);
        }
        $start = $p > 0.0 ? Cdflib::dinvnr($p, 1.0 - $p) : 0.0;

        return Quantile::continuous(
            $p,
            $start,
            -INF,
            INF,
            self::stdLogCdf(...),
            static fn (float $z): float => self::stdLogCdf(-$z),
            self::stdLogpdf(...),
        );
    }

    // ---------------------------------------------------------------- evaluation

    private function z(float $x): float
    {
        return ($x - $this->mean) / $this->sd;
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        if ($this->sd == 1.0) {
            return self::stdPdf($this->z($x));
        }

        return exp($this->logpdfAt($x));
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        return self::stdLogpdf($this->z($x)) - log($this->sd);
    }

    /** @internal */
    public function cdfAt(float $x): float
    {
        return self::stdTails($this->z($x))[0];
    }

    /** @internal */
    public function sfAt(float $x): float
    {
        return self::stdTails($this->z($x))[1];
    }

    /** @internal */
    public function logcdfAt(float $x): float
    {
        return self::stdLogCdf($this->z($x));
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        return self::stdLogCdf(-$this->z($x));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        return $this->mean + $this->sd * self::stdQuantile($p);
    }
}
