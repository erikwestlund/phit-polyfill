<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** Beta distribution with shape parameters α and β. */
final class Beta
{
    private function __construct(private readonly float $alpha, private readonly float $beta) {}

    /** @internal */
    public static function create(float $alpha, float $beta): self
    {
        return new self(Support::positive('alpha', $alpha), Support::positive('beta', $beta));
    }

    /** stats_beta_pdf */
    public static function pdf(float $x, float $alpha, float $beta): float
    {
        return Support::point($x, self::create($alpha, $beta)->pdfAt(...));
    }

    /** stats_beta_pdf_array */
    public static function pdfArray(array $values, float $alpha, float $beta): array
    {
        return Support::points($values, self::create($alpha, $beta)->pdfAt(...));
    }

    /** stats_beta_logpdf */
    public static function logpdf(float $x, float $alpha, float $beta): float
    {
        return Support::point($x, self::create($alpha, $beta)->logpdfAt(...));
    }

    /** stats_beta_logpdf_array */
    public static function logpdfArray(array $values, float $alpha, float $beta): array
    {
        return Support::points($values, self::create($alpha, $beta)->logpdfAt(...));
    }

    /** stats_beta_cdf */
    public static function cdf(float $x, float $alpha, float $beta): float
    {
        return Support::point($x, self::create($alpha, $beta)->cdfAt(...));
    }

    /** stats_beta_cdf_array */
    public static function cdfArray(array $values, float $alpha, float $beta): array
    {
        return Support::points($values, self::create($alpha, $beta)->cdfAt(...));
    }

    /** stats_beta_logcdf */
    public static function logcdf(float $x, float $alpha, float $beta): float
    {
        return Support::point($x, self::create($alpha, $beta)->logcdfAt(...));
    }

    /** stats_beta_logcdf_array */
    public static function logcdfArray(array $values, float $alpha, float $beta): array
    {
        return Support::points($values, self::create($alpha, $beta)->logcdfAt(...));
    }

    /** stats_beta_sf */
    public static function sf(float $x, float $alpha, float $beta): float
    {
        return Support::point($x, self::create($alpha, $beta)->sfAt(...));
    }

    /** stats_beta_sf_array */
    public static function sfArray(array $values, float $alpha, float $beta): array
    {
        return Support::points($values, self::create($alpha, $beta)->sfAt(...));
    }

    /** stats_beta_logsf */
    public static function logsf(float $x, float $alpha, float $beta): float
    {
        return Support::point($x, self::create($alpha, $beta)->logsfAt(...));
    }

    /** stats_beta_logsf_array */
    public static function logsfArray(array $values, float $alpha, float $beta): array
    {
        return Support::points($values, self::create($alpha, $beta)->logsfAt(...));
    }

    /** stats_beta_quantile */
    public static function quantile(float $p, float $alpha, float $beta): float
    {
        return Support::probability($p, self::create($alpha, $beta)->quantileAt(...));
    }

    /** stats_beta_quantile_array */
    public static function quantileArray(array $probabilities, float $alpha, float $beta): array
    {
        return Support::probabilities($probabilities, self::create($alpha, $beta)->quantileAt(...));
    }

    private function lnLower(float $x): float
    {
        return Special::lnBetaLowerFraction($this->alpha, $this->beta, $x, log($x), log1p(-$x));
    }

    private function lnUpper(float $x): float
    {
        return Special::lnBetaLowerFraction($this->beta, $this->alpha, 1.0 - $x, log1p(-$x), log($x));
    }

    private static function exact(float $x): ?array
    {
        if ($x <= 0.0) {
            return [0.0, 1.0];
        }
        if ($x >= 1.0) {
            return [1.0, 0.0];
        }

        return null;
    }

    private function tails(float $x): array
    {
        if (($exact = self::exact($x)) !== null) {
            return $exact;
        }
        [$lower, $upper] = Support::betaInc($this->alpha, $this->beta, $x, 1.0 - $x);

        return [
            Special::refine($lower, fn (): float => $this->lnLower($x)),
            Special::refine($upper, fn (): float => $this->lnUpper($x)),
        ];
    }

    private static function edge(float $shape, float $other): float
    {
        if ($shape < 1.0) {
            return INF;
        }

        return $shape == 1.0 ? log($other) : -INF;
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        return exp($this->logpdfAt($x));
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        $a = $this->alpha;
        $b = $this->beta;
        if (! ($x >= 0.0 && $x <= 1.0)) {
            return -INF;
        }
        if ($x == 0.0) {
            return self::edge($a, $b);
        }
        if ($x == 1.0) {
            return self::edge($b, $a);
        }
        if ($a <= 2.0 || $b <= 2.0) {
            return ($a - 1.0) * log($x) + ($b - 1.0) * log1p(-$x) - Cdflib::betaLog($a, $b);
        }

        return log($a + $b - 1.0) + Special::lnBinomialRaw($a - 1.0, $a + $b - 2.0, $x, 1.0 - $x);
    }

    /** @internal */
    public function cdfAt(float $x): float
    {
        return $this->tails($x)[0];
    }

    /** @internal */
    public function sfAt(float $x): float
    {
        return $this->tails($x)[1];
    }

    /** @internal */
    public function logcdfAt(float $x): float
    {
        if (($exact = self::exact($x)) !== null) {
            return log($exact[0]);
        }
        [$lower, $upper] = $this->tails($x);

        return Special::logTail($lower, $upper, fn (): float => $this->lnLower($x));
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        if (($exact = self::exact($x)) !== null) {
            return log($exact[1]);
        }
        [$lower, $upper] = $this->tails($x);

        return Special::logTail($upper, $lower, fn (): float => $this->lnUpper($x));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        $a = $this->alpha;
        $b = $this->beta;
        $mean = $a / ($a + $b);
        if ($p < 0.01) {
            $start = exp((log($p) + log($a) + Cdflib::betaLog($a, $b)) / $a);
            $start = self::min($start, $mean);
        } elseif ($p > 0.99) {
            $tail = exp((log(1.0 - $p) + log($b) + Cdflib::betaLog($a, $b)) / $b);
            $start = 1.0 - self::min($tail, 1.0 - $mean);
        } else {
            $sd = sqrt(fdiv($a * $b, ($a + $b) * ($a + $b) * ($a + $b + 1.0)));
            $start = $mean + $sd * Normal::stdQuantile($p);
            // f64::clamp(min, max)
            $lo = 1e-3 * $mean;
            $hi = 1.0 - 1e-3 * (1.0 - $mean);
            if ($start < $lo) {
                $start = $lo;
            }
            if ($start > $hi) {
                $start = $hi;
            }
        }

        return Quantile::continuous($p, $start, 0.0, 1.0, $this->logcdfAt(...), $this->logsfAt(...), $this->logpdfAt(...));
    }

    /** f64::min: the non-NaN operand when one is NaN. */
    private static function min(float $a, float $b): float
    {
        if (is_nan($a)) {
            return $b;
        }
        if (is_nan($b)) {
            return $a;
        }

        return $a < $b ? $a : $b;
    }
}
