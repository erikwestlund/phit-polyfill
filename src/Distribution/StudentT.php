<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** Student's t distribution with df degrees of freedom (location 0, scale 1). */
final class StudentT
{
    private function __construct(private readonly float $df) {}

    /** @internal */
    public static function create(float $df): self
    {
        return new self(Support::positive('df', $df));
    }

    /** stats_t_pdf */
    public static function pdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->pdfAt(...));
    }

    /** stats_t_pdf_array */
    public static function pdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->pdfAt(...));
    }

    /** stats_t_logpdf */
    public static function logpdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->logpdfAt(...));
    }

    /** stats_t_logpdf_array */
    public static function logpdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->logpdfAt(...));
    }

    /** stats_t_cdf */
    public static function cdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->cdfAt(...));
    }

    /** stats_t_cdf_array */
    public static function cdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->cdfAt(...));
    }

    /** stats_t_logcdf */
    public static function logcdf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->logcdfAt(...));
    }

    /** stats_t_logcdf_array */
    public static function logcdfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->logcdfAt(...));
    }

    /** stats_t_sf */
    public static function sf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->sfAt(...));
    }

    /** stats_t_sf_array */
    public static function sfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->sfAt(...));
    }

    /** stats_t_logsf */
    public static function logsf(float $x, float $df): float
    {
        return Support::point($x, self::create($df)->logsfAt(...));
    }

    /** stats_t_logsf_array */
    public static function logsfArray(array $values, float $df): array
    {
        return Support::points($values, self::create($df)->logsfAt(...));
    }

    /** stats_t_quantile */
    public static function quantile(float $p, float $df): float
    {
        return Support::probability($p, self::create($df)->quantileAt(...));
    }

    /** stats_t_quantile_array */
    public static function quantileArray(array $probabilities, float $df): array
    {
        return Support::probabilities($probabilities, self::create($df)->quantileAt(...));
    }

    /** stats_t_rng */
    public static function rng(int $n, int $seed, float $df): array
    {
        return Support::sample($n, $seed, self::create($df)->quantileAt(...));
    }

    /**
     * Beta-function arguments for the tail P(T > |x|) = I_h(df/2, 1/2) / 2.
     *
     * @return array{0: float, 1: float, 2: float, 3: float} [h, y, ln h, ln y]
     */
    private function args(float $x): array
    {
        $u = \abs($x) / \sqrt($this->df);
        $r = $u * $u;
        if (\is_finite($r)) {
            return [1.0 / (1.0 + $r), $r / (1.0 + $r), -\log1p($r), -\log1p(\fdiv(1.0, $r))];
        }
        // |x|/√df beyond 1e154: work with ln(|x|/√df) to avoid overflow.
        $lnU = \log(\abs($x)) - 0.5 * \log($this->df);
        $h = \exp(-2.0 * $lnU);

        return [$h, 1.0, -2.0 * $lnU, -$h];
    }

    private function lnSmallTail(array $a): float
    {
        return -M_LN2 + Special::lnBetaLowerFraction(0.5 * $this->df, 0.5, $a[0], $a[2], $a[3]);
    }

    /** (P(T ≤ x), P(T > x)) with the small tail computed directly. */
    private function tails(float $x): array
    {
        $a = $this->args($x);
        [$i, $oneMinusI] = $a[0] == 0.0 ? [0.0, 1.0] : Support::betaInc(0.5 * $this->df, 0.5, $a[0], $a[1]);
        $small = Special::refine(0.5 * $i, fn (): float => $this->lnSmallTail($a));
        $big = $oneMinusI + 0.5 * $i;

        return $x <= 0.0 ? [$small, $big] : [$big, $small];
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        return \exp($this->logpdfAt($x));
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        $half = 0.5 * $this->df;
        // ln Γ((df + 1)/2) - ln Γ(df/2) without cancellation for large df.
        $ratio = $half >= 8.0
            ? -Cdflib::algdiv(0.5, $half)
            : Cdflib::gammaLog($half + 0.5) - Cdflib::gammaLog($half);
        $a = $this->args($x);

        return $ratio - 0.5 * \log(M_PI * $this->df) + ($half + 0.5) * $a[2];
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
        [$lower, $upper] = $this->tails($x);

        return Special::logTail($lower, $upper, fn (): float => $this->lnSmallTail($this->args($x)));
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        [$lower, $upper] = $this->tails($x);

        return Special::logTail($upper, $lower, fn (): float => $this->lnSmallTail($this->args($x)));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        if ($p == 0.5) {
            return 0.0;
        }
        if ($p > 0.5) {
            return -$this->quantileAt(1.0 - $p);
        }
        $start = $p > 0.0 ? Cdflib::dt1($p, 1.0 - $p, $this->df) : 0.0;

        return Quantile::continuous($p, $start, -INF, INF, $this->logcdfAt(...), $this->logsfAt(...), $this->logpdfAt(...));
    }
}
