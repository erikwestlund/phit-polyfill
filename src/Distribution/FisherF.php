<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Distribution;

/** F distribution with numerator df1 and denominator df2 degrees of freedom. */
final class FisherF
{
    private function __construct(private readonly float $df1, private readonly float $df2) {}

    /** @internal */
    public static function create(float $df1, float $df2): self
    {
        return new self(Support::positive('df1', $df1), Support::positive('df2', $df2));
    }

    /** stats_f_pdf */
    public static function pdf(float $x, float $df1, float $df2): float
    {
        return Support::point($x, self::create($df1, $df2)->pdfAt(...));
    }

    /** stats_f_pdf_array */
    public static function pdfArray(array $values, float $df1, float $df2): array
    {
        return Support::points($values, self::create($df1, $df2)->pdfAt(...));
    }

    /** stats_f_logpdf */
    public static function logpdf(float $x, float $df1, float $df2): float
    {
        return Support::point($x, self::create($df1, $df2)->logpdfAt(...));
    }

    /** stats_f_logpdf_array */
    public static function logpdfArray(array $values, float $df1, float $df2): array
    {
        return Support::points($values, self::create($df1, $df2)->logpdfAt(...));
    }

    /** stats_f_cdf */
    public static function cdf(float $x, float $df1, float $df2): float
    {
        return Support::point($x, self::create($df1, $df2)->cdfAt(...));
    }

    /** stats_f_cdf_array */
    public static function cdfArray(array $values, float $df1, float $df2): array
    {
        return Support::points($values, self::create($df1, $df2)->cdfAt(...));
    }

    /** stats_f_logcdf */
    public static function logcdf(float $x, float $df1, float $df2): float
    {
        return Support::point($x, self::create($df1, $df2)->logcdfAt(...));
    }

    /** stats_f_logcdf_array */
    public static function logcdfArray(array $values, float $df1, float $df2): array
    {
        return Support::points($values, self::create($df1, $df2)->logcdfAt(...));
    }

    /** stats_f_sf */
    public static function sf(float $x, float $df1, float $df2): float
    {
        return Support::point($x, self::create($df1, $df2)->sfAt(...));
    }

    /** stats_f_sf_array */
    public static function sfArray(array $values, float $df1, float $df2): array
    {
        return Support::points($values, self::create($df1, $df2)->sfAt(...));
    }

    /** stats_f_logsf */
    public static function logsf(float $x, float $df1, float $df2): float
    {
        return Support::point($x, self::create($df1, $df2)->logsfAt(...));
    }

    /** stats_f_logsf_array */
    public static function logsfArray(array $values, float $df1, float $df2): array
    {
        return Support::points($values, self::create($df1, $df2)->logsfAt(...));
    }

    /** stats_f_quantile */
    public static function quantile(float $p, float $df1, float $df2): float
    {
        return Support::probability($p, self::create($df1, $df2)->quantileAt(...));
    }

    /** stats_f_quantile_array */
    public static function quantileArray(array $probabilities, float $df1, float $df2): array
    {
        return Support::probabilities($probabilities, self::create($df1, $df2)->quantileAt(...));
    }

    /** stats_f_rng */
    public static function rng(int $n, int $seed, float $df1, float $df2): array
    {
        return Support::sample($n, $seed, self::create($df1, $df2)->quantileAt(...));
    }

    /**
     * Beta arguments u = df1 x/(df1 x + df2) and v = 1 - u, with their logs.
     *
     * @return array{0: float, 1: float, 2: float, 3: float} [u, v, ln u, ln v]
     */
    private function args(float $x): array
    {
        $r = $this->df1 / $this->df2 * $x;
        if (is_finite($r)) {
            return [$r / (1.0 + $r), 1.0 / (1.0 + $r), -log1p(fdiv(1.0, $r)), -log1p($r)];
        }
        $v = fdiv($this->df2 / $this->df1, $x);

        return [1.0, $v, -$v, log($this->df2 / $this->df1) - log($x)];
    }

    private function lnLower(array $a): float
    {
        return Special::lnBetaLowerFraction(0.5 * $this->df1, 0.5 * $this->df2, $a[0], $a[2], $a[3]);
    }

    private function lnUpper(array $a): float
    {
        return Special::lnBetaLowerFraction(0.5 * $this->df2, 0.5 * $this->df1, $a[1], $a[3], $a[2]);
    }

    private static function exact(float $x): ?array
    {
        if ($x <= 0.0) {
            return [0.0, 1.0];
        }
        if (is_infinite($x)) {
            return [1.0, 0.0];
        }

        return null;
    }

    private function tails(float $x): array
    {
        if (($exact = self::exact($x)) !== null) {
            return $exact;
        }
        $a = $this->args($x);
        [$lower, $upper] = $a[1] == 0.0
            ? [1.0, 0.0]
            : Support::betaInc(0.5 * $this->df1, 0.5 * $this->df2, $a[0], $a[1]);

        return [
            Special::refine($lower, fn (): float => $this->lnLower($a)),
            Special::refine($upper, fn (): float => $this->lnUpper($a)),
        ];
    }

    /** @internal */
    public function pdfAt(float $x): float
    {
        return exp($this->logpdfAt($x));
    }

    /** @internal */
    public function logpdfAt(float $x): float
    {
        $m = $this->df1;
        $n = $this->df2;
        if ($x < 0.0 || is_infinite($x)) {
            return -INF;
        }
        if ($x == 0.0) {
            return $m < 2.0 ? INF : ($m == 2.0 ? 0.0 : -INF);
        }
        $a = $this->args($x);
        if ($a[1] < 1e-100) {
            // Far right tail: the direct formula has no cancellation here.
            return 0.5 * $m * log($m / $n) + (0.5 * $m - 1.0) * log($x) + 0.5 * ($m + $n) * $a[3]
                - Cdflib::betaLog(0.5 * $m, 0.5 * $n);
        }
        if ($m >= 2.0) {
            return log(0.5 * $m * $a[1])
                + Special::lnBinomialRaw(0.5 * $m - 1.0, 0.5 * ($m + $n) - 1.0, $a[0], $a[1]);
        }

        return log(fdiv($m * $m * $a[1], 2.0 * $a[0] * ($m + $n)))
            + Special::lnBinomialRaw(0.5 * $m, 0.5 * ($m + $n), $a[0], $a[1]);
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

        return Special::logTail($lower, $upper, fn (): float => $this->lnLower($this->args($x)));
    }

    /** @internal */
    public function logsfAt(float $x): float
    {
        if (($exact = self::exact($x)) !== null) {
            return log($exact[1]);
        }
        [$lower, $upper] = $this->tails($x);

        return Special::logTail($upper, $lower, fn (): float => $this->lnUpper($this->args($x)));
    }

    /** @internal */
    public function quantileAt(float $p): float
    {
        return Quantile::continuous($p, 1.0, 0.0, INF, $this->logcdfAt(...), $this->logsfAt(...), $this->logpdfAt(...));
    }
}
