<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * The crate's util.rs and consts.rs: PRIMA's evaluate/moderate/checkexit/xinbd
 * helpers, the x0 revision, interval_max, and the exit codes (infos.f90 values).
 *
 * @internal
 */
final class Util
{
    public const float FUNCMAX = 1.0e30;

    public const float BOUNDMAX = 0.25 * PHP_FLOAT_MAX;

    public const int INFO_DFT = 0;

    public const int SMALL_TR_RADIUS = 0;

    public const int FTARGET_ACHIEVED = 1;

    public const int MAXFUN_REACHED = 3;

    public const int MAXTR_REACHED = 20;

    public const int NAN_INF_X = -1;

    public const int NAN_INF_F = -2;

    public const int NAN_INF_MODEL = -3;

    public const int DAMAGING_ROUNDING = 7;

    public const int FTOL_REACHED = 1000;

    /** PRIMA's single-precision literals widened to f64, as gfortran evaluates them. */
    public const float F32_1E_6 = 9.999999974752427e-7;

    public const float F32_1E_20 = 9.999999682655225e-21;

    public const float F32_1E_2 = 0.009999999776482582;

    public static function moderatex1(float $v): float
    {
        if (\is_nan($v)) {
            $v = 0.0;
        }

        return Linalg::max(-Linalg::REALMAX, Linalg::min(Linalg::REALMAX, $v));
    }

    public static function moderatef(float $f): float
    {
        if (\is_nan($f)) {
            $f = self::FUNCMAX;
        }

        return Linalg::max(-Linalg::REALMAX, Linalg::min(self::FUNCMAX, $f));
    }

    /**
     * PRIMA's moderated-extreme-barrier evaluation.
     *
     * @param  callable(list<float>): (float|int)  $calfun
     * @param  list<float>  $x
     */
    public static function evaluate(callable $calfun, array $x): float
    {
        foreach ($x as $v) {
            if (\is_nan($v)) {
                return Linalg::sum($x);
            }
        }
        $xmod = [];
        foreach ($x as $v) {
            $xmod[] = self::moderatex1($v);
        }

        return self::moderatef((float) $calfun($xmod));
    }

    /**
     * Revise x0 so its distance to each inactive bound is 0 or at least rhobeg.
     *
     * @param  list<float>  $x
     * @param  list<float>  $xl
     * @param  list<float>  $xu
     */
    public static function reviseX0(array &$x, array $xl, array $xu, float $rhobeg): void
    {
        $n = count($x);
        for ($i = 0; $i < $n; $i++) {
            if ($x[$i] <= $xl[$i] + 0.5 * $rhobeg) {
                $x[$i] = $xl[$i];
            } elseif ($x[$i] < $xl[$i] + $rhobeg) {
                $x[$i] = $xl[$i] + $rhobeg;
            }
        }
        for ($i = 0; $i < $n; $i++) {
            if ($x[$i] >= $xu[$i] - 0.5 * $rhobeg) {
                $x[$i] = $xu[$i];
            } elseif ($x[$i] > $xu[$i] - $rhobeg) {
                $x[$i] = $xu[$i] - $rhobeg;
            }
        }
    }

    /** @param  list<float>  $x */
    public static function checkexit(int $maxfun, int $nf, float $f, float $ftarget, array $x): int
    {
        $info = self::INFO_DFT;
        foreach ($x as $v) {
            if (\is_nan($v) || \is_infinite($v)) {
                $info = self::NAN_INF_X;
                break;
            }
        }
        if (\is_nan($f) || (\is_infinite($f) && $f > 0.0)) {
            $info = self::NAN_INF_F;
        }
        if ($f <= $ftarget) {
            $info = self::FTARGET_ACHIEVED;
        }
        if ($nf >= $maxfun) {
            $info = self::MAXFUN_REACHED;
        }

        return $info;
    }

    /**
     * X = XBASE + STEP, projected so that bound-hitting steps land exactly on the bounds.
     *
     * @param  list<float>  $xbase
     * @param  list<float>  $step
     * @param  list<float>  $xl
     * @param  list<float>  $xu
     * @param  list<float>  $sl
     * @param  list<float>  $su
     * @return list<float>
     */
    public static function xinbd(array $xbase, array $step, array $xl, array $xu, array $sl, array $su): array
    {
        $x = [];
        $n = count($xbase);
        for ($i = 0; $i < $n; $i++) {
            $s = Linalg::max($sl[$i], Linalg::min($su[$i], $step[$i]));
            $xi = Linalg::max($xl[$i], Linalg::min($xu[$i], $xbase[$i] + $s));
            if ($s <= $sl[$i]) {
                $xi = $xl[$i];
            }
            if ($s >= $su[$i]) {
                $xi = $xu[$i];
            }
            $x[] = $xi;
        }

        return $x;
    }

    /**
     * PRIMA's interval_max specialised to trsbox's objective: approximate maximiser on
     * [lb, ub] from a grid search refined by one quadratic-interpolation step.
     *
     * @param  array{float, float, float, float, float}  $args
     */
    public static function intervalMax(float $lb, float $ub, array $args, int $gridSize): float
    {
        if ($ub <= $lb) {
            return $lb;
        }
        $xgrid = self::linspace($lb, $ub, $gridSize);
        $fgrid = [];
        $allNan = true;
        foreach ($xgrid as $xg) {
            $fg = self::intervalFunTrsbox($xg, $args);
            $fgrid[] = $fg;
            if (! \is_nan($fg)) {
                $allNan = false;
            }
        }
        if ($allNan) {
            return $lb;
        }

        $kopt = 0;
        $fopt = NAN;
        $found = false;
        foreach ($fgrid as $k => $f) {
            if (! \is_nan($f) && (! $found || $f > $fopt)) {
                $kopt = $k;
                $fopt = $f;
                $found = true;
            }
        }

        if ($kopt === 0) {
            return $lb;
        }
        if ($kopt === $gridSize - 1) {
            return $ub;
        }
        $fprev = $fgrid[$kopt - 1];
        $fnext = $fgrid[$kopt + 1];
        $step = 0.0;
        if (\abs($fprev - $fnext) > 0.0) {
            $step = 0.5 * (($fnext - $fprev) / ($fopt + $fopt - $fprev - $fnext));
        }
        if (\is_finite($step) && \abs($step) > 0.0) {
            return $lb + ($ub - $lb) * ((float) $kopt + $step) / (float) ($gridSize - 1);
        }

        return $xgrid[$kopt];
    }

    /**
     * PRIMA trustregion.f90 interval_fun_trsbox; args = [shs, dhd, dhs, dredg, sredg].
     *
     * @param  array{float, float, float, float, float}  $args
     */
    public static function intervalFunTrsbox(float $hangt, array $args): float
    {
        $f = 0.0;
        if (\abs($hangt) > 0.0) {
            $sth = ($hangt + $hangt) / (1.0 + $hangt * $hangt);
            $f = $args[0] + $hangt * ($hangt * $args[1] - $args[2] - $args[2]);
            $f = $sth * ($hangt * $args[3] - $args[4] - 0.5 * $sth * $f);
        }

        return $f;
    }

    /** @return list<float> */
    private static function linspace(float $xstart, float $xstop, int $n): array
    {
        if ($n === 0) {
            return [];
        }
        $nm = $n - 1;
        $x = [];
        if ($n === 1 || ($xstart <= $xstop && $xstop <= $xstart)) {
            $x = \array_fill(0, $n, $xstop);
        } elseif (\abs($xstart) <= \abs($xstop) && \abs($xstop) <= \abs($xstart)) {
            $xunit = $xstop / (float) $nm;
            for ($idx = 0; $idx < $n; $idx++) {
                $x[] = $xunit * (2.0 * (float) $idx - (float) $nm);
            }
            if ($nm % 2 === 0) {
                $x[\intdiv($nm, 2)] = 0.0;
            }
        } else {
            $xunit = ($xstop - $xstart) / (float) $nm;
            for ($i = 0; $i < $n; $i++) {
                $x[] = $xstart + $xunit * (float) $i;
            }
        }
        $x[0] = $xstart;
        $x[$n - 1] = $xstop;

        return $x;
    }
}
