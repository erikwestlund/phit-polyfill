<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal\Bobyqa;

/**
 * The crate's trustregion.rs: PRIMA's trsbox (truncated CG plus a 2-D boundary
 * search under box bounds) and trrad.
 *
 * @internal
 */
final class TrustRegion
{
    /**
     * Approximately minimize Q(XOPT + D) subject to ||D|| <= DELTA and SL <= XOPT + D <= SU.
     * Writes the step into $d and returns CRVMIN.
     *
     * @param  list<float>  $goptIn
     * @param  list<float>  $hqIn
     * @param  list<float>  $pqIn
     * @param  list<float>  $sl
     * @param  list<float>  $su
     * @param  list<float>  $xopt
     * @param  list<float>  $xpt
     * @param  list<float>  $d
     */
    public static function trsbox(
        float $delta,
        array $goptIn,
        array $hqIn,
        array $pqIn,
        array $sl,
        array $su,
        float $tol,
        array $xopt,
        array $xpt,
        int $n,
        int $npt,
        array &$d,
    ): float {
        $maxAbsGopt = 0.0;
        foreach ($goptIn as $v) {
            $maxAbsGopt = Linalg::max($maxAbsGopt, abs($v));
        }
        if ($maxAbsGopt > 1.0e12) {
            $modscal = Linalg::max(2.0 * Linalg::REALMIN, 1.0 / $maxAbsGopt);
            $gopt = [];
            foreach ($goptIn as $v) {
                $gopt[] = $v * $modscal;
            }
            $pq = [];
            foreach ($pqIn as $v) {
                $pq[] = $v * $modscal;
            }
            $hq = [];
            foreach ($hqIn as $v) {
                $hq[] = $v * $modscal;
            }
            $scaled = true;
        } else {
            $gopt = $goptIn;
            $pq = $pqIn;
            $hq = $hqIn;
            $scaled = false;
            $modscal = 1.0;
        }

        $dredsq = 0.0;
        $ggsav = 0.0;

        $xbdi = array_fill(0, $n, 0);
        for ($i = 0; $i < $n; $i++) {
            if ($xopt[$i] >= $su[$i] && $gopt[$i] <= 0.0) {
                $xbdi[$i] = 1;
            }
        }
        for ($i = 0; $i < $n; $i++) {
            if ($xopt[$i] <= $sl[$i] && $gopt[$i] >= 0.0) {
                $xbdi[$i] = -1;
            }
        }
        $nact = self::countActive($xbdi);

        $d = array_fill(0, $n, 0.0);
        $crvmin = -Linalg::REALMAX;

        $gnew = $gopt;
        $gredsq = self::maskedSumsq($gnew, $xbdi);
        $delsq = $delta * $delta;
        $qred = 0.0;
        $beta = 0.0;
        $itercg = 0;
        $twodSearch = false;
        $s = array_fill(0, $n, 0.0);
        $hs = [];

        $maxiter = min(10000, ($n - $nact) * ($n - $nact));

        for ($iter = 0; $iter < $maxiter; $iter++) {
            $resid = $delsq - self::maskedSumsq($d, $xbdi);
            if ($resid <= 0.0) {
                $twodSearch = true;
                break;
            }

            if ($itercg === 0) {
                for ($i = 0; $i < $n; $i++) {
                    $s[$i] = -$gnew[$i];
                }
            } else {
                for ($i = 0; $i < $n; $i++) {
                    $s[$i] = $beta * $s[$i] - $gnew[$i];
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if ($xbdi[$i] !== 0) {
                    $s[$i] = 0.0;
                }
            }
            $stepsq = Linalg::inprod($s, $s);
            $ds = self::maskedInprod($d, $s, $xbdi);

            // PRIMA's .not.(...) test, spelled out: NaN in any operand also exits.
            $thresholdA = Linalg::EPS * $delsq;
            $lhsB = $gredsq * $delsq;
            $rhsB = ($tol * $qred) * ($tol * $qred);
            if (is_nan($ds) || is_nan($stepsq) || is_nan($thresholdA) || is_nan($lhsB) || is_nan($rhsB)
                || $stepsq <= $thresholdA || $lhsB <= $rhsB) {
                break;
            }

            $sqrtd = Linalg::max(Linalg::max(sqrt($stepsq * $resid + $ds * $ds), sqrt($stepsq * $resid)), abs($ds));

            if ($ds >= 0.0) {
                $bstep = $resid / ($sqrtd + $ds);
            } else {
                $bstep = ($sqrtd - $ds) / $stepsq;
            }
            if ($bstep <= 0.0 || ! is_finite($bstep)) {
                break;
            }

            $hs = PowAlg::hessMul($s, $xpt, $n, $npt, $pq, $hq);
            $shs = self::maskedInprod($s, $hs, $xbdi);
            $stplen = $bstep;
            if ($shs > 0.0) {
                $stplen = Linalg::min($bstep, $gredsq / $shs);
            }

            $xnew = [];
            for ($i = 0; $i < $n; $i++) {
                $xnew[] = $xopt[$i] + $d[$i];
            }
            $xtest = [];
            for ($i = 0; $i < $n; $i++) {
                $xtest[] = $xnew[$i] + $stplen * $s[$i];
            }
            $sbound = array_fill(0, $n, $stplen);
            for ($i = 0; $i < $n; $i++) {
                if ($s[$i] > 0.0 && $xtest[$i] > $su[$i]) {
                    $sbound[$i] = ($su[$i] - $xnew[$i]) / $s[$i];
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if ($s[$i] < 0.0 && $xtest[$i] < $sl[$i]) {
                    $sbound[$i] = ($sl[$i] - $xnew[$i]) / $s[$i];
                }
            }
            $anyBelow = false;
            for ($i = 0; $i < $n; $i++) {
                if (is_nan($sbound[$i])) {
                    $sbound[$i] = $stplen;
                }
                if ($sbound[$i] < $stplen) {
                    $anyBelow = true;
                }
            }
            $iact = null;
            if ($anyBelow) {
                $imin = 0;
                for ($i = 1; $i < $n; $i++) {
                    if ($sbound[$i] < $sbound[$imin]) {
                        $imin = $i;
                    }
                }
                $iact = $imin;
                $stplen = $sbound[$imin];
            }

            $sdec = 0.0;
            if ($stplen > 0.0) {
                $itercg++;
                $rayleighq = $shs / $stepsq;
                if ($iact === null && $rayleighq > 0.0) {
                    if ($crvmin <= -Linalg::REALMAX) {
                        $crvmin = $rayleighq;
                    } else {
                        $crvmin = Linalg::min($crvmin, $rayleighq);
                    }
                }
                $ggsav = $gredsq;
                for ($i = 0; $i < $n; $i++) {
                    $gnew[$i] += $stplen * $hs[$i];
                }
                $gredsq = self::maskedSumsq($gnew, $xbdi);
                $dold = $d;
                for ($i = 0; $i < $n; $i++) {
                    $d[$i] += $stplen * $s[$i];
                }
                if (! is_finite(Linalg::sumAbs($d))) {
                    $d = $dold;
                    break;
                }
                $sdec = Linalg::max($stplen * ($ggsav - 0.5 * $stplen * $shs), 0.0);
                $qred += $sdec;
            } else {
                $ggsav = $gredsq;
            }

            if ($iact !== null) {
                $nact++;
                $xbdi[$iact] = (int) Linalg::copysign(1.0, $s[$iact]);
                if ($nact >= $n) {
                    break;
                }
                $delsq -= $d[$iact] * $d[$iact];
                if ($delsq <= 0.0) {
                    $twodSearch = true;
                    break;
                }
                $beta = 0.0;
                $itercg = 0;
                $gredsq = self::maskedSumsq($gnew, $xbdi);
            } elseif ($stplen < $bstep) {
                if ($itercg >= $n - $nact || $sdec <= $tol * $qred || is_nan($sdec) || is_nan($qred)) {
                    break;
                }
                $beta = $gredsq / $ggsav;
            } else {
                $twodSearch = true;
                break;
            }
        }

        if ($twodSearch) {
            $crvmin = 0.0;
            $maxiter = 10 * ($n - $nact);
        } else {
            $maxiter = 0;
        }

        $nactsav = $nact - 1;
        $hdred = [];
        for ($iter1 = 1; $iter1 <= $maxiter; $iter1++) {
            $xnew = [];
            for ($i = 0; $i < $n; $i++) {
                $xnew[] = $xopt[$i] + $d[$i];
            }

            for ($i = 0; $i < $n; $i++) {
                if ($xbdi[$i] === 0 && $xnew[$i] >= $su[$i]) {
                    $xbdi[$i] = 1;
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if ($xbdi[$i] === 0 && $xnew[$i] <= $sl[$i]) {
                    $xbdi[$i] = -1;
                }
            }
            $nact = self::countActive($xbdi);
            if ($nact >= $n - 1) {
                break;
            }

            $gredsq = self::maskedSumsq($gnew, $xbdi);
            $dredg = self::maskedInprod($d, $gnew, $xbdi);
            if ($iter1 === 1 || $nact > $nactsav) {
                $dredsq = self::maskedSumsq($d, $xbdi);
                $dred = $d;
                for ($i = 0; $i < $n; $i++) {
                    if ($xbdi[$i] !== 0) {
                        $dred[$i] = 0.0;
                    }
                }
                $hdred = PowAlg::hessMul($dred, $xpt, $n, $npt, $pq, $hq);
                $nactsav = $nact;
            }

            $temp = $gredsq * $dredsq - $dredg * $dredg;
            $bound = $tol * $tol * Linalg::max($gredsq * $dredsq, $qred * $qred);
            if (is_nan($temp) || is_nan($bound) || $temp <= $bound) {
                break;
            }
            $temp = sqrt($temp);
            for ($i = 0; $i < $n; $i++) {
                $s[$i] = ($dredg * $d[$i] - $dredsq * $gnew[$i]) / $temp;
            }
            for ($i = 0; $i < $n; $i++) {
                if ($xbdi[$i] !== 0) {
                    $s[$i] = 0.0;
                }
            }
            $sredg = -$temp;

            $ssq = [];
            for ($i = 0; $i < $n; $i++) {
                $ssq[] = $d[$i] * $d[$i] + $s[$i] * $s[$i];
            }
            $tanbd = array_fill(0, $n, 1.0);
            $sqdscr = array_fill(0, $n, -Linalg::REALMAX);
            for ($i = 0; $i < $n; $i++) {
                if ($xbdi[$i] === 0 && $xopt[$i] - $sl[$i] < sqrt($ssq[$i])) {
                    $sqdscr[$i] = sqrt(Linalg::max(0.0, $ssq[$i] - ($xopt[$i] - $sl[$i]) * ($xopt[$i] - $sl[$i])));
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if ($sqdscr[$i] - $s[$i] > 0.0) {
                    $tanbd[$i] = Linalg::min($tanbd[$i], ($xnew[$i] - $sl[$i]) / ($sqdscr[$i] - $s[$i]));
                }
            }
            $sqdscr = array_fill(0, $n, -Linalg::REALMAX);
            for ($i = 0; $i < $n; $i++) {
                if ($xbdi[$i] === 0 && $su[$i] - $xopt[$i] < sqrt($ssq[$i])) {
                    $sqdscr[$i] = sqrt(Linalg::max(0.0, $ssq[$i] - ($su[$i] - $xopt[$i]) * ($su[$i] - $xopt[$i])));
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if ($sqdscr[$i] + $s[$i] > 0.0) {
                    $tanbd[$i] = Linalg::min($tanbd[$i], ($su[$i] - $xnew[$i]) / ($sqdscr[$i] + $s[$i]));
                }
            }
            $anyBelow = false;
            for ($i = 0; $i < $n; $i++) {
                if (is_nan($tanbd[$i])) {
                    $tanbd[$i] = 0.0;
                }
                if ($tanbd[$i] < 1.0) {
                    $anyBelow = true;
                }
            }

            $iact = null;
            $hangtBd = 1.0;
            if ($anyBelow) {
                $imin = 0;
                for ($i = 1; $i < $n; $i++) {
                    if ($tanbd[$i] < $tanbd[$imin]) {
                        $imin = $i;
                    }
                }
                $iact = $imin;
                $hangtBd = $tanbd[$imin];
            }
            if ($hangtBd <= 0.0) {
                break;
            }

            $hs = PowAlg::hessMul($s, $xpt, $n, $npt, $pq, $hq);
            $shs = self::maskedInprod($s, $hs, $xbdi);
            $dhs = self::maskedInprod($d, $hs, $xbdi);
            $dhd = self::maskedInprod($d, $hdred, $xbdi);

            $args = [$shs, $dhd, $dhs, $dredg, $sredg];
            if (is_nan($shs) || is_nan($dhd) || is_nan($dhs) || is_nan($dredg) || is_nan($sredg)) {
                break;
            }
            $gridSize = 2 * (int) self::roundHalfAway(17.0 * $hangtBd + 4.1);
            $hangt = Util::intervalMax(0.0, $hangtBd, $args, $gridSize);
            $sdec = Util::intervalFunTrsbox($hangt, $args);
            if (is_nan($sdec) || $sdec <= 0.0) {
                break;
            }

            $cth = Linalg::min((1.0 - $hangt * $hangt) / (1.0 + $hangt * $hangt), 1.0 - $hangt * $hangt);
            $sth = Linalg::min(($hangt + $hangt) / (1.0 + $hangt * $hangt), $hangt + $hangt);
            for ($i = 0; $i < $n; $i++) {
                $gnew[$i] = $gnew[$i] + ($cth - 1.0) * $hdred[$i] + $sth * $hs[$i];
            }
            $dold = $d;
            for ($i = 0; $i < $n; $i++) {
                if ($xbdi[$i] === 0) {
                    $d[$i] = $cth * $d[$i] + $sth * $s[$i];
                }
            }
            if (! is_finite(Linalg::sumAbs($d))) {
                $d = $dold;
                break;
            }
            for ($i = 0; $i < $n; $i++) {
                $hdred[$i] = $cth * $hdred[$i] + $sth * $hs[$i];
            }
            $qred += $sdec;
            $small = is_nan($sdec) || is_nan($qred) || $sdec <= $tol * $qred;
            if ($iact !== null) {
                if ($hangt >= $hangtBd) {
                    $xbdi[$iact] = (int) Linalg::copysign(1.0, $xopt[$iact] + $d[$iact] - Linalg::midpoint($sl[$iact], $su[$iact]));
                } elseif ($small) {
                    break;
                }
            } elseif ($small) {
                break;
            }
        }

        for ($i = 0; $i < $n; $i++) {
            $xi = Linalg::max(Linalg::min($xopt[$i] + $d[$i], $su[$i]), $sl[$i]);
            if ($xbdi[$i] === -1) {
                $xi = $sl[$i];
            }
            if ($xbdi[$i] === 1) {
                $xi = $su[$i];
            }
            $d[$i] = $xi - $xopt[$i];
        }

        if ($crvmin <= -Linalg::REALMAX || is_nan($crvmin)) {
            $crvmin = 0.0;
        }
        if ($scaled && $crvmin > 0.0) {
            $crvmin /= $modscal;
        }

        return $crvmin;
    }

    /** Update the trust-region radius according to RATIO and DNORM. */
    public static function trrad(float $delta, float $dnorm, float $eta1, float $eta2, float $gamma1, float $gamma2, float $ratio): float
    {
        if ($ratio <= $eta1) {
            return Linalg::min($gamma1 * $delta, $dnorm);
        }
        if ($ratio <= $eta2) {
            return Linalg::max($gamma1 * $delta, $dnorm);
        }

        return Linalg::max($gamma1 * $delta, $gamma2 * $dnorm);
    }

    /** f64::round (half away from zero) for a non-negative argument. */
    private static function roundHalfAway(float $x): float
    {
        $r = floor($x);

        return $x - $r >= 0.5 ? $r + 1.0 : $r;
    }

    /** @param  list<int>  $xbdi */
    private static function countActive(array $xbdi): int
    {
        $count = 0;
        foreach ($xbdi as $v) {
            if ($v !== 0) {
                $count++;
            }
        }

        return $count;
    }

    /**
     * @param  list<float>  $x
     * @param  list<int>  $xbdi
     */
    private static function maskedSumsq(array $x, array $xbdi): float
    {
        $acc = 0.0;
        foreach ($x as $i => $v) {
            if ($xbdi[$i] === 0) {
                $acc += $v * $v;
            }
        }

        return $acc;
    }

    /**
     * @param  list<float>  $x
     * @param  list<float>  $y
     * @param  list<int>  $xbdi
     */
    private static function maskedInprod(array $x, array $y, array $xbdi): float
    {
        $acc = 0.0;
        foreach ($x as $i => $v) {
            if ($xbdi[$i] === 0) {
                $acc += $v * $y[$i];
            }
        }

        return $acc;
    }
}
