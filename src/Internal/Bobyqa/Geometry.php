<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * The crate's geometry.rs: PRIMA's setdrop_tr and geostep.
 *
 * @internal
 */
final class Geometry
{
    /**
     * The interpolation point to drop after a trust-region step, or null to drop none.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @param  list<float>|null  $denPrecomputed
     */
    public static function setdropTr(
        int $kopt,
        bool $ximproved,
        array $bmat,
        array $d,
        float $rho,
        array $xpt,
        array $zmat,
        int $n,
        int $npt,
        ?array $denPrecomputed,
    ): ?int {
        $koff = $kopt * $n;
        $distsq = [];
        for ($k = 0; $k < $npt; $k++) {
            $offset = $k * $n;
            $sq = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $diff = $ximproved
                    ? $xpt[$offset + $i] - ($xpt[$koff + $i] + $d[$i])
                    : $xpt[$offset + $i] - $xpt[$koff + $i];
                $sq += $diff * $diff;
            }
            $distsq[] = $sq;
        }

        $rho2 = $rho * $rho;
        $weight = [];
        foreach ($distsq as $v) {
            $m = Linalg::max($v / $rho2, 1.0);
            $m2 = $m * $m;
            $weight[] = $m2 * $m2;
        }

        $den = $denPrecomputed ?? PowAlg::calden($kopt, $bmat, $d, $xpt, $zmat, $n, $npt);
        $score = [];
        for ($k = 0; $k < $npt; $k++) {
            $score[] = $weight[$k] * $den[$k];
        }
        if (! $ximproved) {
            $score[$kopt] = -1.0;
        }
        $anyAboveOne = false;
        $anyPositive = false;
        for ($k = 0; $k < $npt; $k++) {
            if (\is_nan($score[$k])) {
                $score[$k] = -1.0;
            }
            if ($score[$k] > 1.0) {
                $anyAboveOne = true;
            }
            if ($score[$k] > 0.0) {
                $anyPositive = true;
            }
        }

        $knew = null;
        if ($anyAboveOne || ($ximproved && $anyPositive)) {
            $best = $score[0];
            $knew = 0;
            for ($k = 1; $k < $npt; $k++) {
                if ($score[$k] > $best) {
                    $best = $score[$k];
                    $knew = $k;
                }
            }
        }

        if ($ximproved && $knew === null) {
            $best = $distsq[0];
            $knew = 0;
            for ($k = 1; $k < $npt; $k++) {
                if ($distsq[$k] > $best) {
                    $best = $distsq[$k];
                    $knew = $k;
                }
            }
        }

        return $knew;
    }

    /**
     * A geometry-improving step D from XOPT.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $sl
     * @param  list<float>  $su
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @return list<float>
     */
    public static function geostep(
        int $knew,
        int $kopt,
        array $bmat,
        float $delbar,
        array $sl,
        array $su,
        array $xpt,
        array $zmat,
        int $n,
        int $npt,
    ): array {
        $nz = $npt - $n - 1;
        $zrow = [];
        for ($j = 0; $j < $nz; $j++) {
            $zrow[] = $zmat[$knew + $j * $npt];
        }
        $pqlag = Linalg::matprod21($zmat, $npt, $nz, $zrow);
        $alpha = $pqlag[$knew];

        $xopt = Linalg::col($xpt, $n, $kopt);

        $hm = PowAlg::hessMul($xopt, $xpt, $n, $npt, $pqlag, null);
        $glag = [];
        for ($i = 0; $i < $n; $i++) {
            $glag[] = $bmat[$i + $knew * $n] + $hm[$i];
        }

        if (! \is_finite(Linalg::sumAbs($glag))) {
            return self::fallbackStep($xpt, $xopt, $knew, $delbar, $n);
        }

        $dderiv = Linalg::matprod12($glag, $xpt, $n, $npt);
        $glagDotXopt = Linalg::inprod($glag, $xopt);
        for ($k = 0; $k < $npt; $k++) {
            $dderiv[$k] -= $glagDotXopt;
        }

        $distsq = [];
        for ($k = 0; $k < $npt; $k++) {
            $offset = $k * $n;
            $sq = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $diff = $xpt[$offset + $i] - $xopt[$i];
                $sq += $diff * $diff;
            }
            $distsq[] = $sq;
        }

        // stplen and isbd are 3 x npt; isbd keeps PRIMA's signed 1-based bound indices.
        $stplen = \array_fill(0, 3 * $npt, 0.0);
        $isbd = \array_fill(0, 3 * $npt, 0);
        $xdiff = \array_fill(0, $n, 0.0);
        $lfrac = \array_fill(0, $n, 0.0);
        $ufrac = \array_fill(0, $n, 0.0);

        for ($k = 0; $k < $npt; $k++) {
            if ($k === $kopt || \is_nan($dderiv[$k])) {
                $dderiv[$k] = 0.0;

                continue;
            }

            $subdInit = $delbar / \sqrt($distsq[$k]);
            $subd = $subdInit;
            $slbd = -$subdInit;
            $ilbd = 0;
            $iubd = 0;
            $sumin = Linalg::min(1.0, $subdInit);

            $offset = $k * $n;
            for ($i = 0; $i < $n; $i++) {
                $xdiff[$i] = $xpt[$offset + $i] - $xopt[$i];
                $lfrac[$i] = Linalg::copysign($subdInit, -$xdiff[$i]);
                if ($sl[$i] - $xopt[$i] > -\abs($xdiff[$i]) * $subdInit) {
                    $lfrac[$i] = ($sl[$i] - $xopt[$i]) / $xdiff[$i];
                }
                $ufrac[$i] = Linalg::copysign($subdInit, $xdiff[$i]);
                if ($su[$i] - $xopt[$i] < \abs($xdiff[$i]) * $subdInit) {
                    $ufrac[$i] = ($su[$i] - $xopt[$i]) / $xdiff[$i];
                }
            }

            $slbdTest = \array_fill(0, $n, $slbd);
            $anyBetter = false;
            for ($i = 0; $i < $n; $i++) {
                if ($xdiff[$i] > 0.0) {
                    $slbdTest[$i] = $lfrac[$i];
                } elseif ($xdiff[$i] < 0.0) {
                    $slbdTest[$i] = $ufrac[$i];
                }
                if ($slbdTest[$i] > $slbd) {
                    $anyBetter = true;
                }
            }
            if ($anyBetter) {
                $bestVal = -INF;
                $best = 0;
                for ($i = 0; $i < $n; $i++) {
                    if (! \is_nan($slbdTest[$i]) && $slbdTest[$i] > $bestVal) {
                        $bestVal = $slbdTest[$i];
                        $best = $i;
                    }
                }
                $slbd = $slbdTest[$best];
                $ilbd = -($best + 1) * (int) Linalg::copysign(1.0, $xdiff[$best]);
            }

            $subdTest = \array_fill(0, $n, $subdInit);
            $anyBetter = false;
            for ($i = 0; $i < $n; $i++) {
                if ($xdiff[$i] > 0.0) {
                    $subdTest[$i] = $ufrac[$i];
                } elseif ($xdiff[$i] < 0.0) {
                    $subdTest[$i] = $lfrac[$i];
                }
                if ($subdTest[$i] < $subdInit) {
                    $anyBetter = true;
                }
            }
            if ($anyBetter) {
                $bestVal = INF;
                $best = 0;
                for ($i = 0; $i < $n; $i++) {
                    if (! \is_nan($subdTest[$i]) && $subdTest[$i] < $bestVal) {
                        $bestVal = $subdTest[$i];
                        $best = $i;
                    }
                }
                $subd = Linalg::max($sumin, $subdTest[$best]);
                $iubd = ($best + 1) * (int) Linalg::copysign(1.0, $xdiff[$best]);
            }

            $stpm = 0.5;
            if ($k === $knew) {
                $stpm = $slbd;
                if (\abs(1.0 - $dderiv[$k]) > 0.0) {
                    $stpm = -0.5 * $dderiv[$k] / (1.0 - $dderiv[$k]);
                }
            }
            $stpm = Linalg::max($slbd, Linalg::min($subd, $stpm));

            $stplen[3 * $k] = $slbd;
            $stplen[3 * $k + 1] = $subd;
            $stplen[3 * $k + 2] = $stpm;
            $isbd[3 * $k] = $ilbd;
            $isbd[3 * $k + 1] = $iubd;
        }

        $vlag = [];
        for ($k = 0; $k < $npt; $k++) {
            $dd = $dderiv[$k];
            for ($i = 0; $i < 3; $i++) {
                $t = $stplen[3 * $k + $i];
                $vlag[] = $t * (1.0 - $t) * $dd;
            }
        }
        for ($i = 0; $i < 3; $i++) {
            $t = $stplen[3 * $knew + $i];
            $vlag[3 * $knew + $i] = $t * ($t * (1.0 - $dderiv[$knew]) + $dderiv[$knew]);
        }
        foreach ($vlag as $index => $v) {
            if (\is_nan($v)) {
                $vlag[$index] = 0.0;
            }
        }

        $predsq = [];
        for ($k = 0; $k < $npt; $k++) {
            $dk = $distsq[$k];
            for ($i = 0; $i < 3; $i++) {
                $t = $stplen[3 * $k + $i];
                $u = $t * (1.0 - $t) * $dk;
                $betabd = 0.5 * $u * $u;
                $v = $vlag[3 * $k + $i];
                $v2 = $v * $v;
                $p = $v2 * ($v2 + $alpha * $betabd);
                $predsq[] = \is_nan($p) ? 0.0 : $p;
            }
        }

        $ksqs = [0, 0, 0];
        for ($i = 0; $i < 3; $i++) {
            $best = $predsq[$i];
            $bestK = 0;
            for ($k = 1; $k < $npt; $k++) {
                if ($predsq[3 * $k + $i] > $best) {
                    $best = $predsq[3 * $k + $i];
                    $bestK = $k;
                }
            }
            $ksqs[$i] = $bestK;
        }
        $isq = 0;
        $bestPredsq = $predsq[3 * $ksqs[0]];
        for ($i = 1; $i < 3; $i++) {
            $picked = $predsq[3 * $ksqs[$i] + $i];
            if ($picked > $bestPredsq) {
                $bestPredsq = $picked;
                $isq = $i;
            }
        }
        $ksq = $ksqs[$isq];

        $stpsiz = $stplen[3 * $ksq + $isq];
        $ibd = $isbd[3 * $ksq + $isq];

        $xline = [];
        for ($i = 0; $i < $n; $i++) {
            $xline[] = Linalg::max(Linalg::min($xopt[$i] + $stpsiz * ($xpt[$i + $ksq * $n] - $xopt[$i]), $su[$i]), $sl[$i]);
        }
        if ($ibd < 0) {
            $xline[-$ibd - 1] = $sl[-$ibd - 1];
        }
        if ($ibd > 0) {
            $xline[$ibd - 1] = $su[$ibd - 1];
        }

        $d = [];
        for ($i = 0; $i < $n; $i++) {
            $d[] = $xline[$i] - $xopt[$i];
        }
        $denLine = PowAlg::calden($kopt, $bmat, $d, $xpt, $zmat, $n, $npt);

        if ($delbar > Util::F32_1E_2) {
            return $d;
        }

        $bigstp = $delbar + $delbar;
        $xcauchy = $xopt;
        $vlagsqCauchy = 0.0;
        $maskFree = \array_fill(0, $n, false);
        $maskFixl = \array_fill(0, $n, false);
        $maskFixu = \array_fill(0, $n, false);
        $xtemp = \array_fill(0, $n, 0.0);

        for ($uphill = 0; $uphill <= 1; $uphill++) {
            if ($uphill === 1) {
                for ($i = 0; $i < $n; $i++) {
                    $glag[$i] = -$glag[$i];
                }
            }

            $s = \array_fill(0, $n, 0.0);
            for ($i = 0; $i < $n; $i++) {
                $maskFree[$i] = Linalg::min($xopt[$i] - $sl[$i], $glag[$i]) > 0.0
                    || Linalg::max($xopt[$i] - $su[$i], $glag[$i]) < 0.0;
            }
            for ($i = 0; $i < $n; $i++) {
                if ($maskFree[$i]) {
                    $s[$i] = $bigstp;
                }
            }
            $ggfree = 0.0;
            for ($i = 0; $i < $n; $i++) {
                if ($maskFree[$i]) {
                    $ggfree += $glag[$i] * $glag[$i];
                }
            }
            if ($ggfree <= 0.0 || \is_nan($ggfree)) {
                continue;
            }

            $sfixsq = 0.0;
            $grdstp = 0.0;
            for ($pass = 0; $pass < $n; $pass++) {
                $resis = $delbar * $delbar - $sfixsq;
                if ($resis <= 0.0) {
                    break;
                }
                $ssqsav = $sfixsq;
                $grdstp = \sqrt($resis / $ggfree);
                for ($i = 0; $i < $n; $i++) {
                    $xtemp[$i] = $xopt[$i] - $grdstp * $glag[$i];
                }
                for ($i = 0; $i < $n; $i++) {
                    $maskFixl[$i] = $s[$i] >= $bigstp && $xtemp[$i] <= $sl[$i];
                    $maskFixu[$i] = $s[$i] >= $bigstp && $xtemp[$i] >= $su[$i];
                }
                for ($i = 0; $i < $n; $i++) {
                    $maskFree[$i] = $s[$i] >= $bigstp && ! $maskFixl[$i] && ! $maskFixu[$i];
                }
                for ($i = 0; $i < $n; $i++) {
                    if ($maskFixl[$i]) {
                        $s[$i] = $sl[$i] - $xopt[$i];
                    }
                    if ($maskFixu[$i]) {
                        $s[$i] = $su[$i] - $xopt[$i];
                    }
                }
                for ($i = 0; $i < $n; $i++) {
                    if ($maskFixl[$i] || $maskFixu[$i]) {
                        $sfixsq += $s[$i] * $s[$i];
                    }
                }
                $ggfree = 0.0;
                for ($i = 0; $i < $n; $i++) {
                    if ($maskFree[$i]) {
                        $ggfree += $glag[$i] * $glag[$i];
                    }
                }
                if (\is_nan($sfixsq) || \is_nan($ggfree) || $sfixsq <= $ssqsav || $ggfree <= 0.0) {
                    break;
                }
            }

            $x = \array_fill(0, $n, 0.0);
            for ($i = 0; $i < $n; $i++) {
                if ($glag[$i] > 0.0) {
                    $x[$i] = $sl[$i];
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if ($glag[$i] <= 0.0) {
                    $x[$i] = $su[$i];
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if (\abs($s[$i]) <= 0.0) {
                    $x[$i] = $xopt[$i];
                }
            }
            for ($i = 0; $i < $n; $i++) {
                $xtemp[$i] = Linalg::max(Linalg::min($xopt[$i] - $grdstp * $glag[$i], $su[$i]), $sl[$i]);
            }
            for ($i = 0; $i < $n; $i++) {
                if ($s[$i] >= $bigstp) {
                    $x[$i] = $xtemp[$i];
                }
            }
            for ($i = 0; $i < $n; $i++) {
                if ($s[$i] >= $bigstp) {
                    $s[$i] = -$grdstp * $glag[$i];
                }
            }
            $gs = Linalg::inprod($glag, $s);

            $sxpt = Linalg::matprod12($s, $xpt, $n, $npt);
            $curv = 0.0;
            for ($k = 0; $k < $npt; $k++) {
                $t = $pqlag[$k] * $sxpt[$k];
                $curv += $sxpt[$k] * $t;
            }
            if ($uphill === 1) {
                $curv = -$curv;
            }

            if ($curv > -$gs && $curv < -(1.0 + \sqrt(2.0)) * $gs) {
                $scaling = -$gs / $curv;
                for ($i = 0; $i < $n; $i++) {
                    $x[$i] = Linalg::max(Linalg::min($xopt[$i] + $scaling * $s[$i], $su[$i]), $sl[$i]);
                }
                $halfGsScaling = 0.5 * $gs * $scaling;
                $vlagsq = $halfGsScaling * $halfGsScaling;
            } else {
                $t = $gs + 0.5 * $curv;
                $vlagsq = $t * $t;
            }

            if ($vlagsq > $vlagsqCauchy) {
                $xcauchy = $x;
                $vlagsqCauchy = $vlagsq;
            }
        }

        $sCauchy = [];
        for ($i = 0; $i < $n; $i++) {
            $sCauchy[] = $xcauchy[$i] - $xopt[$i];
        }
        $denCauchy = PowAlg::calden($kopt, $bmat, $sCauchy, $xpt, $zmat, $n, $npt);
        if ($denCauchy[$knew] > Linalg::max($denLine[$knew], 0.0) || \is_nan($denLine[$knew])) {
            $d = $sCauchy;
        }

        $dAbsSum = Linalg::sumAbs($d);
        if ($dAbsSum <= 0.0 || ! \is_finite($dAbsSum)) {
            return self::fallbackStep($xpt, $xopt, $knew, $delbar, $n);
        }

        return $d;
    }

    /**
     * The scaled displacement from XOPT toward XPT(:, knew), PRIMA's non-finite fallback.
     *
     * @param  list<float>  $xpt
     * @param  list<float>  $xopt
     * @return list<float>
     */
    private static function fallbackStep(array $xpt, array $xopt, int $knew, float $delbar, int $n): array
    {
        $d = [];
        for ($i = 0; $i < $n; $i++) {
            $d[] = $xpt[$i + $knew * $n] - $xopt[$i];
        }
        $scale = Linalg::min(0.5, $delbar / Linalg::norm($d));
        for ($i = 0; $i < $n; $i++) {
            $d[$i] *= $scale;
        }

        return $d;
    }
}
