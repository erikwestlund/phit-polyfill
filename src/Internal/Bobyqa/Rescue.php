<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * The crate's rescue.rs: PRIMA's RESCUE (Section 5 of the BOBYQA paper) and
 * its private updateh_rsc. The ip/iq values decoded from PTSID stay 1-based
 * with 0 meaning "no direction", as in the crate.
 *
 * @internal
 */
final class Rescue
{
    /**
     * Replace a few interpolation points to restore the geometry; returns info.
     *
     * @param  callable(list<float>): (float|int)  $calfun
     * @param  list<float>  $xl
     * @param  list<float>  $xu
     * @param  list<float>  $fval
     * @param  list<float>  $gopt
     * @param  list<float>  $hq
     * @param  list<float>  $pq
     * @param  list<float>  $sl
     * @param  list<float>  $su
     * @param  list<float>  $xbase
     * @param  list<float>  $xpt
     * @param  list<float>  $bmat
     * @param  list<float>  $zmat
     */
    public static function rescue(
        callable $calfun,
        int $maxfun,
        float $delta,
        float $ftarget,
        array $xl,
        array $xu,
        int &$kopt,
        int &$nf,
        array &$fval,
        array &$gopt,
        array &$hq,
        array &$pq,
        array &$sl,
        array &$su,
        array &$xbase,
        array &$xpt,
        array &$bmat,
        array &$zmat,
        int $n,
        int $npt,
    ): int {
        $nz = $npt - $n - 1;
        $info = Util::INFO_DFT;

        if ($nf >= $maxfun) {
            $bmat = \array_fill(0, $n * ($npt + $n), 0.0);
            $zmat = \array_fill(0, $npt * $nz, 0.0);

            return Util::MAXFUN_REACHED;
        }

        $xopt = Linalg::col($xpt, $n, $kopt);
        for ($i = 0; $i < $n; $i++) {
            $sl[$i] = Linalg::min($sl[$i] - $xopt[$i], 0.0);
            $su[$i] = Linalg::max($su[$i] - $xopt[$i], 0.0);
            $xbase[$i] = Linalg::min(Linalg::max($xl[$i], $xbase[$i] + $xopt[$i]), $xu[$i]);
        }
        for ($k = 0; $k < $npt; $k++) {
            for ($i = 0; $i < $n; $i++) {
                $xpt[$i + $k * $n] -= $xopt[$i];
            }
        }
        for ($i = 0; $i < $n; $i++) {
            $xpt[$i + $kopt * $n] = 0.0;
        }

        $pqsum = 0.0;
        foreach ($pq as $p) {
            $pqsum += $p;
        }
        $halfsum = 0.5 * $pqsum;
        $v = Linalg::matprod21($xpt, $n, $npt, $pq);
        for ($i = 0; $i < $n; $i++) {
            $v[$i] += $halfsum * $xopt[$i];
        }
        Linalg::r2update($hq, $n, 1.0, $xopt, $v);

        $ptsaux1 = [];
        $ptsaux2 = [];
        for ($i = 0; $i < $n; $i++) {
            $ptsaux1[] = Linalg::min($delta, $su[$i]);
            $ptsaux2[] = Linalg::max(-$delta, $sl[$i]);
        }
        for ($i = 0; $i < $n; $i++) {
            if ($ptsaux1[$i] + $ptsaux2[$i] < 0.0) {
                [$ptsaux1[$i], $ptsaux2[$i]] = [$ptsaux2[$i], $ptsaux1[$i]];
            }
        }
        for ($i = 0; $i < $n; $i++) {
            if (\abs($ptsaux2[$i]) < 0.5 * \abs($ptsaux1[$i])) {
                $ptsaux2[$i] = 0.5 * $ptsaux1[$i];
            }
        }

        $sfrac = 0.5 / ((float) $n + 1.0);
        $ptsid = \array_fill(0, $npt, 0.0);
        $ptsid[0] = $sfrac;
        $bmat = \array_fill(0, $n * ($npt + $n), 0.0);
        $zmat = \array_fill(0, $npt * $nz, 0.0);
        for ($k = 0; $k < $n; $k++) {
            $ptsid[$k + 1] = (float) ($k + 1) + $sfrac;
            if ($k + 1 <= $nz) {
                $ptsid[$k + $n + 1] = (float) ($k + 1) / ((float) $n + 1.0) + $sfrac;
                $temp = 1.0 / ($ptsaux1[$k] - $ptsaux2[$k]);
                $bmat[$k + ($k + 1) * $n] = -$temp + 1.0 / $ptsaux1[$k];
                $bmat[$k + ($k + $n + 1) * $n] = $temp + 1.0 / $ptsaux2[$k];
                $bmat[$k] = -$bmat[$k + ($k + 1) * $n] - $bmat[$k + ($k + $n + 1) * $n];
                $zmat[$k * $npt] = \sqrt(2.0) / \abs($ptsaux1[$k] * $ptsaux2[$k]);
                $zmat[$k + 1 + $k * $npt] = $zmat[$k * $npt] * $ptsaux2[$k] * $temp;
                $zmat[$k + $n + 1 + $k * $npt] = -$zmat[$k * $npt] * $ptsaux1[$k] * $temp;
            } else {
                $bmat[$k] = -1.0 / $ptsaux1[$k];
                $bmat[$k + ($k + 1) * $n] = 1.0 / $ptsaux1[$k];
                $bmat[$k + ($k + $npt) * $n] = -0.5 * $ptsaux1[$k] * $ptsaux1[$k];
            }
        }

        $ij = PowAlg::setij($n, $npt);
        for ($k = 2 * $n + 1; $k < $npt; $k++) {
            [$i0, $j0] = $ij[$k - (2 * $n + 1)];
            $ip = $i0 + 1;
            $iq = $j0 + 1;
            $ptsid[$k] = (float) $ip + (float) $iq / ((float) $n + 1.0) + $sfrac;
            $temp = 1.0 / ($ptsaux1[$ip - 1] * $ptsaux1[$iq - 1]);
            $col = ($k - $n - 1) * $npt;
            $zmat[$col] = $temp;
            $zmat[$k + $col] = $temp;
            $zmat[$ip + $col] = -$temp;
            $zmat[$iq + $col] = -$temp;
        }

        if ($kopt !== 0) {
            for ($i = 0; $i < $n; $i++) {
                $tmp = $bmat[$i];
                $bmat[$i] = $bmat[$i + $kopt * $n];
                $bmat[$i + $kopt * $n] = $tmp;
            }
            for ($j = 0; $j < $nz; $j++) {
                $tmp = $zmat[$j * $npt];
                $zmat[$j * $npt] = $zmat[$kopt + $j * $npt];
                $zmat[$kopt + $j * $npt] = $tmp;
            }
        }
        $ptsid[0] = $ptsid[$kopt];
        $ptsid[$kopt] = 0.0;

        $score = [];
        for ($k = 0; $k < $npt; $k++) {
            $s = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $s += $xpt[$i + $k * $n] * $xpt[$i + $k * $n];
            }
            $score[] = \sqrt($s);
        }
        $score[$kopt] = 0.0;
        $scoreinc = 0.0;
        foreach ($score as $s) {
            $scoreinc = Linalg::max($scoreinc, $s);
        }

        $nprov = $npt - 1;
        $maxiter = $npt * $npt;
        $vlag = \array_fill(0, $npt + $n, 0.0);
        for ($iter = 0; $iter < $maxiter; $iter++) {
            $allNonPositive = true;
            foreach ($score as $s) {
                if (\is_nan($s) || $s > 0.0) {
                    $allNonPositive = false;
                    break;
                }
            }
            if ($allNonPositive || $nprov <= 1) {
                break;
            }

            $korig = 0;
            $found = false;
            for ($k = 0; $k < $npt; $k++) {
                if ($score[$k] > 0.0 && (! $found || $score[$k] < $score[$korig])) {
                    $korig = $k;
                    $found = true;
                }
            }

            $wmv = [];
            $korigOff = $korig * $n;
            for ($k = 0; $k < $npt; $k++) {
                if ($k === $kopt) {
                    $w = 0.0;
                } elseif ($ptsid[$k] <= 0.0) {
                    $w = 0.0;
                    for ($i = 0; $i < $n; $i++) {
                        $w += $xpt[$korigOff + $i] * $xpt[$i + $k * $n];
                    }
                } else {
                    [$ip, $iq] = self::decode($ptsid[$k], $n);
                    if ($ip > 0 && $iq > 0) {
                        $w = $xpt[$ip - 1 + $korigOff] * $ptsaux1[$ip - 1] + $xpt[$iq - 1 + $korigOff] * $ptsaux1[$iq - 1];
                    } elseif ($ip > 0) {
                        $w = $xpt[$ip - 1 + $korigOff] * $ptsaux1[$ip - 1];
                    } elseif ($iq > 0) {
                        $w = $xpt[$iq - 1 + $korigOff] * $ptsaux2[$iq - 1];
                    } else {
                        $w = 0.0;
                    }
                }
                $wmv[] = 0.5 * $w * $w;
            }
            for ($i = 0; $i < $n; $i++) {
                $wmv[] = $xpt[$korigOff + $i];
            }
            $wmvHead = \array_slice($wmv, 0, $npt);
            $wmvTail = \array_slice($wmv, $npt, $n);

            $wmvZ = Linalg::matprod12($wmvHead, $zmat, $npt, $nz);
            $zWmvZ = Linalg::matprod21($zmat, $npt, $nz, $wmvZ);
            for ($k = 0; $k < $npt; $k++) {
                $s = 0.0;
                for ($i = 0; $i < $n; $i++) {
                    $s += $wmvTail[$i] * $bmat[$i + $k * $n];
                }
                $vlag[$k] = $zWmvZ[$k] + $s;
            }
            $bmatWmv = Linalg::matprod21($bmat, $n, $npt + $n, $wmv);
            for ($i = 0; $i < $n; $i++) {
                $vlag[$npt + $i] = $bmatWmv[$i];
            }

            $t = \array_fill(0, $n, 0.0);
            for ($j = 0; $j < $npt; $j++) {
                for ($i = 0; $i < $n; $i++) {
                    $t[$i] += $bmat[$i + $j * $n] * $wmv[$j];
                }
            }
            for ($i = 0; $i < $n; $i++) {
                $t[$i] += $bmatWmv[$i];
            }
            // wmv(1:n) is PRIMA's literal (asymmetric-looking) operand.
            $bsum = Linalg::inprod(\array_slice($wmv, 0, $n), $t);
            $s = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $s += $xpt[$korigOff + $i] * $xpt[$korigOff + $i];
            }
            $wmvzSq = 0.0;
            foreach ($wmvZ as $w) {
                $wmvzSq += $w * $w;
            }
            $beta = 0.5 * ($s * $s) - $wmvzSq - $bsum;

            $vlag[$kopt] += 1.0;

            $den = \array_fill(0, $npt, 0.0);
            for ($k = 0; $k < $npt; $k++) {
                if ($ptsid[$k] > 0.0) {
                    $hdiag = 0.0;
                    for ($j = 0; $j < $nz; $j++) {
                        $z = $zmat[$k + $j * $npt];
                        $hdiag += $z * $z;
                    }
                    $den[$k] = $hdiag * $beta + $vlag[$k] * $vlag[$k];
                }
            }

            $vlagAbsSum = Linalg::sumAbs($vlag);
            $vmax = 0.0;
            for ($k = 0; $k < $npt; $k++) {
                $vmax = Linalg::max($vmax, $vlag[$k] * $vlag[$k]);
            }
            $anyDen = false;
            $threshold = 5.0e-2 * $vmax;
            foreach ($den as $dk) {
                if ($dk > $threshold) {
                    $anyDen = true;
                    break;
                }
            }
            if (! \is_finite($vlagAbsSum) || ! $anyDen) {
                $score[$korig] = -$score[$korig] - $scoreinc;

                continue;
            }

            $kprov = 0;
            $kfound = false;
            for ($k = 0; $k < $npt; $k++) {
                if (! \is_nan($den[$k]) && (! $kfound || $den[$k] > $den[$kprov])) {
                    $kprov = $k;
                    $kfound = true;
                }
            }

            if ($kprov !== $korig) {
                for ($i = 0; $i < $n; $i++) {
                    $tmp = $bmat[$i + $kprov * $n];
                    $bmat[$i + $kprov * $n] = $bmat[$i + $korig * $n];
                    $bmat[$i + $korig * $n] = $tmp;
                }
                for ($j = 0; $j < $nz; $j++) {
                    $tmp = $zmat[$kprov + $j * $npt];
                    $zmat[$kprov + $j * $npt] = $zmat[$korig + $j * $npt];
                    $zmat[$korig + $j * $npt] = $tmp;
                }
                [$vlag[$kprov], $vlag[$korig]] = [$vlag[$korig], $vlag[$kprov]];
            }
            $ptsid[$kprov] = $ptsid[$korig];
            $ptsid[$korig] = 0.0;
            $score[$korig] = 0.0;
            foreach ($score as $k => $s) {
                $score[$k] = \abs($s);
            }

            self::updatehRsc($korig, $beta, $vlag, $bmat, $zmat, $n, $npt);

            $nprov--;
        }

        $kbase = $kopt;
        $fbase = $fval[$kopt];
        if ($nprov > 0) {
            for ($kpt = 0; $kpt < $npt; $kpt++) {
                if ($ptsid[$kpt] <= 0.0) {
                    continue;
                }

                Linalg::r1update($hq, $n, $pq[$kpt], Linalg::col($xpt, $n, $kpt));
                $pq[$kpt] = 0.0;

                [$ip, $iq] = self::decode($ptsid[$kpt], $n);

                $xp = 0.0;
                $xq = 0.0;
                $xnew = \array_fill(0, $n, 0.0);
                if ($ip > 0 && $iq > 0) {
                    $xp = $ptsaux1[$ip - 1];
                    $xnew[$ip - 1] = $xp;
                    $xq = $ptsaux1[$iq - 1];
                    $xnew[$iq - 1] = $xq;
                } elseif ($ip > 0) {
                    $xp = $ptsaux1[$ip - 1];
                    $xnew[$ip - 1] = $xp;
                } elseif ($iq > 0) {
                    $xq = $ptsaux2[$iq - 1];
                    $xnew[$iq - 1] = $xq;
                }

                $diffAbsSum = 0.0;
                $xnewAbsSum = 0.0;
                for ($i = 0; $i < $n; $i++) {
                    $diffAbsSum += \abs($xnew[$i] - $xpt[$i + $kpt * $n]);
                    $xnewAbsSum += \abs($xnew[$i]);
                }
                if ($diffAbsSum <= Util::F32_1E_2 * $delta || ! \is_finite($xnewAbsSum)) {
                    continue;
                }
                for ($i = 0; $i < $n; $i++) {
                    $xpt[$i + $kpt * $n] = $xnew[$i];
                }

                $x = Util::xinbd($xbase, $xnew, $xl, $xu, $sl, $su);
                $f = Util::evaluate($calfun, $x);
                $nf++;

                $fval[$kpt] = $f;
                if ($f < $fval[$kopt]) {
                    $kopt = $kpt;
                }

                $subinfo = Util::checkexit($maxfun, $nf, $f, $ftarget, $x);
                if ($subinfo !== Util::INFO_DFT) {
                    $info = $subinfo;
                    break;
                }

                $vquad = $fbase;
                $xxpt = \array_fill(0, $npt, 0.0);
                if ($ip > 0 && $iq > 0) {
                    $vquad += $xp * ($gopt[$ip - 1] + 0.5 * $xp * $hq[$ip - 1 + ($ip - 1) * $n]);
                    $vquad += $xq * ($gopt[$iq - 1] + 0.5 * $xq * $hq[$iq - 1 + ($iq - 1) * $n]);
                    $vquad += $xp * $xq * $hq[$ip - 1 + ($iq - 1) * $n];
                    for ($k = 0; $k < $npt; $k++) {
                        $xxpt[$k] = $xp * $xpt[$ip - 1 + $k * $n] + $xq * $xpt[$iq - 1 + $k * $n];
                    }
                } elseif ($ip > 0) {
                    $vquad += $xp * ($gopt[$ip - 1] + 0.5 * $xp * $hq[$ip - 1 + ($ip - 1) * $n]);
                    for ($k = 0; $k < $npt; $k++) {
                        $xxpt[$k] = $xp * $xpt[$ip - 1 + $k * $n];
                    }
                } elseif ($iq > 0) {
                    $vquad += $xq * ($gopt[$iq - 1] + 0.5 * $xq * $hq[$iq - 1 + ($iq - 1) * $n]);
                    for ($k = 0; $k < $npt; $k++) {
                        $xxpt[$k] = $xq * $xpt[$iq - 1 + $k * $n];
                    }
                }
                $pqXxpt = [];
                for ($k = 0; $k < $npt; $k++) {
                    $pqXxpt[] = $pq[$k] * $xxpt[$k];
                }
                $vquad += 0.5 * Linalg::inprod($xxpt, $pqXxpt);

                $moderr = $f - $vquad;
                for ($i = 0; $i < $n; $i++) {
                    $gopt[$i] += $moderr * $bmat[$i + $kpt * $n];
                }
                $zrow = [];
                for ($j = 0; $j < $nz; $j++) {
                    $zrow[] = $zmat[$kpt + $j * $npt];
                }
                $zZrow = Linalg::matprod21($zmat, $npt, $nz, $zrow);
                $pqinc = [];
                for ($k = 0; $k < $npt; $k++) {
                    $pqinc[] = $moderr * $zZrow[$k];
                }
                for ($k = 0; $k < $npt; $k++) {
                    if ($ptsid[$k] <= 0.0) {
                        $pq[$k] += $pqinc[$k];
                    }
                }
                for ($k = 0; $k < $npt; $k++) {
                    if ($ptsid[$k] <= 0.0) {
                        continue;
                    }
                    [$ipk, $iqk] = self::decode($ptsid[$k], $n);
                    if ($ipk > 0 && $iqk > 0) {
                        $pp = $ipk - 1 + ($ipk - 1) * $n;
                        $qq = $iqk - 1 + ($iqk - 1) * $n;
                        $pqIndex = $ipk - 1 + ($iqk - 1) * $n;
                        $hq[$pp] += $pqinc[$k] * ($ptsaux1[$ipk - 1] * $ptsaux1[$ipk - 1]);
                        $hq[$qq] += $pqinc[$k] * ($ptsaux1[$iqk - 1] * $ptsaux1[$iqk - 1]);
                        $hq[$pqIndex] += $pqinc[$k] * $ptsaux1[$ipk - 1] * $ptsaux1[$iqk - 1];
                        $hq[$iqk - 1 + ($ipk - 1) * $n] = $hq[$pqIndex];
                    } elseif ($ipk > 0) {
                        $hq[$ipk - 1 + ($ipk - 1) * $n] += $pqinc[$k] * ($ptsaux1[$ipk - 1] * $ptsaux1[$ipk - 1]);
                    } elseif ($iqk > 0) {
                        $hq[$iqk - 1 + ($iqk - 1) * $n] += $pqinc[$k] * ($ptsaux2[$iqk - 1] * $ptsaux2[$iqk - 1]);
                    }
                }
                $ptsid[$kpt] = 0.0;
            }
        }

        if ($kopt !== $kbase) {
            $shift = PowAlg::hessMul(Linalg::col($xpt, $n, $kopt), $xpt, $n, $npt, $pq, $hq);
            for ($i = 0; $i < $n; $i++) {
                $gopt[$i] += $shift[$i];
            }
        }

        return $info;
    }

    /**
     * Decode PTSID into PRIMA's 1-based (ip, iq), 0 meaning none.
     *
     * @return array{int, int}
     */
    private static function decode(float $ptsid, int $n): array
    {
        $ip = (int) \floor($ptsid);
        $iq = (int) \floor(((float) $n + 1.0) * $ptsid - (float) (($n + 1) * $ip));

        return [$ip, $iq];
    }

    /**
     * UPDATEH driven by a precomputed (BETA, VLAG).
     *
     * @param  list<float>  $vlagIn
     * @param  list<float>  $bmat
     * @param  list<float>  $zmat
     */
    private static function updatehRsc(int $knew, float $beta, array $vlagIn, array &$bmat, array &$zmat, int $n, int $npt): int
    {
        $nz = $npt - $n - 1;
        $vlag = $vlagIn;
        $tau = $vlag[$knew];
        $zknewSq = 0.0;
        for ($j = 0; $j < $nz; $j++) {
            $z = $zmat[$knew + $j * $npt];
            $zknewSq += $z * $z;
        }
        $denom = $zknewSq * $beta + $tau * $tau;

        $guard = Linalg::sumAbs($vlag) + \abs($beta);
        if (! \is_finite($guard) || \is_nan($denom) || $denom <= 0.0) {
            return Util::DAMAGING_ROUNDING;
        }

        $vlag[$knew] -= 1.0;

        Update::rotateZmatRow($zmat, $knew, $npt, $nz);

        $hcol = [];
        $zknew0 = $zmat[$knew];
        for ($k = 0; $k < $npt; $k++) {
            $hcol[] = $zknew0 * $zmat[$k];
        }
        for ($i = 0; $i < $n; $i++) {
            $hcol[] = $bmat[$i + $knew * $n];
        }

        $sqrtdn = \sqrt($denom);
        $zknew1 = $zmat[$knew] / $sqrtdn;
        $tauSqrtdn = $tau / $sqrtdn;
        for ($k = 0; $k < $npt; $k++) {
            $zmat[$k] = $tauSqrtdn * $zmat[$k] - $zknew1 * $vlag[$k];
        }
        $zmat[$knew] = $zknew1;

        $alpha = $hcol[$knew];
        $v1 = [];
        $v2 = [];
        for ($i = 0; $i < $n; $i++) {
            $v1[] = ($alpha * $vlag[$npt + $i] - $tau * $hcol[$npt + $i]) / $denom;
            $v2[] = (-$beta * $hcol[$npt + $i] - $tau * $vlag[$npt + $i]) / $denom;
        }
        $total = $npt + $n;
        for ($j = 0; $j < $total; $j++) {
            for ($i = 0; $i < $n; $i++) {
                $index = $i + $j * $n;
                $bmat[$index] = $bmat[$index] + $v1[$i] * $vlag[$j] + $v2[$i] * $hcol[$j];
            }
        }
        for ($j = 0; $j < $n; $j++) {
            for ($i = 0; $i < $j; $i++) {
                $bmat[$i + ($npt + $j) * $n] = $bmat[$j + ($npt + $i) * $n];
            }
        }

        return Util::INFO_DFT;
    }
}
