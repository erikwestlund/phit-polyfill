<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal\Bobyqa;

/**
 * The crate's initialize.rs: PRIMA's initxf, initq and inith.
 *
 * @internal
 */
final class Initialize
{
    /**
     * Initialize the interpolation points and their function values.
     *
     * @param  callable(list<float>): (float|int)  $calfun
     * @param  list<float>  $xl
     * @param  list<float>  $xu
     * @param  list<float>  $x0
     * @param  list<array{int, int}>  $ij
     * @param  list<float>  $fval
     * @param  list<float>  $sl
     * @param  list<float>  $su
     * @param  list<float>  $xbase
     * @param  list<float>  $xpt
     * @return array{int, int, int} [kopt, nf, info]
     */
    public static function initxf(
        callable $calfun,
        int $maxfun,
        float $ftarget,
        float $rhobeg,
        array $xl,
        array $xu,
        array &$x0,
        array &$ij,
        array &$fval,
        array &$sl,
        array &$su,
        array &$xbase,
        array &$xpt,
        int $n,
        int $npt,
    ): array {
        $info = Util::INFO_DFT;

        for ($i = 0; $i < $n; $i++) {
            $sl[$i] = $xl[$i] - $x0[$i];
            $su[$i] = $xu[$i] - $x0[$i];
        }
        for ($i = 0; $i < $n; $i++) {
            if ($sl[$i] < 0.0) {
                $sl[$i] = Linalg::min($sl[$i], -$rhobeg);
            } else {
                $x0[$i] = $xl[$i];
                $sl[$i] = 0.0;
                $su[$i] = $xu[$i] - $xl[$i];
            }
        }
        for ($i = 0; $i < $n; $i++) {
            if ($su[$i] > 0.0) {
                $su[$i] = Linalg::max($su[$i], $rhobeg);
            } else {
                $x0[$i] = $xu[$i];
                $sl[$i] = $xl[$i] - $xu[$i];
                $su[$i] = 0.0;
            }
        }

        $xbase = $x0;
        $evaluated = array_fill(0, $npt, false);
        $fval = array_fill(0, $npt, Linalg::REALMAX);

        $xpt = array_fill(0, $n * $npt, 0.0);
        for ($k = 0; $k < $n; $k++) {
            $xpt[$k + ($k + 1) * $n] = $rhobeg;
            if ($su[$k] <= 0.0) {
                $xpt[$k + ($k + 1) * $n] = -$rhobeg;
            }
        }
        $limit = min($npt - $n - 1, $n);
        for ($k = 0; $k < $limit; $k++) {
            $index = $k + ($k + $n + 1) * $n;
            $xpt[$index] = -$rhobeg;
            if ($sl[$k] >= 0.0) {
                $xpt[$index] = Linalg::min(2.0 * $rhobeg, $su[$k]);
            }
            if ($su[$k] <= 0.0) {
                $xpt[$index] = Linalg::max(-2.0 * $rhobeg, $sl[$k]);
            }
        }

        $limit = min($npt, 2 * $n + 1);
        for ($k = 0; $k < $limit; $k++) {
            $x = Util::xinbd($xbase, Linalg::col($xpt, $n, $k), $xl, $xu, $sl, $su);
            $f = Util::evaluate($calfun, $x);
            $evaluated[$k] = true;
            $fval[$k] = $f;
            $subinfo = Util::checkexit($maxfun, $k + 1, $f, $ftarget, $x);
            if ($subinfo !== Util::INFO_DFT) {
                $info = $subinfo;
                break;
            }
        }

        $limit = min($npt - $n, $n + 1);
        for ($k = 1; $k < $limit; $k++) {
            if ($xpt[$k - 1 + $k * $n] * $xpt[$k - 1 + ($k + $n) * $n] < 0.0 && $fval[$k + $n] < $fval[$k]) {
                [$fval[$k], $fval[$k + $n]] = [$fval[$k + $n], $fval[$k]];
                for ($r = 0; $r < $n; $r++) {
                    $tmp = $xpt[$r + $k * $n];
                    $xpt[$r + $k * $n] = $xpt[$r + ($k + $n) * $n];
                    $xpt[$r + ($k + $n) * $n] = $tmp;
                }
            }
        }

        $ij = PowAlg::setij($n, $npt);
        foreach ($ij as $k => [$i, $j]) {
            for ($r = 0; $r < $n; $r++) {
                $xpt[$r + (2 * $n + 1 + $k) * $n] = $xpt[$r + ($i + 1) * $n] + $xpt[$r + ($j + 1) * $n];
            }
        }

        if ($info === Util::INFO_DFT) {
            for ($k = 2 * $n + 1; $k < $npt; $k++) {
                $x = Util::xinbd($xbase, Linalg::col($xpt, $n, $k), $xl, $xu, $sl, $su);
                $f = Util::evaluate($calfun, $x);
                $evaluated[$k] = true;
                $fval[$k] = $f;
                $subinfo = Util::checkexit($maxfun, $k + 1, $f, $ftarget, $x);
                if ($subinfo !== Util::INFO_DFT) {
                    $info = $subinfo;
                    break;
                }
            }
        }

        $nf = 0;
        $kopt = 0;
        $found = false;
        for ($k = 0; $k < $npt; $k++) {
            if ($evaluated[$k]) {
                $nf++;
                if (! $found || $fval[$k] < $fval[$kopt]) {
                    $kopt = $k;
                    $found = true;
                }
            }
        }

        return [$kopt, $nf, $info];
    }

    /**
     * Initialize the quadratic model [GOPT, HQ, PQ].
     *
     * @param  list<array{int, int}>  $ij
     * @param  list<float>  $fval
     * @param  list<float>  $xpt
     * @param  list<float>  $gopt
     * @param  list<float>  $hq
     * @param  list<float>  $pq
     */
    public static function initq(array $ij, array $fval, array $xpt, array &$gopt, array &$hq, array &$pq, int $n, int $npt): void
    {
        $fbase = $fval[0];
        for ($i = 0; $i < $n; $i++) {
            $gopt[$i] = ($fval[$i + 1] - $fbase) / $xpt[$i + ($i + 1) * $n];
        }

        $ndiag = min($n, $npt - $n - 1);
        $xa = [];
        $xb = [];
        for ($k = 0; $k < $ndiag; $k++) {
            $xa[$k] = $xpt[$k + ($k + 1) * $n];
            $xb[$k] = $xpt[$k + ($n + 1 + $k) * $n];
        }

        for ($k = 0; $k < $ndiag; $k++) {
            $gopt[$k] = ($gopt[$k] * $xb[$k] - (($fval[$n + 1 + $k] - $fbase) / $xb[$k]) * $xa[$k]) / ($xb[$k] - $xa[$k]);
        }

        $hq = array_fill(0, $n * $n, 0.0);
        for ($k = 0; $k < $ndiag; $k++) {
            $hq[$k + $k * $n] = 2.0 * (($fval[$k + 1] - $fbase) / $xa[$k] - ($fval[$n + $k + 1] - $fbase) / $xb[$k])
                / ($xa[$k] - $xb[$k]);
        }

        foreach ($ij as $k => [$i, $j]) {
            $xi = $xpt[$i + (2 * $n + 1 + $k) * $n];
            $xj = $xpt[$j + (2 * $n + 1 + $k) * $n];
            $hq[$i + $j * $n] = ($fbase - $fval[$i + 1] - $fval[$j + 1] + $fval[2 * $n + 1 + $k]) / ($xi * $xj);
            $hq[$j + $i * $n] = $hq[$i + $j * $n];
        }

        $kopt = 0;
        for ($k = 0; $k < $npt; $k++) {
            if ($fval[$k] < $fval[$kopt]) {
                $kopt = $k;
            }
        }
        if ($kopt !== 0) {
            $shift = Linalg::matprod21($hq, $n, $n, Linalg::col($xpt, $n, $kopt));
            for ($i = 0; $i < $n; $i++) {
                $gopt[$i] += $shift[$i];
            }
        }

        $pq = array_fill(0, $npt, 0.0);
    }

    /**
     * Initialize [BMAT, ZMAT], the matrix H of (2.7) of the BOBYQA paper.
     *
     * @param  list<array{int, int}>  $ij
     * @param  list<float>  $xpt
     * @param  list<float>  $bmat
     * @param  list<float>  $zmat
     */
    public static function inith(array $ij, array $xpt, array &$bmat, array &$zmat, int $n, int $npt): void
    {
        $rhobeg = 0.0;
        for ($i = 0; $i < $n; $i++) {
            $rhobeg = Linalg::max($rhobeg, abs($xpt[$i + $n]));
        }
        $rhosq = $rhobeg * $rhobeg;

        $ndiag = min($n, $npt - $n - 1);
        $xa = [];
        $xb = [];
        for ($k = 0; $k < $ndiag; $k++) {
            $xa[$k] = $xpt[$k + ($k + 1) * $n];
            $xb[$k] = $xpt[$k + ($n + 1 + $k) * $n];
        }

        $bmat = array_fill(0, $n * ($npt + $n), 0.0);
        for ($k = 0; $k < $ndiag; $k++) {
            $bmat[$k] = -($xa[$k] + $xb[$k]) / ($xa[$k] * $xb[$k]);
        }
        for ($k = 0; $k < $ndiag; $k++) {
            $bmat[$k + ($k + $n + 1) * $n] = -0.5 / $xpt[$k + ($k + 1) * $n];
            $bmat[$k + ($k + 1) * $n] = -$bmat[$k] - $bmat[$k + ($k + $n + 1) * $n];
        }
        for ($k = $ndiag; $k < $n; $k++) {
            $bmat[$k] = -1.0 / $xpt[$k + ($k + 1) * $n];
            $bmat[$k + ($k + 1) * $n] = -$bmat[$k];
            $bmat[$k + ($npt + $k) * $n] = -0.5 * $rhosq;
        }

        $nz = $npt - $n - 1;
        $zmat = array_fill(0, $npt * $nz, 0.0);
        for ($k = 0; $k < $ndiag; $k++) {
            $zmat[$k * $npt] = sqrt(2.0) / ($xa[$k] * $xb[$k]);
        }
        for ($k = 0; $k < $ndiag; $k++) {
            $zmat[$k + 1 + $k * $npt] = -$zmat[$k * $npt] - sqrt(0.5) / $rhosq;
            $zmat[$k + $n + 1 + $k * $npt] = sqrt(0.5) / $rhosq;
        }
        for ($k = $ndiag; $k < $nz; $k++) {
            $zmat[$k * $npt] = 1.0 / $rhosq;
            $zmat[$k + $n + 1 + $k * $npt] = 1.0 / $rhosq;
            [$i, $j] = $ij[$k - $n];
            $zmat[$i + 1 + $k * $npt] = -1.0 / $rhosq;
            $zmat[$j + 1 + $k * $npt] = -1.0 / $rhosq;
        }
    }
}
