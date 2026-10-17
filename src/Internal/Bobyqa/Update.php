<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal\Bobyqa;

/**
 * The crate's update.rs: PRIMA's updateh, updatexf, updateq and tryqalt, the
 * updates when XPT(:, KNEW) becomes XOPT + D.
 *
 * @internal
 */
final class Update
{
    /**
     * Update BMAT and ZMAT; returns INFO_DFT or DAMAGING_ROUNDING.
     *
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $bmat
     * @param  list<float>  $zmat
     * @param  array{list<float>, float}|null  $precomputed  VLAG (with the +1) and BETA from this iteration
     */
    public static function updateh(
        ?int $knew,
        int $kopt,
        array $d,
        array $xpt,
        array &$bmat,
        array &$zmat,
        int $n,
        int $npt,
        ?array $precomputed,
    ): int {
        if ($knew === null) {
            return Util::INFO_DFT;
        }
        $nz = $npt - $n - 1;

        $zrow = [];
        for ($j = 0; $j < $nz; $j++) {
            $zrow[] = $zmat[$knew + $j * $npt];
        }
        $hcol = Linalg::matprod21($zmat, $npt, $nz, $zrow);
        for ($i = 0; $i < $n; $i++) {
            $hcol[] = $bmat[$i + $knew * $n];
        }

        if ($precomputed !== null) {
            [$vlag, $beta] = $precomputed;
        } else {
            $beta = PowAlg::calbeta($kopt, $bmat, $d, $xpt, $zmat, $n, $npt);
            $vlag = PowAlg::calvlag($kopt, $bmat, $d, $xpt, $zmat, $n, $npt);
        }

        $alpha = $hcol[$knew];
        $tau = $vlag[$knew];
        $denom = $alpha * $beta + $tau * $tau;

        $vlag[$knew] -= 1.0;

        $guard = Linalg::sumAbs($hcol) + Linalg::sumAbs($vlag) + abs($beta);
        if (! is_finite($guard) || is_nan($denom) || $denom <= 0.0) {
            return Util::DAMAGING_ROUNDING;
        }

        $v1 = [];
        $v2 = [];
        for ($i = 0; $i < $n; $i++) {
            $v1[] = ($alpha * $vlag[$npt + $i] - $tau * $hcol[$npt + $i]) / $denom;
            $v2[] = (-$beta * $hcol[$npt + $i] - $tau * $vlag[$npt + $i]) / $denom;
        }

        $total = $npt + $n;
        for ($j = 0; $j < $total; $j++) {
            $vj = $vlag[$j];
            $hj = $hcol[$j];
            $offset = $j * $n;
            for ($i = 0; $i < $n; $i++) {
                $bmat[$offset + $i] = ($bmat[$offset + $i] + $v1[$i] * $vj) + $v2[$i] * $hj;
            }
        }

        for ($j = 0; $j < $n; $j++) {
            for ($i = 0; $i < $j; $i++) {
                $bmat[$i + ($npt + $j) * $n] = $bmat[$j + ($npt + $i) * $n];
            }
        }

        self::rotateZmatRow($zmat, $knew, $npt, $nz);

        $sqrtdn = sqrt($denom);
        $tauQ = $tau / $sqrtdn;
        $zkQ = $zmat[$knew] / $sqrtdn;
        for ($i = 0; $i < $npt; $i++) {
            $zmat[$i] = $tauQ * $zmat[$i] - $zkQ * $vlag[$i];
        }

        return Util::INFO_DFT;
    }

    /**
     * Givens rotations zeroing ZMAT(KNEW, 2:NZ), shared with rescue's updateh_rsc.
     *
     * @param  list<float>  $zmat
     */
    public static function rotateZmatRow(array &$zmat, int $knew, int $npt, int $nz): void
    {
        for ($j = 1; $j < $nz; $j++) {
            $maxAbs = 0.0;
            foreach ($zmat as $v) {
                $maxAbs = Linalg::max($maxAbs, abs($v));
            }
            $joff = $j * $npt;
            if (abs($zmat[$knew + $joff]) > Util::F32_1E_20 * $maxAbs) {
                [$c, $s] = Linalg::planerot($zmat[$knew], $zmat[$knew + $joff]);
                for ($i = 0; $i < $npt; $i++) {
                    $a = $zmat[$i];
                    $b = $zmat[$joff + $i];
                    $zmat[$i] = $a * $c + $b * $s;
                    $zmat[$joff + $i] = $a * -$s + $b * $c;
                }
            }
            $zmat[$knew + $joff] = 0.0;
        }
    }

    /**
     * Update XPT, FVAL and KOPT when XPT(:, KNEW) becomes XNEW.
     *
     * @param  list<float>  $xnew
     * @param  list<float>  $fval
     * @param  list<float>  $xpt
     */
    public static function updatexf(?int $knew, bool $ximproved, float $f, array $xnew, int &$kopt, array &$fval, array &$xpt, int $n): void
    {
        if ($knew === null) {
            return;
        }
        for ($i = 0; $i < $n; $i++) {
            $xpt[$i + $knew * $n] = $xnew[$i];
        }
        $fval[$knew] = $f;
        if ($ximproved) {
            $kopt = $knew;
        }
    }

    /**
     * Update GOPT, HQ and PQ when XPT(:, KNEW) changes.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $d
     * @param  list<float>  $xdrop
     * @param  list<float>  $xosav
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @param  list<float>  $gopt
     * @param  list<float>  $hq
     * @param  list<float>  $pq
     */
    public static function updateq(
        ?int $knew,
        bool $ximproved,
        array $bmat,
        array $d,
        float $moderr,
        array $xdrop,
        array $xosav,
        array $xpt,
        array $zmat,
        array &$gopt,
        array &$hq,
        array &$pq,
        int $n,
        int $npt,
    ): void {
        if ($knew === null) {
            return;
        }
        $nz = $npt - $n - 1;

        Linalg::r1update($hq, $n, $pq[$knew], $xdrop);
        $pq[$knew] = 0.0;

        $zrow = [];
        for ($j = 0; $j < $nz; $j++) {
            $zrow[] = $zmat[$knew + $j * $npt];
        }
        $zmatZrow = Linalg::matprod21($zmat, $npt, $nz, $zrow);
        $pqinc = [];
        for ($k = 0; $k < $npt; $k++) {
            $pqinc[] = $moderr * $zmatZrow[$k];
        }
        for ($k = 0; $k < $npt; $k++) {
            $pq[$k] += $pqinc[$k];
        }

        $hm = PowAlg::hessMul($xosav, $xpt, $n, $npt, $pqinc, null);
        for ($i = 0; $i < $n; $i++) {
            $gopt[$i] += $moderr * $bmat[$i + $knew * $n];
        }
        for ($i = 0; $i < $n; $i++) {
            $gopt[$i] += $hm[$i];
        }

        if ($ximproved) {
            $hm = PowAlg::hessMul($d, $xpt, $n, $npt, $pq, $hq);
            for ($i = 0; $i < $n; $i++) {
                $gopt[$i] += $hm[$i];
            }
        }
    }

    /**
     * Test whether to replace Q with the least Frobenius norm interpolant.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $fval
     * @param  list<float>  $sl
     * @param  list<float>  $su
     * @param  list<float>  $xopt
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @param  list<float>  $gopt
     * @param  list<float>  $hq
     * @param  list<float>  $pq
     */
    public static function tryqalt(
        array $bmat,
        array $fval,
        float $ratio,
        array $sl,
        array $su,
        array $xopt,
        array $xpt,
        array $zmat,
        int &$itest,
        array &$gopt,
        array &$hq,
        array &$pq,
        int $n,
        int $npt,
    ): void {
        if ($ratio > 0.1) {
            $itest = 0;

            return;
        }
        $nz = $npt - $n - 1;

        $pgopt = $gopt;
        for ($i = 0; $i < $n; $i++) {
            if ($xopt[$i] >= $su[$i]) {
                $pgopt[$i] = Linalg::max(0.0, $gopt[$i]);
            }
        }
        for ($i = 0; $i < $n; $i++) {
            if ($xopt[$i] <= $sl[$i]) {
                $pgopt[$i] = Linalg::min(0.0, $gopt[$i]);
            }
        }

        $inner = Linalg::matprod12($fval, $zmat, $npt, $nz);
        $pqalt = Linalg::matprod21($zmat, $npt, $nz, $inner);

        $galt = array_fill(0, $n, 0.0);
        for ($j = 0; $j < $npt; $j++) {
            $fvalj = $fval[$j];
            $offset = $j * $n;
            for ($i = 0; $i < $n; $i++) {
                $galt[$i] += $bmat[$offset + $i] * $fvalj;
            }
        }
        $hm = PowAlg::hessMul($xopt, $xpt, $n, $npt, $pqalt, null);
        for ($i = 0; $i < $n; $i++) {
            $galt[$i] += $hm[$i];
        }

        $pgalt = $galt;
        for ($i = 0; $i < $n; $i++) {
            if ($xopt[$i] >= $su[$i]) {
                $pgalt[$i] = Linalg::max(0.0, $galt[$i]);
            }
        }
        for ($i = 0; $i < $n; $i++) {
            if ($xopt[$i] <= $sl[$i]) {
                $pgalt[$i] = Linalg::min(0.0, $galt[$i]);
            }
        }

        if (Linalg::inprod($pgopt, $pgopt) < 10.0 * Linalg::inprod($pgalt, $pgalt)) {
            $itest = 0;
        } else {
            $itest++;
        }

        if ($itest >= 3) {
            $gopt = $galt;
            $pq = $pqalt;
            $hq = array_fill(0, $n * $n, 0.0);
            $itest = 0;
        }
    }
}
