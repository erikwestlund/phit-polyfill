<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal\Bobyqa;

/**
 * The crate's powalg.rs: PRIMA's powalg.f90 subset (idz = 1). XPT is n x npt,
 * BMAT n x (npt + n), ZMAT npt x (npt - n - 1), HQ n x n, all flat column-major.
 *
 * @internal
 */
final class PowAlg
{
    /**
     * The 0-based (P, Q) pairs for the interpolation points beyond the first 2n + 1.
     *
     * @return list<array{int, int}>
     */
    public static function setij(int $n, int $npt): array
    {
        $ij = [];
        if ($npt > 2 * $n + 1) {
            for ($k = $n; $k <= $npt - $n - 2; $k++) {
                $ell = intdiv($k, $n);
                $i1 = $k - $n * $ell + 1;
                $i2 = ($i1 + $ell - 1) % $n + 1;
                $ij[] = [$i1 - 1, $i2 - 1];
            }
        }

        return $ij;
    }

    /**
     * HESSIAN * x with HESSIAN = HQ (zero when null) + sum_k PQ(k) XPT(:, k) XPT(:, k)^T.
     *
     * @param  list<float>  $x
     * @param  list<float>  $xpt
     * @param  list<float>  $pq
     * @param  list<float>|null  $hq
     * @return list<float>
     */
    public static function hessMul(array $x, array $xpt, int $n, int $npt, array $pq, ?array $hq): array
    {
        $dxpt = Linalg::matprod12($x, $xpt, $n, $npt);
        for ($k = 0; $k < $npt; $k++) {
            $dxpt[$k] *= $pq[$k];
        }
        $y = Linalg::matprod21($xpt, $n, $npt, $dxpt);
        if ($hq !== null) {
            for ($j = 0; $j < $n; $j++) {
                $xj = $x[$j];
                $offset = $j * $n;
                for ($i = 0; $i < $n; $i++) {
                    $y[$i] += $hq[$offset + $i] * $xj;
                }
            }
        }

        return $y;
    }

    /**
     * Q(d) - Q(0) for the model [gq, hq, pq].
     *
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $gq
     * @param  list<float>  $pq
     * @param  list<float>  $hq
     */
    public static function quadinc(array $d, array $xpt, int $n, int $npt, array $gq, array $pq, array $hq): float
    {
        $dxpt = Linalg::matprod12($d, $xpt, $n, $npt);
        $hqd = Linalg::matprod21($hq, $n, $n, $d);
        for ($i = 0; $i < $n; $i++) {
            $hqd[$i] = $gq[$i] + 0.5 * $hqd[$i];
        }
        $pqdxpt = [];
        for ($k = 0; $k < $npt; $k++) {
            $pqdxpt[] = $pq[$k] * $dxpt[$k];
        }

        return Linalg::inprod($d, $hqd) + 0.5 * Linalg::inprod($dxpt, $pqdxpt);
    }

    /**
     * VLAG = H * w for the step d from XPT(:, kref), without the +1 at kref; also
     * returns WCHECK through the reference.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @param  list<float>|null  $wcheck
     * @return list<float>
     */
    public static function calvlagNoadd(int $kref, array $bmat, array $d, array $xpt, array $zmat, int $n, int $npt, ?array &$wcheck): array
    {
        $nz = $npt - $n - 1;
        $xref = Linalg::col($xpt, $n, $kref);
        $wcheck = Linalg::matprod12($d, $xpt, $n, $npt);
        $xrefxpt = Linalg::matprod12($xref, $xpt, $n, $npt);
        for ($k = 0; $k < $npt; $k++) {
            $wcheck[$k] *= 0.5 * $wcheck[$k] + $xrefxpt[$k];
        }
        $xz = Linalg::matprod12($wcheck, $zmat, $npt, $nz);
        $omega = Linalg::matprod21($zmat, $npt, $nz, $xz);
        $vlag = [];
        for ($j = 0; $j < $npt; $j++) {
            $s = 0.0;
            $offset = $j * $n;
            for ($i = 0; $i < $n; $i++) {
                $s += $d[$i] * $bmat[$offset + $i];
            }
            $vlag[] = $omega[$j] + $s;
        }
        $wmv = array_merge($wcheck, $d);
        foreach (Linalg::matprod21($bmat, $n, $npt + $n, $wmv) as $v) {
            $vlag[] = $v;
        }

        return $vlag;
    }

    /**
     * @param  list<float>  $d
     * @param  list<float>  $xref
     * @param  list<float>  $wcheck
     * @param  list<float>  $vlag
     */
    private static function betaFromNoaddVlag(array $d, array $xref, array $wcheck, array $vlag, int $n, int $npt): float
    {
        $dxref = Linalg::inprod($d, $xref);
        $dsq = Linalg::inprod($d, $d);
        $xrefsq = Linalg::inprod($xref, $xref);
        $dvlag = Linalg::inprod($d, array_slice($vlag, $npt, $n));
        $wvlag = Linalg::inprod($wcheck, $vlag);

        return $dxref * $dxref + $dsq * ($xrefsq + $dxref + $dxref + 0.5 * $dsq) - $dvlag - $wvlag;
    }

    /**
     * PRIMA calvlag_lfqint: VLAG = H * w (length npt + n).
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @return list<float>
     */
    public static function calvlag(int $kref, array $bmat, array $d, array $xpt, array $zmat, int $n, int $npt): array
    {
        $wcheck = null;
        $vlag = self::calvlagNoadd($kref, $bmat, $d, $xpt, $zmat, $n, $npt, $wcheck);
        $vlag[$kref] += 1.0;

        return $vlag;
    }

    /**
     * PRIMA calbeta.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     */
    public static function calbeta(int $kref, array $bmat, array $d, array $xpt, array $zmat, int $n, int $npt): float
    {
        $wcheck = null;
        $vlag = self::calvlagNoadd($kref, $bmat, $d, $xpt, $zmat, $n, $npt, $wcheck);

        return self::betaFromNoaddVlag($d, Linalg::col($xpt, $n, $kref), $wcheck, $vlag, $n, $npt);
    }

    /**
     * PRIMA calden: DEN(k) if XPT(:, k) were replaced with XPT(:, kref) + d.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @return list<float>
     */
    public static function calden(int $kref, array $bmat, array $d, array $xpt, array $zmat, int $n, int $npt): array
    {
        $vlag = [];
        $den = [];
        self::calvlagAndDen($kref, $bmat, $d, $xpt, $zmat, $n, $npt, $vlag, $den);

        return $den;
    }

    /**
     * VLAG (with the +1 at kref) and DEN from one kernel; returns BETA.
     *
     * @param  list<float>  $bmat
     * @param  list<float>  $d
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @param  list<float>  $vlag
     * @param  list<float>  $den
     */
    public static function calvlagAndDen(int $kref, array $bmat, array $d, array $xpt, array $zmat, int $n, int $npt, array &$vlag, array &$den): float
    {
        $nz = $npt - $n - 1;
        $hdiag = array_fill(0, $npt, 0.0);
        for ($j = 0; $j < $nz; $j++) {
            $offset = $j * $npt;
            for ($k = 0; $k < $npt; $k++) {
                $z = $zmat[$offset + $k];
                $hdiag[$k] += $z * $z;
            }
        }
        $wcheck = null;
        $vlag = self::calvlagNoadd($kref, $bmat, $d, $xpt, $zmat, $n, $npt, $wcheck);
        $beta = self::betaFromNoaddVlag($d, Linalg::col($xpt, $n, $kref), $wcheck, $vlag, $n, $npt);
        $vlag[$kref] += 1.0;
        $den = [];
        for ($k = 0; $k < $npt; $k++) {
            $den[] = $hdiag[$k] * $beta + $vlag[$k] * $vlag[$k];
        }

        return $beta;
    }
}
