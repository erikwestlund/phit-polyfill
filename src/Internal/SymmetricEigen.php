<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/**
 * nalgebra 0.35.0's SymmetricEigen::new for a dynamically sized f64 matrix:
 * scale by the largest magnitude, Householder tridiagonalization reading only
 * the lower triangle, then implicit QR with a Wilkinson shift and
 * eps = f64::EPSILON, with no iteration limit. Eigenvalues come back in
 * nalgebra's (unsorted) order; eigenvector j is column j.
 *
 * Matrices are column-major: $eigenvectors[$j][$i] is component i of vector j.
 *
 * @internal
 */
final class SymmetricEigen
{
    /** @var list<float> */
    public array $eigenvalues;

    /** @var list<list<float>> */
    public array $eigenvectors;

    /**
     * @param  list<list<float>>  $columns  a square matrix as columns
     * @param  int  $maxIterations  0 iterates until convergence, as SymmetricEigen::new does
     */
    public function __construct(array $columns, float $eps = PHP_FLOAT_EPSILON, int $maxIterations = 0)
    {
        $dim = count($columns);
        if ($dim === 0 || count($columns[0]) !== $dim) {
            throw new \InvalidArgumentException('Unable to compute the eigendecomposition of a non-square or empty matrix.');
        }
        $amax = Dense::camax($columns);
        if ($amax != 0.0) {
            foreach ($columns as $j => $column) {
                foreach ($column as $i => $value) {
                    $columns[$j][$i] = \fdiv($value, $amax);
                }
            }
        }

        [$tri, $offRaw] = self::tridiagonalize($columns, $dim);
        $q = self::assembleQ($tri, $offRaw, $dim);
        $diag = [];
        for ($i = 0; $i < $dim; $i++) {
            $diag[$i] = $tri[$i][$i];
        }
        $off = \array_map(\abs(...), $offRaw);

        if ($dim === 1) {
            $diag[0] *= $amax;
            $this->eigenvalues = $diag;
            $this->eigenvectors = $q;

            return;
        }

        $niter = 0;
        [$start, $end] = self::delimit($diag, $off, $dim - 1, $eps);
        while ($end !== $start) {
            $subdim = $end - $start + 1;
            if ($subdim > 2) {
                $m = $end - 1;
                $n = $end;
                $vx = $diag[$start] - self::wilkinsonShift($diag[$m], $diag[$n], $off[$m]);
                $vy = $off[$start];
                for ($i = $start; $i < $n; $i++) {
                    $j = $i + 1;
                    $rot = Svd::cancelY($vx, $vy);
                    if ($rot === null) {
                        break;
                    }
                    [$c, $s, $norm] = $rot;
                    if ($i > $start) {
                        $off[$i - 1] = $norm;
                    }
                    $mii = $diag[$i];
                    $mjj = $diag[$j];
                    $mij = $off[$i];
                    $cc = $c * $c;
                    $ss = $s * $s;
                    $cs = $c * $s;
                    $b = $cs * 2.0 * $mij;
                    $diag[$i] = ($cc * $mii + $ss * $mjj) - $b;
                    $diag[$j] = ($ss * $mii + $cc * $mjj) + $b;
                    $off[$i] = $cs * ($mii - $mjj) + $mij * ($cc - $ss);
                    if ($i !== $n - 1) {
                        $vx = $off[$i];
                        $vy = (-$s) * $off[$i + 1];
                        $off[$i + 1] *= $c;
                    }
                    // rot.inverse().rotate_rows on columns i and j of Q.
                    $x = $q[$i];
                    $y = $q[$j];
                    foreach ($x as $r => $xr) {
                        $yr = $y[$r];
                        $x[$r] = $xr * $c + (-$s) * $yr;
                        $y[$r] = $s * $xr + $yr * $c;
                    }
                    $q[$i] = $x;
                    $q[$j] = $y;
                }
                if (\abs($off[$m]) <= $eps * (\abs($diag[$m]) + \abs($diag[$n]))) {
                    $end--;
                }
            } elseif ($subdim === 2) {
                $d0 = $diag[$start];
                $d1 = $diag[$start + 1];
                $e = $off[$start];
                // Matrix2::eigenvalues (compute_2x2_eigvals) on [[d0, e], [e, d1]].
                $val = ($d0 - $d1) * 0.5;
                $discriminant = $e * $e + $val * $val;
                if ($discriminant >= 0.0) {
                    $root = \sqrt($discriminant);
                } else {
                    throw new \RuntimeException('Symmetric eigen: complex 2x2 eigenvalues.');
                }
                $halfTrace = ($d0 + $d1) * 0.5;
                $lambda0 = $halfTrace + $root;
                $lambda1 = $halfTrace - $root;
                if (\abs($lambda0 - $d1) > \abs($lambda0 - $d0)) {
                    $bx = $lambda0 - $d1;
                    $by = $e;
                } else {
                    $bx = $e;
                    $by = $lambda0 - $d0;
                }
                $diag[$start] = $lambda0;
                $diag[$start + 1] = $lambda1;
                $rot = Svd::tryGivens($bx, $by, $eps);
                if ($rot !== null) {
                    [$c, $s] = $rot;
                    $x = $q[$start];
                    $y = $q[$start + 1];
                    foreach ($x as $r => $xr) {
                        $yr = $y[$r];
                        $x[$r] = $xr * $c + $s * $yr;
                        $y[$r] = (-$s) * $xr + $yr * $c;
                    }
                    $q[$start] = $x;
                    $q[$start + 1] = $y;
                }
                $end--;
            }
            [$start, $end] = self::delimit($diag, $off, $end, $eps);
            $niter++;
            if ($niter === $maxIterations) {
                throw new \RuntimeException('Symmetric eigen did not converge.');
            }
        }

        for ($i = 0; $i < $dim; $i++) {
            $diag[$i] *= $amax;
        }
        $this->eigenvalues = $diag;
        $this->eigenvectors = $q;
    }

    /** symmetric_eigen::wilkinson_shift. */
    public static function wilkinsonShift(float $tmm, float $tnn, float $tmn): float
    {
        $sq = $tmn * $tmn;
        if ($sq == 0.0) {
            return $tnn;
        }
        $d = ($tmm - $tnn) * 0.5;

        return $tnn - \fdiv($sq, $d + Dense::signum($d) * \sqrt($d * $d + $sq));
    }

    /**
     * SymmetricTridiagonal::new: returns the matrix holding the reflectors below
     * the subdiagonal and the tridiagonal's diagonal, and the signed
     * off-diagonal.
     *
     * @return array{list<list<float>>, list<float>}
     */
    private static function tridiagonalize(array $m, int $dim): array
    {
        $off = [];
        $p = \array_fill(0, \max(0, $dim - 1), 0.0);
        for ($i = 0; $i < $dim - 1; $i++) {
            $first = $i + 1;
            $axis = \array_slice($m[$i], $first);
            [$norm, $notZero] = Dense::reflectionAxis($axis);
            foreach ($axis as $k => $value) {
                $m[$i][$first + $k] = $value;
            }
            $off[$i] = $norm;
            if ($notZero) {
                $length = $dim - $first;
                // Local views: a(r, c) = m[first + c][first + r], x = axis, y = p[i..].
                // hegemv(2, a, x, 0), reading the lower triangle.
                $y = [];
                $col = $m[$first];
                $coefficient = 2.0 * $axis[0];
                for ($r = 0; $r < $length; $r++) {
                    $y[$r] = $coefficient * $col[$first + $r];
                }
                $y[0] += 2.0 * self::dotOffset($col, $first + 1, $axis, 1, $length - 1);
                for ($j = 1; $j < $length; $j++) {
                    $col = $m[$first + $j];
                    $dot = self::dotOffset($col, $first + $j, $axis, $j, $length - $j);
                    $y[$j] += 2.0 * $dot;
                    $coefficient = 2.0 * $axis[$j];
                    for ($r = $j + 1; $r < $length; $r++) {
                        $y[$r] = $coefficient * $col[$first + $r] + $y[$r];
                    }
                }
                for ($r = 0; $r < $length; $r++) {
                    $p[$i + $r] = $y[$r];
                }

                $dot = Dense::dot($axis, $y);
                // hegerc(-1, p, axis, 1), hegerc(-1, axis, p, 1), hegerc(2 dot, axis, axis, 1)
                for ($j = 0; $j < $length; $j++) {
                    $col = $m[$first + $j];
                    $coefficient = -$axis[$j];
                    for ($r = $j; $r < $length; $r++) {
                        $col[$first + $r] = $coefficient * $y[$r] + $col[$first + $r];
                    }
                    $coefficient = -$y[$j];
                    for ($r = $j; $r < $length; $r++) {
                        $col[$first + $r] = $coefficient * $axis[$r] + $col[$first + $r];
                    }
                    $m[$first + $j] = $col;
                }
                $twice = $dot * 2.0;
                for ($j = 0; $j < $length; $j++) {
                    $col = $m[$first + $j];
                    $coefficient = $twice * $axis[$j];
                    for ($r = $j; $r < $length; $r++) {
                        $col[$first + $r] = $coefficient * $axis[$r] + $col[$first + $r];
                    }
                    $m[$first + $j] = $col;
                }
            }
        }

        return [$m, $off];
    }

    /** dotx over $length entries of $x from $xStart and $y from $yStart. */
    private static function dotOffset(array $x, int $xStart, array $y, int $yStart, int $length): float
    {
        $acc0 = $acc1 = $acc2 = $acc3 = $acc4 = $acc5 = $acc6 = $acc7 = 0.0;
        $k = 0;
        while ($length - $k >= 8) {
            $a = $xStart + $k;
            $b = $yStart + $k;
            $acc0 += $x[$a] * $y[$b];
            $acc1 += $x[$a + 1] * $y[$b + 1];
            $acc2 += $x[$a + 2] * $y[$b + 2];
            $acc3 += $x[$a + 3] * $y[$b + 3];
            $acc4 += $x[$a + 4] * $y[$b + 4];
            $acc5 += $x[$a + 5] * $y[$b + 5];
            $acc6 += $x[$a + 6] * $y[$b + 6];
            $acc7 += $x[$a + 7] * $y[$b + 7];
            $k += 8;
        }
        $result = 0.0;
        $result += $acc0 + $acc4;
        $result += $acc1 + $acc5;
        $result += $acc2 + $acc6;
        $result += $acc3 + $acc7;
        for (; $k < $length; $k++) {
            $result += $x[$xStart + $k] * $y[$yStart + $k];
        }

        return $result;
    }

    /** householder::assemble_q: no zero-axis check, unlike the SVD's U. */
    private static function assembleQ(array $tri, array $off, int $dim): array
    {
        $res = [];
        for ($j = 0; $j < $dim; $j++) {
            $column = \array_fill(0, $dim, 0.0);
            $column[$j] = 1.0;
            $res[$j] = $column;
        }
        for ($i = $dim - 2; $i >= 0; $i--) {
            $first = $i + 1;
            $axis = \array_slice($tri[$i], $first);
            $sign = Dense::signum($off[$i]);
            $scale = $sign * -2.0;
            $length = count($axis);
            for ($j = $i; $j < $dim; $j++) {
                $column = $res[$j];
                $factor = self::dotOffset($axis, 0, $column, $first, $length) * $scale;
                for ($k = 0; $k < $length; $k++) {
                    $r = $first + $k;
                    $column[$r] = $factor * $axis[$k] + $sign * $column[$r];
                }
                $res[$j] = $column;
            }
        }

        return $res;
    }

    /**
     * @param  list<float>  $diag
     * @param  list<float>  $off
     * @return array{int, int}
     */
    private static function delimit(array $diag, array &$off, int $end, float $eps): array
    {
        $n = $end;
        while ($n > 0) {
            $m = $n - 1;
            if (\abs($off[$m]) > $eps * (\abs($diag[$n]) + \abs($diag[$m]))) {
                break;
            }
            $n--;
        }
        if ($n === 0) {
            return [0, 0];
        }
        $newStart = $n - 1;
        while ($newStart > 0) {
            $m = $newStart - 1;
            if ($off[$m] == 0.0 || \abs($off[$m]) <= $eps * (\abs($diag[$newStart]) + \abs($diag[$m]))) {
                $off[$m] = 0.0;
                break;
            }
            $newStart--;
        }

        return [$newStart, $n];
    }
}
