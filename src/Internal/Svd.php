<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/**
 * nalgebra 0.35.0's SVD::new for a dynamically sized f64 matrix (DMatrix, so
 * never the special 2x2/3x3 paths): scale by the largest magnitude, Householder
 * bidiagonalization (upper when rows >= columns, lower otherwise; no
 * transposition), implicit-shift QR with eps = 5 * f64::EPSILON and no
 * iteration limit, sign fixing that flips columns of U, then a descending sort
 * of the singular values that permutes U's columns and V^T's rows.
 *
 * Matrices are column-major: $u[$j][$i] = U(i, j), $vT[$j][$i] = V^T(i, j).
 *
 * @internal
 */
final class Svd
{
    /** @var list<float> */
    public array $singularValues;

    /** @var list<list<float>>|null rows x min(rows, cols), as columns */
    public ?array $u;

    /** @var list<list<float>>|null min(rows, cols) x cols, as columns */
    public ?array $vT;

    public readonly int $rows;

    public readonly int $cols;

    /** @var list<float> */
    private array $d;

    /** @var list<float> */
    private array $f;

    private bool $upper;

    /**
     * @param  list<list<float>>  $columns  rows x cols, as columns
     * @param  int  $maxIterations  0 iterates until convergence, as SVD::new does; otherwise
     *                              nalgebra's try_new limit (a RuntimeException when reached)
     */
    public function __construct(
        array $columns,
        bool $computeU,
        bool $computeV,
        bool $sorted = true,
        float $eps = 5.0 * PHP_FLOAT_EPSILON,
        int $maxIterations = 0,
    ) {
        $cols = count($columns);
        $rows = $cols > 0 ? count($columns[0]) : 0;
        if ($rows === 0 || $cols === 0) {
            throw new \InvalidArgumentException('Cannot compute the SVD of an empty matrix.');
        }
        $this->rows = $rows;
        $this->cols = $cols;
        $dim = min($rows, $cols);

        $amax = Dense::camax($columns);
        if ($amax != 0.0) {
            foreach ($columns as $j => $column) {
                foreach ($column as $i => $value) {
                    $columns[$j][$i] = fdiv($value, $amax);
                }
            }
        }

        [$uv, $diagonal, $offDiagonal, $upper] = self::bidiagonalize($columns, $rows, $cols, $dim);
        $this->upper = $upper;
        $this->u = $computeU ? self::assembleU($uv, $diagonal, $offDiagonal, $upper, $rows, $dim) : null;
        $this->vT = $computeV ? self::assembleVt($uv, $diagonal, $offDiagonal, $upper, $cols, $dim) : null;
        $this->d = array_map(abs(...), $diagonal);
        $this->f = array_map(abs(...), $offDiagonal);

        $niter = 0;
        [$start, $end] = $this->delimit($dim - 1, $eps);
        while ($end !== $start) {
            $subdim = $end - $start + 1;
            if ($subdim > 2) {
                $this->qrStep($start, $end);
            } elseif ($subdim === 2) {
                $this->solve2x2($start, $computeU, $computeV);
                $end--;
            }
            [$start, $end] = $this->delimit($end, $eps);
            $niter++;
            if ($niter === $maxIterations) {
                throw new \RuntimeException('SVD did not converge.');
            }
        }

        $d = $this->d;
        for ($i = 0; $i < $dim; $i++) {
            $d[$i] *= $amax;
        }
        for ($i = 0; $i < $dim; $i++) {
            if ($d[$i] < 0.0) {
                $d[$i] = -$d[$i];
                if ($this->u !== null) {
                    foreach ($this->u[$i] as $r => $value) {
                        $this->u[$i][$r] = -$value;
                    }
                }
            }
        }
        $this->singularValues = $d;

        if ($sorted) {
            $this->sortBySingularValues();
        }
    }

    /**
     * nalgebra's SVD::solve: V (S^+ (U^T b)) with singular values at or below
     * $eps treated as zero.
     *
     * @param  list<list<float>>  $b  rows x k, as columns
     * @return list<list<float>> cols x k, as columns
     */
    public function solve(array $b, float $eps): array
    {
        if ($eps < 0.0) {
            throw new \InvalidArgumentException('SVD solve: the epsilon must be non-negative.');
        }
        if ($this->u === null || $this->vT === null) {
            throw new \LogicException('SVD solve: U and V^t have not been computed.');
        }
        $utb = Dense::transposeMultiply($this->u, $b);
        foreach ($utb as $j => $column) {
            foreach ($this->singularValues as $i => $value) {
                $utb[$j][$i] = $value > $eps ? fdiv($column[$i], $value) : 0.0;
            }
        }

        return Dense::transposeMultiply($this->vT, $utb);
    }

    /**
     * @param  list<float>  $b
     * @return list<float>
     */
    public function solveVector(array $b, float $eps): array
    {
        return $this->solve([$b], $eps)[0];
    }

    /**
     * Bidiagonal::new. Returns the packed reflectors, the signed diagonal and
     * off-diagonal, and whether the bidiagonal is upper.
     *
     * @return array{list<list<float>>, list<float>, list<float>, bool}
     */
    private static function bidiagonalize(array $m, int $rows, int $cols, int $dim): array
    {
        $diagonal = array_fill(0, $dim, 0.0);
        $offDiagonal = array_fill(0, max(0, $dim - 1), 0.0);
        $upper = $rows >= $cols;
        if ($upper) {
            for ($i = 0; $i < $dim - 1; $i++) {
                $diagonal[$i] = self::clearColumn($m, $i, 0, $rows, $cols);
                $offDiagonal[$i] = self::clearRow($m, $i, 1, $rows, $cols);
            }
            $diagonal[$dim - 1] = self::clearColumn($m, $dim - 1, 0, $rows, $cols);
        } else {
            for ($i = 0; $i < $dim - 1; $i++) {
                $diagonal[$i] = self::clearRow($m, $i, 0, $rows, $cols);
                $offDiagonal[$i] = self::clearColumn($m, $i, 1, $rows, $cols);
            }
            $diagonal[$dim - 1] = self::clearRow($m, $dim - 1, 0, $rows, $cols);
        }

        return [$m, $diagonal, $offDiagonal, $upper];
    }

    /** householder::clear_column_unchecked without the bilateral work vector. */
    private static function clearColumn(array &$m, int $icol, int $shift, int $rows, int $cols): float
    {
        $first = $icol + $shift;
        $axis = array_slice($m[$icol], $first);
        [$norm, $notZero] = Dense::reflectionAxis($axis);
        foreach ($axis as $r => $value) {
            $m[$icol][$first + $r] = $value;
        }
        if ($notZero) {
            $sign = Dense::signum($norm);
            self::reflectColumns($m, $axis, $first, $icol + 1, $cols, $sign);
        }

        return $norm;
    }

    /** householder::clear_row_unchecked. */
    private static function clearRow(array &$m, int $irow, int $shift, int $rows, int $cols): float
    {
        $first = $irow + $shift;
        $axis = [];
        for ($c = $first; $c < $cols; $c++) {
            $axis[] = $m[$c][$irow];
        }
        [$norm, $notZero] = Dense::reflectionAxis($axis);
        if ($notZero) {
            self::reflectRows($m, $axis, $irow + 1, $rows, $first, Dense::signum($norm));
        }
        foreach ($axis as $k => $value) {
            $m[$first + $k][$irow] = $value;
        }

        return $norm;
    }

    /**
     * Reflection::reflect_with_sign on columns [$fromCol, $toCol) of $m, rows
     * from $firstRow on, with the axis covering those rows.
     */
    private static function reflectColumns(array &$m, array $axis, int $firstRow, int $fromCol, int $toCol, float $sign): void
    {
        $length = count($axis);
        $scale = $sign * -2.0;
        for ($j = $fromCol; $j < $toCol; $j++) {
            $column = $m[$j];
            $acc0 = $acc1 = $acc2 = $acc3 = $acc4 = $acc5 = $acc6 = $acc7 = 0.0;
            $k = 0;
            while ($length - $k >= 8) {
                $r = $firstRow + $k;
                $acc0 += $axis[$k] * $column[$r];
                $acc1 += $axis[$k + 1] * $column[$r + 1];
                $acc2 += $axis[$k + 2] * $column[$r + 2];
                $acc3 += $axis[$k + 3] * $column[$r + 3];
                $acc4 += $axis[$k + 4] * $column[$r + 4];
                $acc5 += $axis[$k + 5] * $column[$r + 5];
                $acc6 += $axis[$k + 6] * $column[$r + 6];
                $acc7 += $axis[$k + 7] * $column[$r + 7];
                $k += 8;
            }
            $dot = 0.0;
            $dot += $acc0 + $acc4;
            $dot += $acc1 + $acc5;
            $dot += $acc2 + $acc6;
            $dot += $acc3 + $acc7;
            for (; $k < $length; $k++) {
                $dot += $axis[$k] * $column[$firstRow + $k];
            }
            $factor = $dot * $scale;
            for ($k = 0; $k < $length; $k++) {
                $r = $firstRow + $k;
                $column[$r] = $factor * $axis[$k] + $sign * $column[$r];
            }
            $m[$j] = $column;
        }
    }

    /**
     * Reflection::reflect_rows_with_sign on rows [$fromRow, $toRow) of $m,
     * columns from $firstCol on, with the axis covering those columns: the
     * work vector is the gemv product of that block with the axis, then a
     * rank-one update column by column.
     */
    private static function reflectRows(array &$m, array $axis, int $fromRow, int $toRow, int $firstCol, float $sign): void
    {
        $length = count($axis);
        if ($toRow <= $fromRow) {
            return;
        }
        $work = [];
        $column = $m[$firstCol];
        $value = $axis[0];
        for ($r = $fromRow; $r < $toRow; $r++) {
            $work[$r] = $column[$r] * $value;
        }
        for ($k = 1; $k < $length; $k++) {
            $column = $m[$firstCol + $k];
            $value = $axis[$k];
            for ($r = $fromRow; $r < $toRow; $r++) {
                $work[$r] = $column[$r] * $value + $work[$r];
            }
        }
        $scale = $sign * -2.0;
        for ($k = 0; $k < $length; $k++) {
            $c = $firstCol + $k;
            $column = $m[$c];
            $coefficient = $scale * $axis[$k];
            for ($r = $fromRow; $r < $toRow; $r++) {
                $column[$r] = $coefficient * $work[$r] + $sign * $column[$r];
            }
            $m[$c] = $column;
        }
    }

    /** Bidiagonal::u. */
    private static function assembleU(array $uv, array $diagonal, array $offDiagonal, bool $upper, int $rows, int $dim): array
    {
        $res = [];
        for ($j = 0; $j < $dim; $j++) {
            $column = array_fill(0, $rows, 0.0);
            $column[$j] = 1.0;
            $res[$j] = $column;
        }
        $shift = $upper ? 0 : 1;
        for ($i = $dim - $shift - 1; $i >= 0; $i--) {
            $first = $i + $shift;
            $axis = array_slice($uv[$i], $first);
            if (0.0 + Dense::dot($axis, $axis) == 0.0) {
                continue;
            }
            $sign = Dense::signum($upper ? $diagonal[$i] : $offDiagonal[$i]);
            self::reflectColumns($res, $axis, $first, $i, $dim, $sign);
        }

        return $res;
    }

    /** Bidiagonal::v_t. */
    private static function assembleVt(array $uv, array $diagonal, array $offDiagonal, bool $upper, int $cols, int $dim): array
    {
        $res = [];
        for ($j = 0; $j < $cols; $j++) {
            $column = array_fill(0, $dim, 0.0);
            if ($j < $dim) {
                $column[$j] = 1.0;
            }
            $res[$j] = $column;
        }
        $shift = $upper ? 1 : 0;
        for ($i = $dim - $shift - 1; $i >= 0; $i--) {
            $first = $i + $shift;
            $axis = [];
            for ($c = $first; $c < $cols; $c++) {
                $axis[] = $uv[$c][$i];
            }
            if (0.0 + Dense::dot($axis, $axis) == 0.0) {
                continue;
            }
            $sign = Dense::signum($upper ? $offDiagonal[$i] : $diagonal[$i]);
            self::reflectRows($res, $axis, $i, $dim, $first, $sign);
        }

        return $res;
    }

    /** One implicit-shift QR sweep over the subproblem [$start, $end] (subdim > 2). */
    private function qrStep(int $start, int $end): void
    {
        $d = &$this->d;
        $f = &$this->f;
        $m = $end - 1;
        $n = $end;
        $dm = $d[$m];
        $dn = $d[$n];
        $fm = $f[$m];
        $tmm = $dm * $dm + $f[$m - 1] * $f[$m - 1];
        $tmn = $dm * $fm;
        $tnn = $dn * $dn + $fm * $fm;
        $shift = SymmetricEigen::wilkinsonShift($tmm, $tnn, $tmn);
        $vx = $d[$start] * $d[$start] - $shift;
        $vy = $d[$start] * $f[$start];

        for ($k = $start; $k < $n; $k++) {
            $m12 = $k === $n - 1 ? 0.0 : $f[$k + 1];
            // subm = [[s00, s01, s02], [s10, s11, s12]]
            $s00 = $d[$k];
            $s01 = $f[$k];
            $s02 = 0.0;
            $s10 = 0.0;
            $s11 = $d[$k + 1];
            $s12 = $m12;

            $rot1 = self::cancelY($vx, $vy);
            if ($rot1 === null) {
                break;
            }
            [$c1, $sn1, $norm1] = $rot1;
            // rot1.inverse().rotate_rows on the first two columns.
            [$s00, $s01] = [$s00 * $c1 + (-$sn1) * $s01, $sn1 * $s00 + $s01 * $c1];
            [$s10, $s11] = [$s10 * $c1 + (-$sn1) * $s11, $sn1 * $s10 + $s11 * $c1];

            if ($k > $start) {
                $f[$k - 1] = $norm1;
            }

            $rot2 = self::cancelY($s00, $s10) ?? [1.0, 0.0, $s00];
            [$c2, $sn2, $norm2] = $rot2;
            // rot2.rotate on the last two columns.
            [$s01, $s11] = [$s01 * $c2 - $sn2 * $s11, $sn2 * $s01 + $s11 * $c2];
            [$s02, $s12] = [$s02 * $c2 - $sn2 * $s12, $sn2 * $s02 + $s12 * $c2];
            $s00 = $norm2;

            if ($this->vT !== null) {
                if ($this->upper) {
                    $this->rotateVt($c1, $sn1, $k, $k + 1);
                } else {
                    $this->rotateVt($c2, $sn2, $k, $k + 1);
                }
            }
            if ($this->u !== null) {
                if ($this->upper) {
                    $this->rotateU($c2, -$sn2, $k, $k + 1);
                } else {
                    $this->rotateU($c1, -$sn1, $k, $k + 1);
                }
            }

            $d[$k] = $s00;
            $d[$k + 1] = $s11;
            $f[$k] = $s01;
            if ($k !== $n - 1) {
                $f[$k + 1] = $s12;
            }
            $vx = $s01;
            $vy = $s02;
        }
    }

    /** The remaining 2x2 subproblem at [$start, $start + 1]. */
    private function solve2x2(int $start, bool $computeU, bool $computeV): void
    {
        $upper = $this->upper;
        [$u2, $v1, $v2, $v2rot] = self::uptrig2x2(
            $this->d[$start],
            $this->f[$start],
            $this->d[$start + 1],
            ($computeU && $upper) || ($computeV && ! $upper),
            ($computeV && $upper) || ($computeU && ! $upper),
        );
        $this->d[$start] = $v1;
        $this->d[$start + 1] = $v2;
        $this->f[$start] = 0.0;
        if ($this->u !== null) {
            [$c, $s] = $upper ? $u2 : $v2rot;
            $this->rotateU($c, $s, $start, $start + 1);
        }
        if ($this->vT !== null) {
            [$c, $s] = $upper ? $v2rot : $u2;
            $this->rotateVt($c, -$s, $start, $start + 1);
        }
    }

    /**
     * compute_2x2_uptrig_svd (Qiao and Wang). Returns the left rotation, the two
     * values, and the right rotation; rotations are [c, s] or null.
     *
     * @return array{array{float, float}|null, float, float, array{float, float}|null}
     */
    private static function uptrig2x2(float $m11, float $m12, float $m22, bool $computeU, bool $computeV): array
    {
        $denom = hypot($m11 + $m22, $m12) + hypot($m11 - $m22, $m12);
        $v1 = fdiv($m11 * $m22 * 2.0, $denom);
        $v2 = 0.5 * $denom;
        $u = null;
        $vt = null;
        if ($computeU || $computeV) {
            [$cv, $sv, $sgnV] = self::givens($m11 * $m12, $v1 * $v1 - $m11 * $m11);
            $v1 *= $sgnV;
            $v2 *= $sgnV;
            if ($computeV) {
                $vt = [$cv, $sv];
            }
            $cu = fdiv($m11 * $cv + $m12 * $sv, $v1);
            $su = fdiv($m22 * $sv, $v1);
            [$cu, $su, $sgnU] = self::givens($cu, $su);
            $v1 *= $sgnU;
            $v2 *= $sgnU;
            if ($computeU) {
                $u = [$cu, $su];
            }
        }

        return [$u, $v1, $v2, $vt];
    }

    /**
     * GivensRotation::new(c, s): try_new with eps 0, else the identity and 0.
     *
     * @return array{float, float, float} c, s, and the signed norm
     */
    public static function givens(float $c, float $s, float $eps = 0.0): array
    {
        return self::tryGivens($c, $s, $eps) ?? [1.0, 0.0, 0.0];
    }

    /** @return array{float, float, float}|null */
    public static function tryGivens(float $c, float $s, float $eps): ?array
    {
        [$mod0, $sign0] = $c >= 0.0 ? [$c, 1.0] : [-$c, -1.0];
        $denom = sqrt($mod0 * $mod0 + $s * $s);
        if ($denom > $eps) {
            $norm = $sign0 * $denom;

            return [fdiv($mod0, $denom), fdiv($s, $norm), $norm];
        }

        return null;
    }

    /**
     * GivensRotation::cancel_y: the rotation zeroing y, or null when y is zero.
     *
     * @return array{float, float, float}|null c, s, r
     */
    public static function cancelY(float $x, float $y): ?array
    {
        if ($y == 0.0) {
            return null;
        }
        [$mod0, $sign0] = $x >= 0.0 ? [$x, 1.0] : [-$x, -1.0];
        $denom = sqrt($mod0 * $mod0 + $y * $y);
        $c = fdiv($mod0, $denom);
        $s = fdiv(-$y, $sign0 * $denom);

        return [$c, $s, $sign0 * $denom];
    }

    /**
     * GivensRotation::cancel_x: the rotation zeroing x, or null when x is zero.
     *
     * @return array{float, float, float}|null c, s, r
     */
    private static function cancelX(float $x, float $y): ?array
    {
        if ($x == 0.0) {
            return null;
        }
        [$mod1, $sign1] = $y >= 0.0 ? [$y, 1.0] : [-$y, -1.0];
        $denom = sqrt($mod1 * $mod1 + $x * $x);
        $c = fdiv($mod1, $denom);
        $s = fdiv($x * $sign1, $denom);

        return [$c, $s, $sign1 * $denom];
    }

    /** GivensRotation (c, s) rotate_rows on columns $a and $b of U. */
    private function rotateU(float $c, float $s, int $a, int $b): void
    {
        $x = $this->u[$a];
        $y = $this->u[$b];
        foreach ($x as $r => $xr) {
            $yr = $y[$r];
            $x[$r] = $xr * $c + $s * $yr;
            $y[$r] = (-$s) * $xr + $yr * $c;
        }
        $this->u[$a] = $x;
        $this->u[$b] = $y;
    }

    /** GivensRotation (c, s) rotate on rows $a and $b of V^T. */
    private function rotateVt(float $c, float $s, int $a, int $b): void
    {
        foreach ($this->vT as $j => $column) {
            $xr = $column[$a];
            $yr = $column[$b];
            $column[$a] = $xr * $c - $s * $yr;
            $column[$b] = $s * $xr + $yr * $c;
            $this->vT[$j] = $column;
        }
    }

    /** @return array{int, int} */
    private function delimit(int $end, float $eps): array
    {
        $d = &$this->d;
        $f = &$this->f;
        $n = $end;
        while ($n > 0) {
            $m = $n - 1;
            if ($f[$m] == 0.0 || abs($f[$m]) <= $eps * (abs($d[$n]) + abs($d[$m]))) {
                $f[$m] = 0.0;
            } elseif (abs($d[$m]) <= $eps) {
                $d[$m] = 0.0;
                $this->cancelHorizontal($m, $m + 1);
                if ($m !== 0) {
                    $this->cancelVertical($m - 1);
                }
            } elseif (abs($d[$n]) <= $eps) {
                $d[$n] = 0.0;
                $this->cancelVertical($m);
            } else {
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
            if (abs($f[$m]) <= $eps * (abs($d[$newStart]) + abs($d[$m]))) {
                $f[$m] = 0.0;
                break;
            } elseif (abs($d[$m]) <= $eps) {
                $d[$m] = 0.0;
                $this->cancelHorizontal($m, $n);
                if ($m !== 0) {
                    $this->cancelVertical($m - 1);
                }
                break;
            }
            $newStart--;
        }

        return [$newStart, $n];
    }

    private function cancelHorizontal(int $i, int $end): void
    {
        $d = &$this->d;
        $f = &$this->f;
        $vx = $f[$i];
        $vy = $d[$i + 1];
        $f[$i] = 0.0;
        for ($k = $i; $k < $end; $k++) {
            $rot = self::cancelX($vx, $vy);
            if ($rot === null) {
                break;
            }
            [$c, $s, $norm] = $rot;
            $d[$k + 1] = $norm;
            if ($this->upper) {
                if ($this->u !== null) {
                    $this->rotateU($c, -$s, $i, $k + 1);
                }
            } elseif ($this->vT !== null) {
                $this->rotateVt($c, $s, $i, $k + 1);
            }
            if ($k + 1 !== $end) {
                $vx = (-$s) * $f[$k + 1];
                $vy = $d[$k + 2];
                $f[$k + 1] *= $c;
            }
        }
    }

    private function cancelVertical(int $i): void
    {
        $d = &$this->d;
        $f = &$this->f;
        $vx = $d[$i];
        $vy = $f[$i];
        $f[$i] = 0.0;
        for ($k = $i; $k >= 0; $k--) {
            $rot = self::cancelY($vx, $vy);
            if ($rot === null) {
                break;
            }
            [$c, $s, $norm] = $rot;
            $d[$k] = $norm;
            if ($this->upper) {
                if ($this->vT !== null) {
                    $this->rotateVt($c, $s, $k, $i + 1);
                }
            } elseif ($this->u !== null) {
                $this->rotateU($c, -$s, $k, $i + 1);
            }
            if ($k > 0) {
                $vx = $d[$k - 1];
                $vy = $s * $f[$k - 1];
                $f[$k - 1] *= $c;
            }
        }
    }

    /**
     * SVD::sort_by_singular_values: Rust's sort_unstable_by on (value, index)
     * pairs in descending order (exact ties land where Rust puts them), then
     * the cycle-walking swaps of nalgebra's PermutationSequence applied to U's
     * columns and V^T's rows.
     */
    private function sortBySingularValues(): void
    {
        $pairs = [];
        foreach ($this->singularValues as $i => $value) {
            // Rust's sort compares every element of a slice of two or more,
            // and partial_cmp(..).expect() panics on NaN.
            if (is_nan($value) && count($this->singularValues) > 1) {
                throw new \RuntimeException('Singular value was NaN');
            }
            $pairs[] = [$value, $i];
        }
        $pairs = RustSort::sortUnstableBy($pairs, fn (array $x, array $y): bool => $y[0] < $x[0]);
        $this->singularValues = array_column($pairs, 0);

        $processed = -1;
        $swaps = [];
        $count = count($pairs);
        for ($i = 0; $i < $count; $i++) {
            $index1 = $i;
            $index2 = $pairs[$i][1];
            while ($index2 !== $processed && $pairs[$index2][1] !== $processed) {
                if ($index1 !== $index2) {
                    $swaps[] = [$index1, $index2];
                }
                $pairs[$index1][1] = $processed;
                $index1 = $index2;
                $index2 = $pairs[$index1][1];
            }
        }
        foreach ($swaps as [$a, $b]) {
            if ($this->u !== null) {
                [$this->u[$a], $this->u[$b]] = [$this->u[$b], $this->u[$a]];
            }
            if ($this->vT !== null) {
                foreach ($this->vT as $j => $column) {
                    [$column[$a], $column[$b]] = [$column[$b], $column[$a]];
                    $this->vT[$j] = $column;
                }
            }
        }
    }
}
