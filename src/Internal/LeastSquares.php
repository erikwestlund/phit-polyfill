<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/**
 * Least squares through Householder QR of the n×p design, as the extension
 * does since php-stats 91ee538. The singular values of the small p×p factor R
 * (computed by a one-sided Jacobi sweep) are those of the design and drive the
 * extension's rank rule s <= EPS * max(n, p) * s_max. Solutions and (R'R)^-1
 * come from R by back-substitution, never from singular vectors or X'X.
 *
 * Matrices are column-major: $columns[$j][$i] is row i of column j.
 *
 * @internal
 */
final class LeastSquares
{
    /** @var list<float> singular values, descending */
    public array $singular;

    public readonly float $tolerance;

    /** @var list<list<float>> R as columns: $r[$j][$i] = R(i, j), zero below the diagonal */
    private array $r;

    /** @var list<float>|null the first p entries of Q'b for the supplied right-hand side */
    private ?array $qtb;

    /**
     * @param  list<list<float>>  $columns
     * @param  list<float>|null  $rhs
     */
    public function __construct(array $columns, int $n, private readonly int $p, ?array $rhs = null)
    {
        // Householder QR exactly as nalgebra's QR::new / q_tr_mul compute it (same
        // operation order, 8-way dot products, reflections applied "with sign"),
        // so R matches the extension's to the last bit where possible: near the
        // convergence threshold of an ill-conditioned fit, rounding decides.
        $qtb = $rhs;
        $diagonal = array_fill(0, $p, 0.0);
        for ($k = 0; $k < $p; $k++) {
            $ck = &$columns[$k];
            $sq = self::dot($ck, $ck, $k, $n);
            $norm = sqrt($sq);
            $x0 = $ck[$k];
            $sign = ($x0 < 0.0 || ($x0 == 0.0 && fdiv(1.0, $x0) < 0.0)) ? -1.0 : 1.0;
            $signedNorm = $sign * $norm;
            $factor = ($sq + abs($x0) * $norm) * 2.0;
            $ck[$k] += $signedNorm;
            if ($factor == 0.0) {
                $diagonal[$k] = $signedNorm;
                unset($ck);

                continue;
            }
            $root = sqrt($factor);
            for ($i = $k; $i < $n; $i++) {
                $ck[$i] /= $root;
            }
            $length = sqrt(self::dot($ck, $ck, $k, $n));
            for ($i = $k; $i < $n; $i++) {
                $ck[$i] /= $length;
            }
            $reflectionNorm = -$signedNorm;
            $diagonal[$k] = $reflectionNorm;
            $reflectSign = $reflectionNorm < 0.0 || ($reflectionNorm == 0.0 && fdiv(1.0, $reflectionNorm) < 0.0) ? -1.0 : 1.0;
            $scale = $reflectSign * -2.0;
            for ($j = $k + 1; $j < $p; $j++) {
                $cj = &$columns[$j];
                $f = (self::dot($ck, $cj, $k, $n) - 0.0) * $scale;
                for ($i = $k; $i < $n; $i++) {
                    $cj[$i] = $f * $ck[$i] + $reflectSign * $cj[$i];
                }
                unset($cj);
            }
            if ($qtb !== null) {
                $f = (self::dot($ck, $qtb, $k, $n) - 0.0) * $scale;
                for ($i = $k; $i < $n; $i++) {
                    $qtb[$i] = $f * $ck[$i] + $reflectSign * $qtb[$i];
                }
            }
            unset($ck);
        }

        // R: the upper triangle, with |diagonal| on the diagonal (nalgebra's r()).
        $r = [];
        for ($j = 0; $j < $p; $j++) {
            $column = array_fill(0, $p, 0.0);
            for ($i = 0; $i < $j; $i++) {
                $column[$i] = $columns[$j][$i];
            }
            $column[$j] = abs($diagonal[$j]);
            $r[$j] = $column;
        }
        $this->r = $r;
        $this->qtb = $qtb === null ? null : array_slice($qtb, 0, $p);
        $this->singular = self::singularValues($r, $p);
        $this->tolerance = PHP_FLOAT_EPSILON * max($n, $p) * $this->singular[0];
    }

    /**
     * Dot product of $x and $y over indices [$start, $n), accumulated as
     * nalgebra's dotx does for dynamic vectors: eight interleaved partial sums
     * over blocks of eight, combined pairwise, then the remainder in order.
     */
    public static function dot(array $x, array $y, int $start, int $n): float
    {
        $acc0 = $acc1 = $acc2 = $acc3 = $acc4 = $acc5 = $acc6 = $acc7 = 0.0;
        $i = $start;
        while ($n - $i >= 8) {
            $acc0 += $x[$i] * $y[$i];
            $acc1 += $x[$i + 1] * $y[$i + 1];
            $acc2 += $x[$i + 2] * $y[$i + 2];
            $acc3 += $x[$i + 3] * $y[$i + 3];
            $acc4 += $x[$i + 4] * $y[$i + 4];
            $acc5 += $x[$i + 5] * $y[$i + 5];
            $acc6 += $x[$i + 6] * $y[$i + 6];
            $acc7 += $x[$i + 7] * $y[$i + 7];
            $i += 8;
        }
        $result = 0.0;
        $result += $acc0 + $acc4;
        $result += $acc1 + $acc5;
        $result += $acc2 + $acc6;
        $result += $acc3 + $acc7;
        for (; $i < $n; $i++) {
            $result += $x[$i] * $y[$i];
        }

        return $result;
    }

    /**
     * Singular values of a square matrix given as columns, sorted descending.
     *
     * @param  list<list<float>>  $columns
     * @return list<float>
     */
    public static function singularValuesOf(array $columns, int $p): array
    {
        return self::singularValues($columns, $p);
    }

    /** Singular values of R by a one-sided Jacobi sweep, sorted descending. */
    private static function singularValues(array $a, int $p): array
    {
        for ($sweep = 0; $sweep < 60; $sweep++) {
            $rotated = false;
            for ($i = 0; $i < $p - 1; $i++) {
                for ($j = $i + 1; $j < $p; $j++) {
                    $ai = $a[$i];
                    $aj = $a[$j];
                    $alpha = 0.0;
                    $beta = 0.0;
                    $gamma = 0.0;
                    for ($k = 0; $k < $p; $k++) {
                        $alpha += $ai[$k] * $ai[$k];
                        $beta += $aj[$k] * $aj[$k];
                        $gamma += $ai[$k] * $aj[$k];
                    }
                    // sqrt(alpha) * sqrt(beta): the product alpha * beta overflows once the
                    // column norms pass about 1e77, and an infinite threshold would stop
                    // every rotation (heavily weighted IRLS steps reach that range).
                    if ($gamma == 0.0 || abs($gamma) <= PHP_FLOAT_EPSILON * sqrt($alpha) * sqrt($beta)) {
                        continue;
                    }
                    $rotated = true;
                    $zeta = ($beta - $alpha) / (2.0 * $gamma);
                    $t = ($zeta >= 0.0 ? 1.0 : -1.0) / (abs($zeta) + sqrt(1.0 + $zeta * $zeta));
                    $cs = 1.0 / sqrt(1.0 + $t * $t);
                    $sn = $cs * $t;
                    for ($k = 0; $k < $p; $k++) {
                        $x = $ai[$k];
                        $y = $aj[$k];
                        $ai[$k] = $cs * $x - $sn * $y;
                        $aj[$k] = $sn * $x + $cs * $y;
                    }
                    $a[$i] = $ai;
                    $a[$j] = $aj;
                }
            }
            if (! $rotated) {
                break;
            }
        }
        $s = [];
        foreach ($a as $column) {
            $norm2 = 0.0;
            foreach ($column as $value) {
                $norm2 += $value * $value;
            }
            $s[] = sqrt($norm2);
        }
        rsort($s);

        return $s;
    }

    /** Whether any singular value is at or below the rank tolerance. */
    public function rankDeficient(): bool
    {
        foreach ($this->singular as $s) {
            if ($s <= $this->tolerance) {
                return true;
            }
        }

        return false;
    }

    /**
     * Solve R x = b in place, as nalgebra's solve_upper_triangular does: from
     * the last row up, divide by the diagonal and subtract that column from the
     * rows above. Null when a diagonal entry is zero.
     *
     * @param  list<float>  $b
     * @return list<float>|null
     */
    private function backSubstitute(array $b): ?array
    {
        for ($i = $this->p - 1; $i >= 0; $i--) {
            $column = $this->r[$i];
            $diagonal = $column[$i];
            if ($diagonal == 0.0) {
                return null;
            }
            $b[$i] /= $diagonal;
            $value = $b[$i];
            for ($k = 0; $k < $i; $k++) {
                $b[$k] -= $column[$k] * $value;
            }
        }

        return $b;
    }

    /**
     * Least-squares solution for the right-hand side given to the constructor,
     * or null when R has a zero on its diagonal.
     *
     * @return list<float>|null
     */
    public function solution(): ?array
    {
        return $this->backSubstitute($this->qtb ?? throw new \LogicException('no right-hand side was supplied'));
    }

    /**
     * (A'A)^-1 = R^-1 R^-T with the column scales undone and an optional
     * factor, or null when R has a zero on its diagonal.
     *
     * @return list<list<float>>|null
     */
    public function inverseGram(array $scales, float $factor = 1.0): ?array
    {
        $p = $this->p;
        // Columns of R^-1: solve R x = e_j.
        $inverse = [];
        for ($j = 0; $j < $p; $j++) {
            $unit = array_fill(0, $p, 0.0);
            $unit[$j] = 1.0;
            $column = $this->backSubstitute($unit);
            if ($column === null) {
                return null;
            }
            $inverse[$j] = $column;
        }
        $out = [];
        for ($i = 0; $i < $p; $i++) {
            $row = [];
            for ($j = 0; $j < $p; $j++) {
                // (R^-1 R^-T)(i, j) = sum_k R^-1(i, k) R^-1(j, k)
                $sum = 0.0;
                for ($k = 0; $k < $p; $k++) {
                    $sum += $inverse[$k][$i] * $inverse[$k][$j];
                }
                $row[$j] = $sum / $scales[$i] / $scales[$j] * $factor;
            }
            $out[$i] = $row;
        }

        return $out;
    }
}
