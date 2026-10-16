<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/**
 * Least squares through Householder QR of the n×p design followed by a
 * one-sided Jacobi SVD of the small p×p factor R. The singular values of R
 * are those of the design, so the extension's rank rule
 * (s <= EPS * max(n, p) * s_max) gives the same decisions without an n×p SVD.
 *
 * Matrices are column-major: $columns[$j][$i] is row i of column j.
 *
 * @internal
 */
final class LeastSquares
{
    /** @var list<float> singular values, descending */
    public array $singular;

    /** @var list<list<float>> right singular vectors: $v[$i][$k] = V(i, k) */
    public array $v;

    /** @var list<float> U' Q' b for the supplied right-hand side */
    private array $uqtb = [];

    public readonly float $tolerance;

    /**
     * @param  list<list<float>>  $columns
     * @param  list<float>|null  $rhs
     */
    public function __construct(array $columns, int $n, int $p, ?array $rhs = null)
    {
        $qtb = $rhs;
        for ($k = 0; $k < $p; $k++) {
            $ck = &$columns[$k];
            $norm2 = 0.0;
            for ($i = $k; $i < $n; $i++) {
                $norm2 += $ck[$i] * $ck[$i];
            }
            if ($norm2 == 0.0) {
                unset($ck);

                continue;
            }
            $norm = sqrt($norm2);
            $alpha = $ck[$k] >= 0.0 ? -$norm : $norm;
            // v = x - alpha e1 stored in place of the column below the diagonal.
            $ck[$k] -= $alpha;
            $vnorm2 = $norm2 - ($ck[$k] + $alpha) * ($ck[$k] + $alpha) + $ck[$k] * $ck[$k];
            $tau = 2.0 / $vnorm2;
            for ($j = $k + 1; $j < $p; $j++) {
                $cj = &$columns[$j];
                $dot = 0.0;
                for ($i = $k; $i < $n; $i++) {
                    $dot += $ck[$i] * $cj[$i];
                }
                $f = $tau * $dot;
                for ($i = $k; $i < $n; $i++) {
                    $cj[$i] -= $f * $ck[$i];
                }
                unset($cj);
            }
            if ($qtb !== null) {
                $dot = 0.0;
                for ($i = $k; $i < $n; $i++) {
                    $dot += $ck[$i] * $qtb[$i];
                }
                $f = $tau * $dot;
                for ($i = $k; $i < $n; $i++) {
                    $qtb[$i] -= $f * $ck[$i];
                }
            }
            $ck[$k] = $alpha;
            unset($ck);
        }

        // R as columns of a p×p upper-triangular matrix.
        $r = [];
        for ($j = 0; $j < $p; $j++) {
            $column = array_fill(0, $p, 0.0);
            for ($i = 0; $i <= $j; $i++) {
                $column[$i] = $columns[$j][$i];
            }
            $r[$j] = $column;
        }
        $this->jacobi($r, $p, $qtb === null ? null : array_slice($qtb, 0, $p));
        $this->tolerance = PHP_FLOAT_EPSILON * max($n, $p) * $this->singular[0];
    }

    /** One-sided Jacobi: R V = U S with orthogonal V, then sort by s descending. */
    private function jacobi(array $a, int $p, ?array $c): void
    {
        $v = [];
        for ($j = 0; $j < $p; $j++) {
            $v[$j] = array_fill(0, $p, 0.0);
            $v[$j][$j] = 1.0;
        }
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
                    if ($gamma == 0.0 || abs($gamma) <= PHP_FLOAT_EPSILON * sqrt($alpha * $beta)) {
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
                    $vi = $v[$i];
                    $vj = $v[$j];
                    for ($k = 0; $k < $p; $k++) {
                        $x = $vi[$k];
                        $y = $vj[$k];
                        $vi[$k] = $cs * $x - $sn * $y;
                        $vj[$k] = $sn * $x + $cs * $y;
                    }
                    $v[$i] = $vi;
                    $v[$j] = $vj;
                }
            }
            if (! $rotated) {
                break;
            }
        }
        $s = [];
        for ($j = 0; $j < $p; $j++) {
            $norm2 = 0.0;
            foreach ($a[$j] as $value) {
                $norm2 += $value * $value;
            }
            $s[$j] = sqrt($norm2);
        }
        $order = range(0, $p - 1);
        usort($order, fn (int $x, int $y): int => $s[$y] <=> $s[$x]);
        $this->singular = [];
        $this->v = array_fill(0, $p, []);
        foreach ($order as $k => $j) {
            $this->singular[$k] = $s[$j];
            // $v[$j] is column j of V, i.e. V(i, j) = $v[$j][$i].
            for ($i = 0; $i < $p; $i++) {
                $this->v[$i][$k] = $v[$j][$i];
            }
            if ($c !== null) {
                $dot = 0.0;
                if ($s[$j] > 0.0) {
                    for ($i = 0; $i < $p; $i++) {
                        $dot += $a[$j][$i] / $s[$j] * $c[$i];
                    }
                }
                $this->uqtb[$k] = $dot;
            }
        }
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

    /** Least-squares solution for the right-hand side given to the constructor. */
    public function solution(): array
    {
        $p = count($this->singular);
        $x = array_fill(0, $p, 0.0);
        for ($k = 0; $k < $p; $k++) {
            if ($this->singular[$k] <= $this->tolerance) {
                continue;
            }
            $w = $this->uqtb[$k] / $this->singular[$k];
            for ($i = 0; $i < $p; $i++) {
                $x[$i] += $this->v[$i][$k] * $w;
            }
        }

        return $x;
    }

    /**
     * (A'A)^-1 with the column scales undone: V diag(1/s^2) V' / (scale_i scale_j).
     *
     * @return list<list<float>>
     */
    public function inverseGram(array $scales, float $factor = 1.0): array
    {
        $p = count($this->singular);
        $out = [];
        for ($i = 0; $i < $p; $i++) {
            $row = [];
            for ($j = 0; $j < $p; $j++) {
                $sum = 0.0;
                for ($k = 0; $k < $p; $k++) {
                    $sum += $this->v[$i][$k] * $this->v[$j][$k] / ($this->singular[$k] * $this->singular[$k]);
                }
                $row[$j] = $sum / $scales[$i] / $scales[$j] * $factor;
            }
            $out[$i] = $row;
        }

        return $out;
    }
}
