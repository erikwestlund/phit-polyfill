<?php

declare(strict_types=1);

namespace Phit\Polyfill;

use Phit\Polyfill\Internal\Dense;
use Phit\Polyfill\Internal\LeastSquares;

/**
 * Posterior-draw kernels: Phit\draws_moments, Phit\draws_covariance_map, and
 * Phit\draws_affine (stats-core draws.rs and the extension's draws.rs glue).
 * The factorizations follow nalgebra's operation order.
 */
final class Draws
{
    /** Phit\draws_moments */
    public static function moments(array $draws, ?array $weights = null, string $normalization = 'population'): array
    {
        $unbiased = match ($normalization) {
            'population' => false,
            'unbiased' => true,
            default => throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "normalization must be 'population' or 'unbiased'"),
        };
        $draws = self::matrix($draws);
        $w = $weights === null ? null : self::vector($weights);

        [$n, $p] = self::shape($draws);
        if ($w === null) {
            $w = array_fill(0, $n, 1.0 / $n);
        } else {
            if (count($w) !== $n) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'weights must have one element per draw');
            }
            foreach ($w as $v) {
                if ($v < 0.0) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'weights must be nonnegative');
                }
            }
            $largest = 0.0;
            foreach ($w as $v) {
                // f64::max
                $largest = $v > $largest ? $v : $largest;
            }
            if ($largest == 0.0) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'weights must have positive total weight');
            }
            $total = self::sum(array_map(fn (float $v): float => $v / $largest, $w));
            $w = array_map(fn (float $v): float => ($v / $largest) / $total, $w);
        }
        $denominator = $unbiased ? 1.0 - self::sum(array_map(fn (float $v): float => $v * $v, $w)) : 1.0;
        if ($denominator <= 0.0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'unbiased covariance normalization is undefined; need at least two effective draws');
        }
        $anchor = 0;
        while ($w[$anchor] <= 0.0) {
            $anchor++;
        }
        $mean = [];
        $centered = array_fill(0, $p, array_fill(0, $n, 0.0));
        for ($j = 0; $j < $p; $j++) {
            $origin = $draws[$anchor][$j];
            $terms = [];
            for ($i = 0; $i < $n; $i++) {
                if ($w[$i] > 0.0) {
                    $terms[] = $w[$i] * ($draws[$i][$j] - $origin);
                }
            }
            $m = self::finite($origin + self::sum($terms));
            $terms = [];
            for ($i = 0; $i < $n; $i++) {
                if ($w[$i] > 0.0) {
                    $terms[] = $w[$i] * ($draws[$i][$j] - $m);
                }
            }
            $m = self::finite($m + self::sum($terms));
            $mean[] = $m;
            for ($i = 0; $i < $n; $i++) {
                if ($w[$i] > 0.0) {
                    $centered[$j][$i] = self::finite($draws[$i][$j] - $m);
                }
            }
        }
        $covariance = array_fill(0, $p, array_fill(0, $p, 0.0));
        for ($j = 0; $j < $p; $j++) {
            for ($k = 0; $k <= $j; $k++) {
                $terms = [];
                for ($i = 0; $i < $n; $i++) {
                    $terms[] = ($w[$i] * $centered[$j][$i]) * $centered[$k][$i];
                }
                $c = self::finite(self::sum($terms) / $denominator);
                $covariance[$j][$k] = $c;
                $covariance[$k][$j] = $c;
            }
        }

        return ['mean' => $mean, 'covariance' => $covariance];
    }

    /** Phit\draws_covariance_map */
    public static function covarianceMap(array $source_covariance, array $target_covariance): array
    {
        $source = self::covarianceFactor(self::matrix($source_covariance));
        $target = self::covarianceFactor(self::matrix($target_covariance));
        $p = count($source);
        if ($p !== count($target)) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'covariances must have equal dimensions');
        }
        // L_source^T A^T = L_target^T, one column of L_target^T (a row of
        // L_target) at a time, as nalgebra's solve_upper_triangular does.
        $a = [];
        for ($c = 0; $c < $p; $c++) {
            $b = $target[$c];
            for ($i = $p - 1; $i >= 0; $i--) {
                // U = L_source^T, so U(k, i) = L_source(i, k).
                $diagonal = $source[$i][$i];
                if ($diagonal == 0.0) {
                    throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'covariance triangular solve failed');
                }
                $coeff = $b[$i] / $diagonal;
                $b[$i] = $coeff;
                for ($k = 0; $k < $i; $k++) {
                    $b[$k] = -$coeff * $source[$i][$k] + $b[$k];
                }
            }
            // Column c of A^T is row c of A.
            $a[$c] = $b;
        }
        foreach ($a as $row) {
            foreach ($row as $v) {
                self::finite($v);
            }
        }
        $terms = [];
        for ($j = 0; $j < $p; $j++) {
            $terms[] = log($target[$j][$j]) - log($source[$j][$j]);
        }

        return ['matrix' => $a, 'log_abs_det' => self::finite(self::sum($terms))];
    }

    /** Phit\draws_affine */
    public static function affine(array $draws, array $matrix, array $shift): array
    {
        $draws = self::matrix($draws);
        $transform = self::matrix($matrix);
        $shift = self::vector($shift);
        [$n, $p] = self::shape($draws);
        [$r, $c] = self::shape($transform);
        if ($r !== $c || $c !== $p || count($shift) !== $p) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'affine matrix must be square and match draws and shift dimensions');
        }
        $logAbsDet = self::affineLogDet($transform);
        // nalgebra multiplies with matrixmultiply's fused kernel once every
        // dimension exceeds 5, and with plain column updates (gemv) otherwise.
        $fused = $n > 5 && $p > 5;
        $output = [];
        foreach ($draws as $i => $row) {
            $out = [];
            for ($j = 0; $j < $p; $j++) {
                $aj = $transform[$j];
                if ($fused) {
                    $total = 0.0;
                    // k is split into blocks of 256 (KC); each block is accumulated
                    // from zero and then added to the running result.
                    for ($start = 0; $start < $p; $start += 256) {
                        $acc = 0.0;
                        $end = min($p, $start + 256);
                        for ($k = $start; $k < $end; $k++) {
                            $acc = Dense::fma($row[$k], $aj[$k], $acc);
                        }
                        $total = $start === 0 ? $acc : $acc + $total;
                    }
                    $value = $total;
                } else {
                    $value = $row[0] * $aj[0];
                    for ($k = 1; $k < $p; $k++) {
                        $value = $row[$k] * $aj[$k] + $value;
                    }
                }
                $out[] = $value + $shift[$j];
            }
            $output[] = $out;
        }
        foreach ($output as $row) {
            foreach ($row as $v) {
                self::finite($v);
            }
        }

        return ['draws' => $output, 'log_abs_det' => $logAbsDet];
    }

    /**
     * Lower Cholesky factor of a covariance, computed on the diagonal-
     * standardized matrix (symmetrized within 64 EPS) and scaled back.
     *
     * @return list<list<float>> rows of L
     */
    private static function covarianceFactor(array $covariance): array
    {
        [$n, $p] = self::shape($covariance);
        if ($n !== $p) {
            throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'covariance matrices must be square');
        }
        for ($j = 0; $j < $p; $j++) {
            if ($covariance[$j][$j] <= 0.0) {
                throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'covariance must be positive definite');
            }
        }
        $scales = [];
        for ($j = 0; $j < $p; $j++) {
            $scales[] = sqrt($covariance[$j][$j]);
        }
        // Column-major: $m[$col][$row].
        $m = [];
        for ($j = 0; $j < $p; $j++) {
            $column = array_fill(0, $p, 0.0);
            $column[$j] = 1.0;
            $m[$j] = $column;
        }
        for ($i = 0; $i < $p; $i++) {
            for ($j = 0; $j < $i; $j++) {
                $a = self::finite($covariance[$i][$j] / $scales[$i] / $scales[$j]);
                $b = self::finite($covariance[$j][$i] / $scales[$i] / $scales[$j]);
                if (abs($a - $b) > 64.0 * PHP_FLOAT_EPSILON * max(abs($a), abs($b), 1.0)) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'covariance matrices must be symmetric within tolerance');
                }
                $value = 0.5 * $a + 0.5 * $b;
                $m[$j][$i] = $value;
                $m[$i][$j] = $value;
            }
        }
        // nalgebra Cholesky::new, column by column.
        for ($j = 0; $j < $p; $j++) {
            for ($k = 0; $k < $j; $k++) {
                $factor = -$m[$k][$j];
                for ($i = $j; $i < $p; $i++) {
                    $m[$j][$i] = $factor * $m[$k][$i] + $m[$j][$i];
                }
            }
            $diagonal = $m[$j][$j];
            if ($diagonal == 0.0 || is_nan($diagonal) || $diagonal < 0.0) {
                throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'covariance must be positive definite');
            }
            $denominator = sqrt($diagonal);
            $m[$j][$j] = $denominator;
            for ($i = $j + 1; $i < $p; $i++) {
                $m[$j][$i] /= $denominator;
            }
        }
        $lower = [];
        for ($i = 0; $i < $p; $i++) {
            $row = array_fill(0, $p, 0.0);
            for ($j = 0; $j <= $i; $j++) {
                $row[$j] = $m[$j][$i];
            }
            $lower[$i] = $row;
        }
        for ($i = 0; $i < $p; $i++) {
            if ($lower[$i][$i] * $lower[$i][$i] <= 64.0 * PHP_FLOAT_EPSILON * $p) {
                throw ErrorCode::exception(ErrorCode::FACTORIZATION_FAILURE, 'covariance is numerically singular within tolerance');
            }
            for ($j = 0; $j <= $i; $j++) {
                $lower[$i][$j] = self::finite($lower[$i][$j] * $scales[$i]);
            }
        }

        return $lower;
    }

    /**
     * log|det A| after row and column equilibration, from the pivots of a
     * partially pivoted LU; numerical rank is checked on the singular values.
     */
    private static function affineLogDet(array $a): float
    {
        $p = count($a);
        $logScale = 0.0;
        for ($i = 0; $i < $p; $i++) {
            $scale = 0.0;
            foreach ($a[$i] as $v) {
                $scale = self::fmax($scale, abs($v));
            }
            if ($scale == 0.0) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'affine matrix must be invertible');
            }
            for ($j = 0; $j < $p; $j++) {
                $a[$i][$j] /= $scale;
            }
            $logScale += log($scale);
        }
        for ($j = 0; $j < $p; $j++) {
            $scale = 0.0;
            for ($i = 0; $i < $p; $i++) {
                $scale = self::fmax($scale, abs($a[$i][$j]));
            }
            if ($scale == 0.0) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'affine matrix must be invertible');
            }
            for ($i = 0; $i < $p; $i++) {
                $a[$i][$j] /= $scale;
            }
            $logScale += log($scale);
        }
        $columns = [];
        for ($j = 0; $j < $p; $j++) {
            $columns[$j] = array_column($a, $j);
        }
        $singular = LeastSquares::singularValuesOf($columns, $p);
        if (end($singular) <= 64.0 * PHP_FLOAT_EPSILON * $p * $singular[0]) {
            throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'affine matrix is numerically singular within tolerance');
        }
        // nalgebra LU::new with partial pivoting; $m[$col][$row].
        $m = $columns;
        for ($i = 0; $i < $p; $i++) {
            $piv = $i;
            $max = abs($m[$i][$i]);
            for ($r = $i + 1; $r < $p; $r++) {
                if (abs($m[$i][$r]) > $max) {
                    $max = abs($m[$i][$r]);
                    $piv = $r;
                }
            }
            $diagonal = $m[$i][$piv];
            if ($diagonal == 0.0) {
                continue;
            }
            if ($piv !== $i) {
                for ($c = 0; $c < $p; $c++) {
                    [$m[$c][$i], $m[$c][$piv]] = [$m[$c][$piv], $m[$c][$i]];
                }
            }
            $inverse = 1.0 / $diagonal;
            for ($r = $i + 1; $r < $p; $r++) {
                $m[$i][$r] *= $inverse;
            }
            for ($c = $i + 1; $c < $p; $c++) {
                $factor = -$m[$c][$i];
                for ($r = $i + 1; $r < $p; $r++) {
                    $m[$c][$r] = $factor * $m[$i][$r] + $m[$c][$r];
                }
            }
        }
        $terms = [];
        for ($j = 0; $j < $p; $j++) {
            $terms[] = log(abs($m[$j][$j]));
        }

        return self::finite($logScale + self::sum($terms));
    }

    /** Rust's f64::max: NaN loses to a number. */
    private static function fmax(float $a, float $b): float
    {
        if (is_nan($a)) {
            return $b;
        }

        return $b > $a ? $b : $a;
    }

    /** Neumaier summation, as in draws.rs. */
    private static function sum(array $values): float
    {
        $total = 0.0;
        $correction = 0.0;
        foreach ($values as $v) {
            $t = $total + $v;
            $correction += abs($total) >= abs($v) ? ($total - $t) + $v : ($v - $t) + $total;
            $total = $t;
        }

        return $total + $correction;
    }

    private static function finite(float $value): float
    {
        if (is_finite($value)) {
            return $value;
        }
        throw ErrorCode::exception(ErrorCode::NONFINITE, 'draw operation exceeded floating-point range');
    }

    /** @return array{0: int, 1: int} */
    private static function shape(array $rows): array
    {
        if ($rows === [] || $rows[0] === []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'matrices must have at least one row and column');
        }
        $p = count($rows[0]);
        foreach ($rows as $row) {
            if (count($row) !== $p) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, 'matrices must be rectangular');
            }
        }

        return [count($rows), $p];
    }

    /** Integers or finite floats only; null and other types are rejected. */
    private static function number(mixed $value): float
    {
        if (is_int($value)) {
            return (float) $value;
        }
        if (! is_float($value)) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'draw inputs must contain integers or floats; null is not allowed');
        }
        if (! is_finite($value)) {
            throw ErrorCode::exception(ErrorCode::NONFINITE, 'draw inputs must contain finite numbers');
        }

        return $value;
    }

    /** @return list<float> */
    private static function vector(array $input): array
    {
        $out = [];
        foreach ($input as $value) {
            $out[] = self::number($value);
        }

        return $out;
    }

    /** @return list<list<float>> */
    private static function matrix(array $input): array
    {
        $out = [];
        foreach ($input as $row) {
            if (! is_array($row)) {
                throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'each matrix row must be an array');
            }
            $out[] = self::vector($row);
        }

        return $out;
    }
}
