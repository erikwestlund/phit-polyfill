<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

use PhpStats\Polyfill\ErrorCode;

/**
 * The extension's damped Newton / Fisher-scoring loop for the logistic,
 * probit and Poisson models. Each iteration solves a weighted least-squares
 * problem for the step (QR, then SVD of R), halves the step until the
 * likelihood does not decrease, and stops only after a small full step.
 * Families supply the vectorised likelihood, weights and score.
 *
 * @internal
 */
abstract class Glm
{
    /** Name used in messages: logistic, probit or poisson. */
    abstract protected function name(): string;

    abstract protected function solver(): string;

    /** Validate outcomes and weights; return canonical [y, weights, offset]. */
    abstract protected function prepareOutcomes(array $y, array $weights, array $offset, int $n): array;

    /** Starting theta (in scaled units) or null for zeros. */
    abstract protected function start(array $y, array $weights, array $offset, array $design, array $scales, int $n, int $p): ?array;

    abstract protected function likelihood(array $eta, array $y, array $weights): float;

    /** sqrt of the working weight per row (throws for non-finite where the extension does). */
    abstract protected function rootWeights(array $eta): array;

    /** Per-row score d loglik / d eta (before prior weights). */
    abstract protected function scores(array $eta, array $y): array;

    abstract protected function rankMessage(): string;

    abstract protected function overflowMessage(): string;

    abstract protected function nonconvergenceMessage(int $max): string;

    /** @return array{0: float, 1: float} step and likelihood convergence tolerances */
    abstract protected function tolerances(): array;

    protected function separated(array $linear, array $y): bool
    {
        return false;
    }

    /** sqrt working weights at the final estimate. */
    abstract protected function finalRootWeights(array $eta): array;

    /** A scaled design entry times its final root weight, in the family's operation order. */
    protected function finalWeighted(float $value, float $root, float $eta): float
    {
        return $value * $root;
    }

    abstract protected function fittedValues(array $eta): array;

    /** Model-specific likelihood statistics: log_likelihood, deviance, null_deviance, aic. */
    abstract protected function statistics(
        array $eta, array $y, array $weights, array $offset, float $likelihood,
        bool $intercept, bool $computeNull, int $n, int $p, int $maxIterations,
    ): array;

    protected static function error(int $code, string $message): \Exception
    {
        return ErrorCode::exception($code, $message);
    }

    /** Linear predictor sum_j column_j * theta_j + offset, in nalgebra's accumulation order. */
    private static function linear(array $columns, array $theta, array $offset, int $n): array
    {
        $eta = [];
        $c0 = $columns[0];
        $t0 = $theta[0];
        for ($i = 0; $i < $n; $i++) {
            $eta[$i] = $c0[$i] * $t0;
        }
        $p = count($theta);
        for ($j = 1; $j < $p; $j++) {
            $cj = $columns[$j];
            $tj = $theta[$j];
            for ($i = 0; $i < $n; $i++) {
                $eta[$i] += $cj[$i] * $tj;
            }
        }
        for ($i = 0; $i < $n; $i++) {
            $eta[$i] += $offset[$i];
        }

        return $eta;
    }

    /**
     * @param  list<list<float>>  $x  rows including any intercept column
     * @return array<string, mixed>
     */
    public function fitData(array $x, array $y, array $weights, array $offset, int $maxIterations = 100, bool $computeNull = true): array
    {
        $n = count($x);
        if ($n === 0 || $x[0] === []) {
            throw self::error(ErrorCode::INVALID_INPUT, 'x must have at least one row and column');
        }
        $p = count($x[0]);
        if (count($y) !== $n) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, 'x must be rectangular with one y value per row');
        }
        foreach ($x as $row) {
            if (count($row) !== $p) {
                throw self::error(ErrorCode::DIMENSION_MISMATCH, 'x must be rectangular with one y value per row');
            }
        }
        if ($n <= $p) {
            throw self::error(ErrorCode::DIMENSION_MISMATCH, $this->name().' regression requires n > p');
        }
        foreach ($x as $row) {
            foreach ($row as $value) {
                if (! is_finite($value)) {
                    throw self::error(ErrorCode::NONFINITE, 'x and y must be finite');
                }
            }
        }
        foreach ($y as $value) {
            if (! is_finite($value)) {
                throw self::error(ErrorCode::NONFINITE, 'x and y must be finite');
            }
        }
        [$y, $weights, $offset] = $this->prepareOutcomes($y, $weights, $offset, $n);

        $design = [];
        $columns = [];
        $scales = [];
        $intercept = false;
        for ($j = 0; $j < $p; $j++) {
            $column = array_column($x, $j);
            $ones = true;
            foreach ($column as $value) {
                $ones = $ones && $value == 1.0;
            }
            $intercept = $intercept || $ones;
            // nalgebra's column norm, summed the same way.
            $scale = sqrt(LeastSquares::dot($column, $column, 0, $n));
            if (! is_finite($scale) || $scale == 0.0) {
                throw self::error(ErrorCode::RANK_DEFICIENCY, 'design is rank deficient or exceeds floating-point range');
            }
            $design[$j] = $column;
            foreach ($column as $i => $value) {
                $column[$i] = $value / $scale;
            }
            $columns[$j] = $column;
            $scales[$j] = $scale;
        }
        $rootPrior = array_map(sqrt(...), $weights);

        $theta = array_fill(0, $p, 0.0);
        $eta = $offset;
        $started = $this->start($y, $weights, $offset, $design, $scales, $n, $p);
        if ($started !== null) {
            $theta = $started;
            $eta = self::linear($columns, $theta, $offset, $n);
        }
        $likelihood = $this->likelihood($eta, $y, $weights);
        [$stepTolerance, $likelihoodTolerance] = $this->tolerances();
        $convergedAt = null;
        for ($iteration = 1; $iteration <= $maxIterations; $iteration++) {
            $root = $this->rootWeights($eta);
            $weighted = [];
            for ($j = 0; $j < $p; $j++) {
                $column = $columns[$j];
                for ($i = 0; $i < $n; $i++) {
                    $column[$i] = $column[$i] * $root[$i] * $rootPrior[$i];
                }
                $weighted[$j] = $column;
            }
            $scores = $this->scores($eta, $y);
            $rhs = [];
            for ($i = 0; $i < $n; $i++) {
                $rhs[$i] = ($root[$i] == 0.0 && $scores[$i] == 0.0)
                    ? 0.0
                    : fdiv($rootPrior[$i] * $scores[$i], $root[$i]);
            }
            $ls = new LeastSquares($weighted, $n, $p, $rhs);
            if ($ls->rankDeficient()) {
                throw self::error(ErrorCode::NONCONVERGENCE, $this->rankMessage());
            }
            foreach ($rhs as $value) {
                if (! is_finite($value)) {
                    throw self::error(ErrorCode::NONCONVERGENCE, $this->overflowMessage());
                }
            }
            $step = $ls->solution()
                ?? throw self::error(ErrorCode::FACTORIZATION_FAILURE, 'least-squares factorization failed');
            foreach ($step as $value) {
                if (! is_finite($value)) {
                    throw self::error(ErrorCode::NONFINITE, $this->name().' step exceeds floating-point range');
                }
            }
            // Concave likelihood: halve an overlarge step until it does not decrease.
            $fraction = 1.0;
            $accepted = null;
            $floor = $likelihood - 16.0 * PHP_FLOAT_EPSILON * (1.0 + abs($likelihood));
            for ($attempt = 0; $attempt < 31; $attempt++) {
                $candidate = [];
                for ($j = 0; $j < $p; $j++) {
                    $candidate[$j] = $theta[$j] + $fraction * $step[$j];
                }
                $candidateEta = self::linear($columns, $candidate, $offset, $n);
                $candidateLikelihood = $this->likelihood($candidateEta, $y, $weights);
                if (is_finite($candidateLikelihood) && $candidateLikelihood >= $floor) {
                    $accepted = [$candidate, $candidateEta, $candidateLikelihood];
                    break;
                }
                $fraction *= 0.5;
            }
            if ($accepted === null) {
                throw self::error(ErrorCode::NONCONVERGENCE, $this->name().' line search failed to improve the likelihood');
            }
            [$nextTheta, $nextEta, $nextLikelihood] = $accepted;
            $linear = [];
            for ($i = 0; $i < $n; $i++) {
                $linear[$i] = $nextEta[$i] - $offset[$i];
            }
            if ($this->separated($linear, $y)) {
                throw self::error(ErrorCode::NONCONVERGENCE, 'complete or quasi-complete separation: no finite unpenalized '.$this->name().' MLE');
            }
            // Require a small full step, not merely a damped one or a flat likelihood.
            $smallStep = true;
            for ($j = 0; $j < $p; $j++) {
                if (abs($step[$j]) > $stepTolerance * (1.0 + abs($nextTheta[$j]))) {
                    $smallStep = false;
                    break;
                }
            }
            $smallChange = abs($nextLikelihood - $likelihood) <= $likelihoodTolerance * (1.0 + abs($likelihood));
            $theta = $nextTheta;
            $eta = $nextEta;
            $likelihood = $nextLikelihood;
            if ($smallStep && $smallChange) {
                $convergedAt = $iteration;
                break;
            }
        }
        if ($convergedAt === null) {
            throw self::error(ErrorCode::NONCONVERGENCE, $this->nonconvergenceMessage($maxIterations));
        }
        $coefficients = [];
        for ($j = 0; $j < $p; $j++) {
            $coefficients[$j] = $theta[$j] / $scales[$j];
        }
        // Information at the final estimate, not the preceding iterate.
        $root = $this->finalRootWeights($eta);
        $weighted = [];
        for ($j = 0; $j < $p; $j++) {
            $column = $columns[$j];
            for ($i = 0; $i < $n; $i++) {
                $column[$i] = $this->finalWeighted($column[$i], $root[$i], $eta[$i]) * $rootPrior[$i];
            }
            $weighted[$j] = $column;
        }
        $ls = new LeastSquares($weighted, $n, $p);
        if ($ls->rankDeficient()) {
            throw self::error(ErrorCode::RANK_DEFICIENCY, 'final '.$this->name().' information is numerically rank deficient');
        }
        $covariance = $ls->inverseGram($scales)
            ?? throw self::error(ErrorCode::FACTORIZATION_FAILURE, 'triangular factor could not be inverted');
        $critical = Statrs::Z_975;
        $stdErrors = $zStatistics = $pValues = $confInt = [];
        $finite = true;
        for ($j = 0; $j < $p; $j++) {
            $se = sqrt($covariance[$j][$j]);
            $z = fdiv($coefficients[$j], $se);
            $stdErrors[$j] = $se;
            $zStatistics[$j] = $z;
            $pValues[$j] = 2.0 * Statrs::normalSf(abs($z));
            $confInt[$j] = [$coefficients[$j] - $critical * $se, $coefficients[$j] + $critical * $se];
            $finite = $finite && is_finite($coefficients[$j]) && is_finite($se) && is_finite($z)
                && is_finite($pValues[$j]) && is_finite($confInt[$j][0]) && is_finite($confInt[$j][1]);
            foreach ($covariance[$j] as $value) {
                $finite = $finite && is_finite($value);
            }
        }
        if (! $finite) {
            throw self::error(ErrorCode::NONFINITE, $this->name().' inference exceeds floating-point range');
        }
        $fitted = $this->fittedValues($eta);
        $residuals = [];
        foreach ($y as $i => $value) {
            $residuals[$i] = $value - $fitted[$i];
        }
        $stats = $this->statistics($eta, $y, $weights, $offset, $likelihood, $intercept, $computeNull, $n, $p, $maxIterations);

        return [
            'coefficients' => $coefficients,
            'covariance' => $covariance,
            'std_errors' => $stdErrors,
            'z_statistics' => $zStatistics,
            'p_values' => $pValues,
            'conf_int' => $confInt,
            'fitted' => $fitted,
            'residuals' => $residuals,
            'log_likelihood' => $stats['log_likelihood'],
            'deviance' => $stats['deviance'],
            'null_deviance' => $stats['null_deviance'],
            'aic' => $stats['aic'],
            'df_residual' => $n - $p,
            'iterations' => $convergedAt,
            'converged' => true,
        ];
    }

    /** Public entry shared by the stats_* model functions. */
    public function run(array $x, array $y, string $missing, ?array $weights, ?array $offset, bool $intercept, string $outcome, string $function): array
    {
        $data = ModelData::prepare($x, $y, $missing, $outcome, $weights, $offset, $intercept);

        return $data->finish($this->fitData($data->x, $data->y, $data->weights, $data->offset), $function, $this->solver());
    }

    /** Validate prior weights as the extension's observations::prepare does. */
    protected static function checkWeights(array $weights, bool $binary): void
    {
        foreach ($weights as $w) {
            if ($w <= 0.0 || ($binary && fmod($w, 1.0) != 0.0)) {
                throw self::error(ErrorCode::INVALID_INPUT, 'weights must be positive (integer trials for binary outcomes)');
            }
        }
        if ($binary) {
            foreach ($weights as $w) {
                if ($w > Limits::MAX_TRIALS) {
                    throw self::error(ErrorCode::UNSUPPORTED, 'binomial trial weights must not exceed 1000000; larger-count likelihood arithmetic is not verified');
                }
            }
        }
    }

    protected static function softplus(float $t): float
    {
        return max($t, 0.0) + log1p(exp(-abs($t)));
    }
}
