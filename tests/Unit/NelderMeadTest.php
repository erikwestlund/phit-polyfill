<?php

declare(strict_types=1);

use PhpStats\Polyfill\Distribution\RustPanic;
use PhpStats\Polyfill\Internal\NelderMead;
use PhpStats\Polyfill\Internal\NelderMeadError;
use PhpStats\Polyfill\Internal\NelderMeadTermination;

// Bit-for-bit comparison with argmin 0.11.0's NelderMead + Executor, from the
// harness in tests/reference/neldermead. Floats are IEEE-754 bit patterns in hex.

$fixture = loadFixture('neldermead-reference.json');

$bits = static fn (float $v): string => bin2hex(pack('E', $v));
$float = static fn (?string $h): ?float => $h === null ? null : unpack('E', hex2bin($h))[1];

/**
 * The harness's cost function: a weighted quadratic or Rosenbrock, optionally
 * floored to a plateau grid, with lmm.rs-style bounds (+INF) and evaluation
 * budget (an Err), and NaN / -INF regions.
 */
$solve = static function (array $spec) use ($bits, $float): array {
    $floats = static fn (?array $l): ?array => $l === null ? null : array_map($float, $l);
    $centre = $floats($spec['centre']);
    $weights = $floats($spec['weights']);
    $lower = $floats($spec['lower']);
    $nanAbove = $float($spec['nan_above']);
    $negInfBelow = $float($spec['neginf_below']);
    $plateau = $float($spec['plateau']);
    $evaluations = 0;
    $trace = [];
    $cost = static function (array $x) use (&$evaluations, &$trace, $spec, $centre, $weights, $lower, $nanAbove, $negInfBelow, $plateau, $bits): float {
        if ($spec['budget'] !== null && $evaluations >= $spec['budget']) {
            $trace[] = [array_map($bits, $x), null];
            throw new NelderMeadError($spec['error_message']);
        }
        $evaluations++;
        $outside = false;
        foreach ($lower ?? [] as $i => $l) {
            $outside = $outside || is_nan($x[$i]) || is_infinite($x[$i]) || $x[$i] < $l || $x[$i] > 1e4;
        }
        if ($outside) {
            $c = INF;
        } elseif ($nanAbove !== null && $x[0] > $nanAbove) {
            $c = NAN;
        } elseif ($negInfBelow !== null && $x[0] < $negInfBelow) {
            $c = -INF;
        } else {
            $y = $plateau === null ? $x : array_map(static fn (float $v): float => floor($v / $plateau) * $plateau, $x);
            $n = count($y);
            $c = 0.0;
            if ($spec['kind'] === 'quad') {
                for ($i = 0; $i < $n; $i++) {
                    $d = $y[$i] - $centre[$i];
                    $c += $weights[$i] * $d * $d;
                }
            } else {
                for ($i = 0; $i < $n - 1; $i++) {
                    $a = $y[$i + 1] - $y[$i] * $y[$i];
                    $b = 1.0 - $y[$i];
                    $c += 100.0 * $a * $a + $b * $b;
                }
                $b = $centre[0] - $y[$n - 1];
                $c += $b * $b;
            }
        }
        $trace[] = [array_map($bits, $x), $bits($c)];

        return $c;
    };

    try {
        $solver = NelderMead::new(array_map($floats, $spec['simplex']))->withSdTolerance($float($spec['sd_tolerance']));
    } catch (NelderMeadError $e) {
        return [['status' => 'invalid', 'message' => $e->getMessage()], 0, []];
    }
    try {
        $r = $solver->run($cost, $spec['max_iters'] ?? PHP_INT_MAX);
        $result = [
            'status' => 'ok',
            'termination' => $r->termination->value,
            'iterations' => $r->iterations,
            'best_param' => $r->bestParam === null ? null : array_map($bits, $r->bestParam),
            'best_cost' => $bits($r->bestCost),
        ];
    } catch (NelderMeadError $e) {
        $result = ['status' => 'error', 'message' => $e->getMessage()];
    } catch (RustPanic $e) {
        $result = ['status' => 'panic', 'message' => substr($e->getMessage(), strlen('Rust panic: '))];
    }

    return [$result, $evaluations, $trace];
};

it('reproduces argmin bit for bit', function () use ($fixture, $solve) {
    $statuses = [];
    foreach ($fixture['problems'] as $p) {
        [$result, $evaluations, $trace] = $solve($p['problem']);
        expect($result)->toBe($p['result'], $p['name'])
            ->and($evaluations)->toBe($p['evaluations'], "{$p['name']}: evaluations")
            // The fixture keeps only the first 120 evaluated points of each run.
            ->and(array_slice($trace, 0, count($p['trace'])))->toBe($p['trace'], "{$p['name']}: evaluated points");
        $statuses[$result['status']] = true;
    }
    expect(array_keys($statuses))->toEqualCanonicalizing(['ok', 'error', 'panic', 'invalid'])
        ->and(count($fixture['problems']))->toBeGreaterThan(70);
});

it('validates sd_tolerance like with_sd_tolerance', function () {
    $simplex = [[0.0], [1.0]];
    expect(fn () => NelderMead::new($simplex)->withSdTolerance(-1e-9))
        ->toThrow(NelderMeadError::class, 'Invalid parameter: "`Nelder-Mead`: sd_tolerance must be >= 0."');
    $cost = static fn (array $x): float => $x[0] * $x[0];
    expect(NelderMead::new($simplex)->withSdTolerance(NAN)->run($cost, 5)->termination)
        ->toBe(NelderMeadTermination::MaxItersReached)
        ->and(NelderMead::new($simplex)->run($cost)->termination)->toBe(NelderMeadTermination::SolverConverged);
});

it('maps the lmm.rs call site', function () {
    // lmm.rs: start, then start + 0.1 in each coordinate; the budget error ends the run.
    $start = [0.5, 0.0];
    $simplex = [$start, [0.6, 0.0], [0.5, 0.1]];
    $calls = 0;
    $cost = function (array $theta) use (&$calls): float {
        if ($calls >= 10) {
            throw new NelderMeadError('Nelder-Mead evaluation budget exhausted');
        }
        $calls++;

        return $theta[0] < 0.0 ? INF : ($theta[0] - 0.3) ** 2 + ($theta[1] + 0.2) ** 2;
    };
    try {
        NelderMead::new($simplex)->withSdTolerance(1e-10)->run($cost, 10);
        $message = null;
    } catch (NelderMeadError $e) {
        $message = 'mixed-model Nelder-Mead failed: '.$e->getMessage();
    }
    expect($message)->toBe('mixed-model Nelder-Mead failed: Nelder-Mead evaluation budget exhausted');

    $calls = 0;
    $tight = NelderMead::new([$start, [0.6, 0.0], [0.5, 0.1]])->withSdTolerance(1e-10);
    expect(fn () => $tight->run(function (array $x) use (&$calls): float {
        if (++$calls > 2) {
            throw new NelderMeadError('Nelder-Mead evaluation budget exhausted');
        }

        return 0.0;
    }))->toThrow(RustPanic::class, 'Rust panic: called `Result::unwrap()` on an `Err` value: Nelder-Mead evaluation budget exhausted');
});
