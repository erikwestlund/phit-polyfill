<?php

declare(strict_types=1);

// Accuracy stress tests beyond the fixtures: tail probabilities down to ~1e-300
// and log-scale tails far beyond the double range, against closed forms or
// independently summed series; extreme parameters; quantile/CDF round trips.
// Tolerances follow the extension suite: 1e-12 relative for probabilities,
// 1e-12 relative to max(1, |value|) on the log scale.

if (! function_exists('stressRel')) {
    function stressRel(float $actual, float $expected): float
    {
        if ($actual === $expected) {
            return 0.0;
        }

        return abs($actual - $expected) / abs($expected);
    }

    function stressProb(string $fn, array $args, float $expected, float $rel = 1e-12): void
    {
        $actual = pf($fn, ...$args);
        expect(stressRel($actual, $expected) <= $rel)->toBeTrue(sprintf(
            '%s(%s) = %.17g, expected %.17g (rel %.3g)', $fn, json_encode($args), $actual, $expected, stressRel($actual, $expected),
        ));
    }

    function stressLog(string $fn, array $args, float $expected, float $rel = 1e-12): void
    {
        $actual = pf($fn, ...$args);
        $error = abs($actual - $expected) / max(1.0, abs($expected));
        expect($error <= $rel)->toBeTrue(sprintf(
            '%s(%s) = %.17g, expected %.17g (err %.3g)', $fn, json_encode($args), $actual, $expected, $error,
        ));
    }

    /** ln n! - ((n + 1/2) ln n - n + ln √(2π)) by the Stirling series (n ≥ 100). */
    function stressStirlingCorrection(float $n): float
    {
        $n2 = $n * $n;

        return 1 / (12 * $n) - 1 / (360 * $n * $n2) + 1 / (1260 * $n * $n2 * $n2) - 1 / (1680 * $n * $n2 ** 3);
    }

    /** ln n! exactly summed for small n, Stirling with five correction terms beyond. */
    function stressLnFactorial(int $n): float
    {
        if ($n < 100) {
            $s = 0.0;
            for ($i = 2; $i <= $n; $i++) {
                $s += log($i);
            }

            return $s;
        }
        $x = $n + 1.0;
        $x2 = $x * $x;

        return ($x - 0.5) * log($x) - $x + 0.5 * log(2 * M_PI)
            + 1 / (12 * $x) - 1 / (360 * $x * $x2) + 1 / (1260 * $x * $x2 * $x2) - 1 / (1680 * $x * $x2 ** 3);
    }

    /** ln Q(k, x) = ln(e^-x Σ_{j<k} x^j/j!) for integer k and x > k (terms summed downward). */
    function stressLnUpperGamma(int $k, float $x): float
    {
        $sum = 0.0;
        $term = 1.0;
        for ($m = 0; $m < $k; $m++) {
            $sum += $term;
            $term *= ($k - 1 - $m) / $x;
            if ($term < 1e-17 * $sum) {
                break;
            }
        }

        // -x + n ln x - ln n! with n = k - 1, written as -(x - n) + n log1p((x - n)/n) - ...
        // so the ~n ln n terms cancel analytically, not numerically.
        $n = $k - 1;
        if ($n < 100) {
            return -$x + $n * log($x) - stressLnFactorial($n) + log($sum);
        }

        return -($x - $n) + $n * log1p(($x - $n) / $n) - 0.5 * log(2 * M_PI * $n) - stressStirlingCorrection($n) + log($sum);
    }

    /** ln P(k, x) = ln(x^k e^-x/k! Σ_j x^j/((k+1)...(k+j))) for integer k and x < k. */
    function stressLnLowerGamma(int $k, float $x): float
    {
        $sum = 0.0;
        $term = 1.0;
        for ($j = 1; $j < 100000 && $term >= 1e-17 * $sum; $j++) {
            $sum += $term;
            $term *= $x / ($k + $j);
        }

        if ($k < 100) {
            return $k * log($x) - $x - stressLnFactorial($k) + log($sum);
        }

        return $k * log1p(($x - $k) / $k) - ($x - $k) - 0.5 * log(2 * M_PI * $k) - stressStirlingCorrection($k) + log($sum);
    }

    /** ln of the standard normal upper tail for z ≥ 20 from the asymptotic Mills-ratio series. */
    function stressLnNormalTail(float $z): float
    {
        $s = 1.0;
        $term = 1.0;
        $inv = 1.0 / ($z * $z);
        for ($n = 1; $n <= 12; $n++) {
            $term *= -(2 * $n - 1) * $inv;
            $s += $term;
        }

        return -0.5 * $z * $z - 0.5 * log(2 * M_PI) - log($z) + log($s);
    }
}

it('matches the normal tail asymptotics down to 1e-300 and far beyond on the log scale', function () {
    foreach ([20.0, 25.0, 30.0, 35.0, 37.0, 37.5] as $z) {
        $ln = stressLnNormalTail($z);
        stressProb('stats_normal_sf', [$z], exp($ln));
        stressProb('stats_normal_cdf', [-$z], exp($ln));
        $x = 5.0 - 1e-8 * $z; // the input itself rounds, so use the z that x - mean actually forms
        stressProb('stats_normal_cdf', [$x, 5.0, 1e-8], exp(stressLnNormalTail((5.0 - $x) / 1e-8)));
        stressLog('stats_normal_logsf', [$z], $ln);
    }
    foreach ([38.5, 40.0, 100.0, 1e3, 1e10, 1e100, 1e150] as $z) {
        stressLog('stats_normal_logsf', [$z], stressLnNormalTail($z));
        stressLog('stats_normal_logcdf', [-$z], stressLnNormalTail($z));
        stressLog('stats_normal_logsf', [3.0 + 4.0 * $z, 3.0, 4.0], stressLnNormalTail($z));
    }
    // Log-normal tails are normal tails in ln x.
    foreach ([30.0, 40.0, 300.0] as $z) {
        $x = exp(2.0 + 0.5 * $z);
        $zz = (log($x) - 2.0) / 0.5;
        stressLog('stats_lognormal_logsf', [$x, 2.0, 0.5], stressLnNormalTail($zz));
    }
});

it('matches closed-form Student t tails for df = 1 and df = 2', function () {
    foreach ([10.0, 1e5, 1e50, 1e150, 1e300] as $x) {
        // Cauchy: S(x) = atan(1/x)/π.
        stressProb('stats_t_sf', [$x, 1.0], atan2(1.0, $x) / M_PI);
        stressProb('stats_t_cdf', [-$x, 1.0], atan2(1.0, $x) / M_PI);
        stressLog('stats_t_logsf', [$x, 1.0], log(atan2(1.0, $x)) - log(M_PI));
    }
    foreach ([10.0, 1e5, 1e50, 1e100, 7e149] as $x) {
        // df = 2: S(x) = 1/(s(s + x)), s = √(2 + x²).
        $s = sqrt(2.0 + $x * $x);
        stressProb('stats_t_sf', [$x, 2.0], 1.0 / ($s * ($s + $x)));
        stressProb('stats_t_cdf', [-$x, 2.0], 1.0 / ($s * ($s + $x)));
    }
    foreach ([1e10, 1e200, 1e300] as $x) {
        stressLog('stats_t_logsf', [$x, 2.0], -M_LN2 - 2.0 * log($x));
        stressLog('stats_t_logcdf', [-$x, 2.0], -M_LN2 - 2.0 * log($x));
    }
});

it('approaches the normal distribution for very large df', function () {
    foreach ([-5.0, -2.0, -0.5, 0.5, 3.0] as $x) {
        stressProb('stats_t_cdf', [$x, 1e15], pf('stats_normal_cdf', $x), 1e-12);
        stressProb('stats_t_pdf', [$x, 1e15], pf('stats_normal_pdf', $x), 1e-12);
    }
});

it('matches closed-form chi-square, exponential and gamma tails', function () {
    foreach ([1e-300, 1e-100, 1e-10, 0.5, 30.0, 600.0, 1370.0] as $x) {
        // df = 2 is the exponential with rate 1/2.
        stressProb('stats_chisquare_cdf', [$x, 2.0], -expm1(-$x / 2));
        stressProb('stats_chisquare_sf', [$x, 2.0], exp(-$x / 2));
        stressProb('stats_exponential_sf', [$x, 0.5], exp(-$x / 2));
        stressLog('stats_chisquare_logcdf', [$x, 2.0], log(-expm1(-$x / 2)));
    }
    foreach ([1e4, 1e10, 1e300] as $x) {
        stressLog('stats_chisquare_logsf', [$x, 2.0], -$x / 2);
        stressLog('stats_exponential_logsf', [$x, 1e-3], -1e-3 * $x);
        stressLog('stats_gamma_logsf', [$x, 1.0, 4.0], -$x / 4);
    }
    // df = 1 near zero: F(x) = √(2x/π)(1 - x/6 + ...).
    foreach ([1e-300, 1e-200, 1e-30] as $x) {
        stressProb('stats_chisquare_cdf', [$x, 1.0], sqrt(2 * $x / M_PI));
    }
});

it('matches exactly summed incomplete gamma tails for integer shapes (gamma, chi-square, Poisson)', function () {
    foreach ([[3, 500.0], [3, 700.0], [10, 750.0], [50, 1e3], [1000, 3000.0], [100000, 120000.0], [1000000, 1010000.0]] as [$k, $x]) {
        $ln = stressLnUpperGamma($k, $x);
        stressLog('stats_gamma_logsf', [$x * 8.0, (float) $k, 8.0], $ln);
        stressLog('stats_poisson_logcdf', [$k - 1, $x], $ln);
        if ($ln > log(1e-300)) {
            stressProb('stats_gamma_sf', [$x, (float) $k], exp($ln), 1e-11);
            stressProb('stats_poisson_cdf', [$k - 1, $x], exp($ln), 1e-11);
        }
    }
    foreach ([[3, 1e-100], [10, 1e-30], [10, 1e-31], [50, 1e-5], [200, 1.0], [1000, 500.0], [100000, 90000.0]] as [$k, $x]) {
        $ln = stressLnLowerGamma($k, $x);
        stressLog('stats_gamma_logcdf', [$x, (float) $k], $ln);
        stressLog('stats_poisson_logsf', [$k - 1, $x], $ln);
        stressLog('stats_chisquare_logcdf', [2 * $x, 2.0 * $k], $ln);
        if ($ln > log(1e-300)) {
            stressProb('stats_gamma_cdf', [$x, (float) $k], exp($ln), 1e-11);
            stressProb('stats_poisson_sf', [$k - 1, $x], exp($ln), 1e-11);
        }
    }
});

it('matches closed-form beta and F tails, including tiny and huge shapes', function () {
    foreach ([0.5, 3.0, 1e3, 1e6] as $b) {
        foreach ([1e-300, 1e-10, 0.3, 0.999] as $x) {
            // Beta(1, b): S(x) = (1 - x)^b.
            $ln = $b * log1p(-$x);
            stressLog('stats_beta_logsf', [$x, 1.0, $b], $ln);
            if ($ln > log(1e-300)) {
                stressProb('stats_beta_sf', [$x, 1.0, $b], exp($ln));
            }
            // Beta(b, 1): F(x) = x^b.
            $ln = $b * log($x);
            stressLog('stats_beta_logcdf', [$x, $b, 1.0], $ln);
            if ($ln > log(1e-300)) {
                stressProb('stats_beta_cdf', [$x, $b, 1.0], exp($ln));
            }
        }
    }
    foreach ([0.5, 3.0, 40.0, 1e6] as $d) {
        foreach ([1e-250, 1e-3, 2.0, 1e5, 1e100] as $x) {
            // F(2, d): S(x) = (1 + 2x/d)^(-d/2).
            $ln = -($d / 2) * log1p(2 * $x / $d);
            stressLog('stats_f_logsf', [$x, 2.0, $d], $ln);
            if ($ln > log(1e-300)) {
                stressProb('stats_f_sf', [$x, 2.0, $d], exp($ln));
            }
            // F(d, 2): F(x) = (1 + 2/(d x))^(-d/2).
            $ln = -($d / 2) * log1p(2 / ($d * $x));
            stressLog('stats_f_logcdf', [$x, $d, 2.0], $ln);
            if ($ln > log(1e-300)) {
                stressProb('stats_f_cdf', [$x, $d, 2.0], exp($ln));
            }
        }
    }
});

it('matches closed-form discrete tails at the support edges', function () {
    foreach ([[900, 0.5], [1000000, 0.5], [2 ** 53, 1e-10], [1000, 0.999], [10, 1e-31]] as [$n, $p]) {
        $ln = $n * log1p(-$p); // P(X = 0) = P(X <= 0) = q^n
        stressLog('stats_binomial_logcdf', [0, $n, $p], $ln);
        stressLog('stats_binomial_logpmf', [0, $n, $p], $ln);
        $ln = $n * log($p); // P(X = n) = P(X > n - 1) = p^n
        stressLog('stats_binomial_logsf', [$n - 1, $n, $p], $ln);
        if ($ln > log(1e-300)) {
            stressProb('stats_binomial_sf', [$n - 1, $n, $p], exp($ln));
        }
    }
    foreach ([[3.0, 1e-100], [1e6, 0.5], [0.01, 1e-300]] as [$size, $prob]) {
        stressLog('stats_nbinom_logcdf', [0, $size, $prob], $size * log($prob));
    }
    foreach ([[3.0, 1e100], [1e6, 1e7], [1e-3, 1e300]] as [$size, $mu]) {
        stressLog('stats_nbinom_mu_logcdf', [0, $size, $mu], -$size * log1p($mu / $size));
    }
    stressLog('stats_poisson_logpmf', [0, 1e300], -1e300);
});

it('round-trips continuous quantiles through the CDF from 1e-300 to 1 - 2^-52', function (string $name, array $params) {
    $ps = [1e-300, 1e-200, 1e-100, 1e-50, 1e-20, 1e-8, 0.01, 0.3, 0.5, 0.7, 0.99];
    foreach ($ps as $p) {
        $q = pf("stats_{$name}_quantile", $p, ...$params);
        if (is_infinite($q) || ($q <= 0.0 && in_array($name, ['chisquare', 'f', 'gamma', 'beta', 'lognormal', 'exponential'], true))) {
            continue; // the quantile lies beyond the double range
        }
        if ($name === 'beta' && 1.0 - $q <= 4 * PHP_FLOAT_EPSILON) {
            continue; // e.g. Beta(0.05, 0.05) at 0.99 is 1 - 1e-34: not representable below 1
        }
        $back = $p < 1e-280 ? exp(pf("stats_{$name}_logcdf", $q, ...$params)) : pf("stats_{$name}_cdf", $q, ...$params);
        // The extension's log-scale tolerance, 1e-12 max(1, |ln p|), plus the rounding of Q
        // itself: a relative change eps in Q moves F by |Q| f(Q)/F(Q) eps.
        $sensitivity = abs($q) * exp(pf("stats_{$name}_logpdf", $q, ...$params) - log($p));
        $tolerance = 1e-12 * max(1.0, abs(log($p))) + 8 * PHP_FLOAT_EPSILON * $sensitivity;
        expect(abs($back - $p) <= $tolerance * $p)->toBeTrue("$name(".json_encode($params).") F(Q($p)) = $back");
    }
    foreach ([1, 10, 30, 52] as $k) {
        $tail = 2.0 ** -$k;
        $q = pf("stats_{$name}_quantile", 1.0 - $tail, ...$params);
        if (is_infinite($q) || ($name === 'beta' && 1.0 - $q <= 4 * PHP_FLOAT_EPSILON)) {
            continue; // beyond the double range, or not representable below 1
        }
        $back = pf("stats_{$name}_sf", $q, ...$params);
        $sensitivity = abs($q) * exp(pf("stats_{$name}_logpdf", $q, ...$params) - log($tail));
        $tolerance = 1e-12 * max(1.0, abs(log($tail))) + 8 * PHP_FLOAT_EPSILON * $sensitivity;
        expect(abs($back - $tail) <= $tolerance * $tail)->toBeTrue("$name(".json_encode($params).") S(Q(1 - 2^-$k)) = $back");
    }
})->with([
    ['normal', [0.0, 1.0]], ['normal', [1e6, 1e-8]], ['normal', [-2.0, 0.5]], ['t', [0.2]], ['t', [1.0]], ['t', [3.0]], ['t', [1e8]],
    ['chisquare', [0.01]], ['chisquare', [3.0]], ['chisquare', [1e6]], ['f', [0.5, 0.5]], ['f', [5.0, 1e6]],
    ['f', [1e4, 3.0]], ['gamma', [1e-3, 2.0]], ['gamma', [50.0, 1e-5]], ['gamma', [1e5, 1.0]],
    ['beta', [0.05, 0.05]], ['beta', [2.0, 5.0]], ['beta', [1e4, 3.0]], ['beta', [1e5, 1e5]],
    ['exponential', [1e-50]], ['exponential', [1e50]], ['lognormal', [0.0, 3.0]], ['lognormal', [-50.0, 0.1]],
    ['logistic', [0.0, 1.0]], ['logistic', [1e3, 1e-3]],
]);

it('round-trips discrete CDFs through the quantile across large and small parameters', function (string $name, array $params, array $ks) {
    // Q(p) = k for every p in (F(k - 1), F(k)]; probe the midpoint, since at p = F(k) itself
    // the upper-half search compares an independently computed S(k) with 1 - p.
    foreach ($ks as $k) {
        $below = $k > 0 ? pf("stats_{$name}_cdf", $k - 1, ...$params) : 0.0;
        $at = pf("stats_{$name}_cdf", $k, ...$params);
        if ($at - $below < 1e-9 * $at) {
            continue;
        }
        $p = $below + 0.5 * ($at - $below);
        expect(pf("stats_{$name}_quantile", $p, ...$params))->toBe((float) $k, "$name(".json_encode($params).") Q($p)");
    }
})->with([
    ['binomial', [10, 0.3], range(0, 9)],
    ['binomial', [1000000, 1e-4], [0, 50, 100, 150, 200]],
    ['binomial', [2 ** 40, 0.5], [2 ** 39 - 2000000, 2 ** 39, 2 ** 39 + 1000000]],
    ['poisson', [1e-8], [0, 1, 2]],
    ['poisson', [4.5], range(0, 25)],
    ['poisson', [1e9], [1e9 - 1e5, 1e9, 1e9 + 2e5]],
    ['nbinom', [0.01, 0.5], [0, 1, 10, 1000]],
    ['nbinom', [1e6, 0.999], [800, 1000, 1200]],
    ['nbinom_mu', [0.5, 1e6], [0, 10, 1000, 1000000, 10000000]],
    ['nbinom_mu', [1e8, 3.0], [0, 1, 3, 10]],
]);

it('keeps densities accurate for extreme parameters', function () {
    // Gamma(k, 1) density at its mode for huge k: (k-1)^(k-1) e^-(k-1)/Γ(k) ≈ 1/√(2π(k-1)) (1 - 1/(12(k-1))).
    foreach ([1e8, 1e12] as $k) {
        $n = $k - 1;
        stressProb('stats_gamma_pdf', [$n, $k], (1 - 1 / (12 * $n)) / sqrt(2 * M_PI * $n), 1e-12);
    }
    // Poisson pmf at the mean for huge λ: ≈ 1/√(2πλ) (1 - 1/(12λ)).
    foreach ([1e10, 1e14] as $lambda) {
        stressProb('stats_poisson_pmf', [$lambda, $lambda], (1 - 1 / (12 * $lambda)) / sqrt(2 * M_PI * $lambda), 1e-12);
    }
    // Binomial(n, 1/2) at n/2: C(n, n/2) 2^-n ≈ √(2/(πn)) (1 - 1/(4n)).
    foreach ([2 ** 30, 2 ** 52] as $n) {
        stressProb('stats_binomial_pmf', [$n / 2, $n, 0.5], sqrt(2 / (M_PI * $n)) * (1 - 1 / (4 * $n)), 1e-12);
    }
    // Beta(a, a) at 1/2 is 2Γ(a + 1/2)/(Γ(a)√π) ≈ 2√(a/π)(1 - 1/(8a)).
    foreach ([1e6, 1e10] as $a) {
        stressProb('stats_beta_pdf', [0.5, $a, $a], 2.0 * sqrt($a / M_PI) * (1 - 1 / (8 * $a)), 1e-12);
    }
});

it('documents the precision of deep log-scale incomplete gamma tails at huge shapes', function () {
    // Below 1e-280 the extension (and so the polyfill) evaluates ln P(a, x) as
    // a ln x - x - ln Γ(a) + ln(series): absolute error grows like eps·a·ln x,
    // so relative to |ln P| it exceeds 1e-12 once a ≳ 1e6 (5e-10 at a = 1e8).
    // Assert that bound, and the 1e-12 tolerance where it is the tighter one.
    foreach ([[1e5, 9.0e4], [1e6, 9.6e5], [1e7, 9.88e6], [1e8, 9.96e7]] as [$a, $x]) {
        $expected = stressLnLowerGamma((int) $a, $x);
        $actual = pf('stats_gamma_logcdf', $x, $a);
        $bound = max(1e-12 * abs($expected), 4 * PHP_FLOAT_EPSILON * $a * log($x));
        expect(abs($actual - $expected))->toBeLessThanOrEqual($bound);
    }
    expect(abs(pf('stats_gamma_logcdf', 9.0e4, 1e5) - stressLnLowerGamma(100000, 9.0e4)))
        ->toBeLessThanOrEqual(1e-12 * 540.4);
});