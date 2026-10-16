<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Formula;

use PhpStats\Polyfill\ErrorCode;

/**
 * Evaluates formula variables over all input rows, before missing-value
 * deletion, as R's model.frame does (stats-core formula::eval). NaN results
 * are missing.
 *
 * Values are arrays keyed by 't':
 *   num    ['v' => list<?float>]
 *   mat    ['cols' => list<list<?float>>, 'names' => list<string>]
 *   lgl    ['v' => list<?bool>]
 *   text   ['v' => list<?string>]
 *   factor ['levels' => list<string>, 'codes' => list<?int>, 'ordered' => bool]
 *   str    ['s' => string]
 *
 * @internal
 */
final class Evaluator
{
    public static function missing(array $value, int $i): bool
    {
        return match ($value['t']) {
            'num', 'lgl', 'text' => $value['v'][$i] === null,
            'factor' => $value['codes'][$i] === null,
            'mat' => (function () use ($value, $i) {
                foreach ($value['cols'] as $c) {
                    if ($c[$i] === null) {
                        return true;
                    }
                }

                return false;
            })(),
            default => false,
        };
    }

    private static function clean(float $v): ?float
    {
        return is_nan($v) ? null : $v;
    }

    /** @return list<?float> */
    private static function numeric(array $value, string $what): array
    {
        if ($value['t'] === 'num') {
            return $value['v'];
        }
        if ($value['t'] === 'lgl') {
            return array_map(fn (?bool $b) => $b === null ? null : ($b ? 1.0 : 0.0), $value['v']);
        }
        throw Parser::error("'$what' requires numeric or logical values, not a factor or string");
    }

    private static function map1(array $v, callable $f): array
    {
        $out = [];
        foreach ($v as $x) {
            $out[] = $x === null ? null : self::clean($f($x));
        }

        return ['t' => 'num', 'v' => $out];
    }

    /** LLVM's powi (compiler-rt __powidf2), which Rust's f64::powi calls. */
    public static function powi(float $a, int $b): float
    {
        $recip = $b < 0;
        $r = 1.0;
        while (true) {
            if ($b & 1) {
                $r *= $a;
            }
            $b = intdiv($b, 2);
            if ($b === 0) {
                break;
            }
            $a *= $a;
        }

        return $recip ? 1.0 / $r : $r;
    }

    /** Sorted unique labels in C-locale (byte) order. @return list<string> */
    public static function sortedLevels(array $values): array
    {
        $levels = array_values(array_unique($values, SORT_STRING));
        sort($levels, SORT_STRING);

        return $levels;
    }

    /** Rust total_cmp order for non-NaN floats: -0 before +0. */
    private static function totalCompare(float $a, float $b): int
    {
        if ($a < $b) {
            return -1;
        }
        if ($a > $b) {
            return 1;
        }
        $sa = $a === 0.0 && fdiv(1.0, $a) < 0;
        $sb = $b === 0.0 && fdiv(1.0, $b) < 0;

        return $sb <=> $sa;
    }

    /** Sort with total_cmp, then dedup equal neighbours (keeping the first). @return list<float> */
    private static function distinct(array $values): array
    {
        usort($values, self::totalCompare(...));
        $out = [];
        foreach ($values as $v) {
            if ($out === [] || $out[count($out) - 1] != $v) {
                $out[] = $v;
            }
        }

        return $out;
    }

    private static function toFactor(array $value, bool $ordered): array
    {
        switch ($value['t']) {
            case 'num':
                $present = array_values(array_filter($value['v'], fn ($x) => $x !== null));
                $distinct = self::distinct($present);
                $levels = array_map(Parser::formatNumber(...), $distinct);
                $codes = [];
                foreach ($value['v'] as $x) {
                    $codes[] = $x === null ? null : self::partitionPoint($distinct, $x);
                }

                return ['t' => 'factor', 'levels' => $levels, 'codes' => $codes, 'ordered' => $ordered];
            case 'lgl':
                $levels = [];
                if (in_array(false, $value['v'], true)) {
                    $levels[] = 'FALSE';
                }
                if (in_array(true, $value['v'], true)) {
                    $levels[] = 'TRUE';
                }
                $codes = [];
                foreach ($value['v'] as $b) {
                    $codes[] = $b === null ? null : (int) (array_search($b ? 'TRUE' : 'FALSE', $levels, true) ?: 0);
                }

                return ['t' => 'factor', 'levels' => $levels, 'codes' => $codes, 'ordered' => $ordered];
            case 'text':
                $levels = self::sortedLevels(array_values(array_filter($value['v'], fn ($x) => $x !== null)));
                $index = array_flip($levels);
                $codes = [];
                foreach ($value['v'] as $x) {
                    $codes[] = $x === null ? null : ($index[$x] ?? null);
                }

                return ['t' => 'factor', 'levels' => $levels, 'codes' => $codes, 'ordered' => $ordered];
            case 'factor':
                // factor(f) drops unused levels and keeps the remaining order.
                [$levels, $codes] = self::dropUnused($value['levels'], $value['codes']);

                return ['t' => 'factor', 'levels' => $levels, 'codes' => $codes, 'ordered' => $ordered];
            default:
                throw Parser::error('factor() requires a single column');
        }
    }

    /** Number of sorted values strictly below $x. */
    private static function partitionPoint(array $sorted, float $x): int
    {
        $lo = 0;
        $hi = count($sorted);
        while ($lo < $hi) {
            $mid = intdiv($lo + $hi, 2);
            if ($sorted[$mid] < $x) {
                $lo = $mid + 1;
            } else {
                $hi = $mid;
            }
        }

        return $lo;
    }

    /** @return array{0: list<string>, 1: list<?int>} */
    public static function dropUnused(array $levels, array $codes): array
    {
        $used = array_fill(0, count($levels), false);
        foreach ($codes as $c) {
            if ($c !== null) {
                $used[$c] = true;
            }
        }
        $map = [];
        $kept = [];
        foreach ($levels as $i => $level) {
            if ($used[$i]) {
                $map[$i] = count($kept);
                $kept[] = $level;
            }
        }
        $out = [];
        foreach ($codes as $c) {
            $out[] = $c === null ? null : $map[$c];
        }

        return [$kept, $out];
    }

    /**
     * Orthonormal polynomial basis (degree columns, constant removed) of
     * centered x, R's poly()/contr.poly construction by two Gram-Schmidt passes.
     *
     * @param  list<float>  $x
     * @return list<list<float>>
     */
    public static function orthonormalPoly(array $x, int $degree): array
    {
        if ($degree >= count(self::distinct($x))) {
            throw Parser::error("'degree' must be less than number of unique points");
        }
        $sum = 0.0;
        foreach ($x as $v) {
            $sum += $v;
        }
        $mean = $sum / count($x);
        $centered = array_map(fn ($v) => $v - $mean, $x);
        $basis = [];
        for ($k = 0; $k <= $degree; $k++) {
            $v = array_map(fn ($c) => self::powi($c, $k), $centered);
            for ($pass = 0; $pass < 2; $pass++) {
                foreach ($basis as $q) {
                    $dot = 0.0;
                    foreach ($q as $i => $qi) {
                        $dot += $qi * $v[$i];
                    }
                    foreach ($q as $i => $qi) {
                        $v[$i] -= $dot * $qi;
                    }
                }
            }
            $norm = 0.0;
            foreach ($v as $a) {
                $norm += $a * $a;
            }
            $norm = sqrt($norm);
            if ($norm == 0.0 || ! is_finite($norm)) {
                throw ErrorCode::exception(ErrorCode::RANK_DEFICIENCY, 'polynomial basis is numerically rank deficient');
            }
            $basis[] = array_map(fn ($a) => $a / $norm, $v);
        }
        array_shift($basis);

        return $basis;
    }

    /** @param list<array{0: ?string, 1: array}> $args */
    private static function arg(array $args, int $position, string $name): ?array
    {
        foreach ($args as [$n, $value]) {
            if ($n === $name) {
                return $value;
            }
        }
        $k = 0;
        foreach ($args as [$n, $value]) {
            if ($n === null) {
                if ($k === $position) {
                    return $value;
                }
                $k++;
            }
        }

        return null;
    }

    private static function constant(array $expr): ?float
    {
        return match (true) {
            $expr[0] === 'num' => $expr[1],
            $expr[0] === 'unary' && $expr[1] === '-' => ($c = self::constant($expr[2])) === null ? null : -$c,
            $expr[0] === 'paren' => self::constant($expr[1]),
            default => null,
        };
    }

    public static function evaluate(array $expr, Data $data): array
    {
        $n = $data->rows;
        switch ($expr[0]) {
            case 'num':
                return ['t' => 'num', 'v' => array_fill(0, $n, $expr[1])];
            case 'bool':
                return ['t' => 'lgl', 'v' => array_fill(0, $n, $expr[1])];
            case 'str':
                return ['t' => 'str', 's' => $expr[1]];
            case 'dot':
                throw Parser::error("'.' is only valid as a formula term");
            case 'paren':
                return self::evaluate($expr[1], $data);
            case 'name':
                $column = $data->column($expr[1]) ?? throw Parser::error("unknown column '{$expr[1]}' in formula");
                if ($column['t'] === 'invalid') {
                    throw Parser::error("column '{$expr[1]}': {$column['message']}");
                }

                return $column;
            case 'unary':
                $v = self::evaluate($expr[2], $data);
                switch ($expr[1]) {
                    case '-':
                        return self::map1(self::numeric($v, Parser::deparse($expr)), fn ($x) => -$x);
                    case '+':
                        return ['t' => 'num', 'v' => self::numeric($v, Parser::deparse($expr))];
                    case '!':
                        if ($v['t'] === 'lgl') {
                            return ['t' => 'lgl', 'v' => array_map(fn ($b) => $b === null ? null : ! $b, $v['v'])];
                        }
                        if ($v['t'] === 'num') {
                            return ['t' => 'lgl', 'v' => array_map(fn ($x) => $x === null ? null : $x == 0.0, $v['v'])];
                        }
                        throw Parser::error("'!' requires logical values");
                    default:
                        throw Parser::error("cannot evaluate '".Parser::deparse($expr)."'");
                }
            case 'binary':
                $left = self::evaluate($expr[2], $data);
                $right = self::evaluate($expr[3], $data);

                return self::binary($expr[1], $left, $right, $expr);
            default:
                return self::call($expr[1], $expr[2], $data, $expr);
        }
    }

    private static function labelsOf(array $value): ?array
    {
        if ($value['t'] === 'text') {
            return $value['v'];
        }
        if ($value['t'] === 'factor') {
            return array_map(fn ($c) => $c === null ? null : $value['levels'][$c], $value['codes']);
        }

        return null;
    }

    private static function binary(string $op, array $left, array $right, array $expr): array
    {
        $text = Parser::deparse($expr);
        if ($op === '==' || $op === '!=') {
            $pair = null;
            if ($left['t'] === 'str' || $right['t'] === 'str') {
                [$s, $other] = $left['t'] === 'str' ? [$left['s'], $right] : [$right['s'], $left];
                $labels = self::labelsOf($other) ?? throw Parser::error("'$text' compares a string with a non-string");
                $pair = [$labels, $s];
            }
            if ($pair !== null) {
                $eq = $op === '==';

                return ['t' => 'lgl', 'v' => array_map(fn ($v) => $v === null ? null : (($v === $pair[1]) === $eq), $pair[0])];
            }
        }
        if ($op === '&' || $op === '|') {
            $a = self::numeric($left, $text);
            $b = self::numeric($right, $text);
            $out = [];
            foreach ($a as $i => $x) {
                $y = $b[$i];
                $x = $x === null ? null : $x != 0.0;
                $y = $y === null ? null : $y != 0.0;
                if ($op === '&') {
                    $out[] = match (true) {
                        $x === false || $y === false => false,
                        $x === true && $y === true => true,
                        default => null,
                    };
                } else {
                    $out[] = match (true) {
                        $x === true || $y === true => true,
                        $x === false && $y === false => false,
                        default => null,
                    };
                }
            }

            return ['t' => 'lgl', 'v' => $out];
        }
        $a = self::numeric($left, $text);
        $b = self::numeric($right, $text);
        $f = match ($op) {
            '+' => fn ($x, $y) => $x + $y,
            '-' => fn ($x, $y) => $x - $y,
            '*' => fn ($x, $y) => $x * $y,
            '/' => fn ($x, $y) => fdiv($x, $y),
            // R: 1^y and x^0 are 1 even for NaN/Inf operands.
            '^' => fn ($x, $y) => ($x == 1.0 || $y == 0.0) ? 1.0 : $x ** $y,
            '<', '>', '<=', '>=', '==', '!=' => null,
            default => throw Parser::error("operator '$text' is not supported inside I() or a transform"),
        };
        $out = [];
        if ($f !== null) {
            foreach ($a as $i => $x) {
                $y = $b[$i];
                $out[] = ($x === null || $y === null) ? null : self::clean((float) $f($x, $y));
            }

            return ['t' => 'num', 'v' => $out];
        }
        foreach ($a as $i => $x) {
            $y = $b[$i];
            $out[] = ($x === null || $y === null) ? null : match ($op) {
                '<' => $x < $y,
                '>' => $x > $y,
                '<=' => $x <= $y,
                '>=' => $x >= $y,
                '==' => $x == $y,
                '!=' => $x != $y,
            };
        }

        return ['t' => 'lgl', 'v' => $out];
    }

    private static function call(string $name, array $args, Data $data, array $expr): array
    {
        $text = Parser::deparse($expr);
        $first = fn () => self::arg($args, 0, 'x') ?? throw Parser::error("$name() needs an argument");
        $unary = function (callable $f) use ($args, $name, $first, $data, $text): array {
            if (count($args) !== 1) {
                throw Parser::error("$name() takes one argument");
            }

            return self::map1(self::numeric(self::evaluate($first(), $data), $text), $f);
        };
        switch ($name) {
            case 'I':
            case 'offset':
                if (count($args) !== 1) {
                    throw Parser::error("$name() takes one argument");
                }
                $v = self::evaluate($first(), $data);
                if ($name === 'offset') {
                    return ['t' => 'num', 'v' => self::numeric($v, $text)];
                }

                return $v;
            case 'log':
                $x = self::numeric(self::evaluate($first(), $data), $text);
                $base = self::arg($args, 1, 'base');
                if ($base === null && count($args) === 1) {
                    return self::map1($x, fn ($v) => log($v));
                }
                if ($base !== null && count($args) === 2) {
                    $b = self::constant($base) ?? throw Parser::error('log() base must be a number');
                    $lb = $b > 0.0 ? log($b) : ($b == 0.0 ? -INF : NAN);

                    return self::map1($x, fn ($v) => fdiv(log($v), $lb));
                }
                throw Parser::error('log() takes x and an optional base');
            case 'log2':
                return $unary(fn ($v) => $v > 0.0 ? log($v, 2) : log($v));
            case 'log10':
                return $unary(log10(...));
            case 'log1p':
                return $unary(log1p(...));
            case 'exp':
                return $unary(exp(...));
            case 'expm1':
                return $unary(expm1(...));
            case 'sqrt':
                return $unary(sqrt(...));
            case 'abs':
                return $unary(fn ($v) => abs($v));
            case 'sin':
                return $unary(sin(...));
            case 'cos':
                return $unary(cos(...));
            case 'tan':
                return $unary(tan(...));
            case 'as.numeric':
            case 'as.double':
                if (count($args) !== 1) {
                    throw Parser::error("$name() takes one argument");
                }
                $v = self::evaluate($first(), $data);
                if ($v['t'] === 'factor') {
                    return ['t' => 'num', 'v' => array_map(fn ($c) => $c === null ? null : $c + 1.0, $v['codes'])];
                }

                return ['t' => 'num', 'v' => self::numeric($v, $text)];
            case 'factor':
            case 'as.factor':
            case 'ordered':
                if (count($args) !== 1) {
                    throw Parser::error("$name() supports one argument; set levels in a factor descriptor");
                }

                return self::toFactor(self::evaluate($first(), $data), $name === 'ordered');
            case 'poly':
                return self::poly($args, $data, $text, $first);
            case 'cbind':
                $cols = [];
                $names = [];
                foreach ($args as [$argName, $value]) {
                    $cols[] = self::numeric(self::evaluate($value, $data), $text);
                    $names[] = $argName ?? Parser::deparse($value);
                }
                if ($cols === []) {
                    throw Parser::error('cbind() needs arguments');
                }

                return ['t' => 'mat', 'cols' => $cols, 'names' => $names];
            default:
                throw ErrorCode::exception(ErrorCode::UNSUPPORTED, "unsupported function '$name()' in formula");
        }
    }

    private static function poly(array $args, Data $data, string $text, \Closure $first): array
    {
        $x = self::numeric(self::evaluate($first(), $data), $text);
        $degreeExpr = null;
        foreach ($args as [$n, $value]) {
            if ($n === 'degree') {
                $degreeExpr = $value;
                break;
            }
        }
        if ($degreeExpr === null) {
            $k = 0;
            foreach ($args as [$n, $value]) {
                if ($n === null && $k++ === 1) {
                    $degreeExpr = $value;
                    break;
                }
            }
        }
        $degree = $degreeExpr === null ? 1.0
            : (self::constant($degreeExpr) ?? throw Parser::error('poly() degree must be a number'));
        if ($degree < 1.0 || fmod($degree, 1.0) != 0.0 || $degree > 25.0) {
            throw Parser::error('poly() degree must be an integer from 1 to 25');
        }
        $degree = (int) $degree;
        $unnamed = 0;
        foreach ($args as [$n]) {
            if ($n !== null && ! in_array($n, ['x', 'degree', 'raw'], true)) {
                throw Parser::error("unsupported poly() argument '$n'");
            }
        }
        foreach ($args as [$n]) {
            $unnamed += $n === null ? 1 : 0;
        }
        if ($unnamed > 2) {
            throw Parser::error('poly() of several variables is not supported');
        }
        $rawExpr = null;
        foreach ($args as [$n, $value]) {
            if ($n === 'raw') {
                $rawExpr = $value;
                break;
            }
        }
        if ($rawExpr !== null && $rawExpr[0] !== 'bool') {
            throw Parser::error('poly() raw must be TRUE or FALSE');
        }
        $raw = $rawExpr !== null && $rawExpr[1];
        $names = array_map('strval', range(1, $degree));
        if ($raw) {
            $cols = [];
            for ($k = 1; $k <= $degree; $k++) {
                $cols[] = array_map(fn ($v) => $v === null ? null : self::powi($v, $k), $x);
            }

            return ['t' => 'mat', 'cols' => $cols, 'names' => $names];
        }
        if (in_array(null, $x, true)) {
            throw Parser::error('missing values are not allowed in orthogonal poly(); use raw = TRUE or remove them');
        }

        return ['t' => 'mat', 'cols' => self::orthonormalPoly($x, $degree), 'names' => $names];
    }
}
