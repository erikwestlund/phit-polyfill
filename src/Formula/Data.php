<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Formula;

use PhpStats\Polyfill\ErrorCode;

/**
 * Named data columns for a formula, converted from PHP values as the
 * extension's formula glue does (crates/php-extension/src/formula.rs).
 *
 * Columns use the evaluator's value shapes (num, lgl, text, factor) plus
 * ['t' => 'invalid', 'len' => int, 'message' => string] for columns that are
 * an error only when a formula references them.
 *
 * @internal
 */
final class Data
{
    /**
     * @param  list<string>  $names
     * @param  list<array>  $columns
     */
    private function __construct(
        public readonly array $names,
        private readonly array $columns,
        public readonly int $rows,
    ) {}

    public function column(string $name): ?array
    {
        $i = \array_search($name, $this->names, true);

        return $i === false ? null : $this->columns[$i];
    }

    private static function invalid(string $message): \Exception
    {
        return ErrorCode::exception(ErrorCode::INVALID_INPUT, $message);
    }

    /** Whether ext-php-rs reads $value as a Rust &str (a valid UTF-8 string). */
    public static function utf8(mixed $value): bool
    {
        return \is_string($value) && \mb_check_encoding($value, 'UTF-8');
    }

    /** An array of columns, or an object with toModelData() or toColumns(). */
    public static function from(mixed $data): self
    {
        if (\is_array($data)) {
            return self::fromColumns($data);
        }
        if (\is_object($data)) {
            foreach (['toModelData', 'toColumns'] as $method) {
                if (! \method_exists($data, $method)) {
                    continue;
                }
                $result = new \ReflectionMethod($data, $method)->invoke($data);
                if (! \is_array($result)) {
                    throw self::invalid("$method() must return an array of columns");
                }

                return self::fromColumns($result);
            }
        }
        throw self::invalid('data must be an array of columns or an object with toModelData() or toColumns()');
    }

    private static function fromColumns(array $table): self
    {
        $columns = [];
        foreach ($table as $key => $value) {
            if (! \is_string($key)) {
                throw self::invalid('data must be an associative array of columns (name => list of values), not a list of rows');
            }
            $columns[] = [$key, self::convert($key, $value)];
        }
        $n = $columns === [] ? 0 : self::length($columns[0][1]);
        $names = [];
        $values = [];
        foreach ($columns as [$name, $column]) {
            $length = self::length($column);
            if ($length !== $n) {
                throw ErrorCode::exception(ErrorCode::DIMENSION_MISMATCH, "column '$name' has $length values; every column needs $n");
            }
            if ($column['t'] === 'factor' && count(\array_unique($column['levels'], SORT_STRING)) !== count($column['levels'])) {
                throw self::invalid("factor column '$name' has duplicate levels");
            }
            $names[] = $name;
            $values[] = $column;
        }

        return new self($names, $values, $n);
    }

    private static function length(array $column): int
    {
        return match ($column['t']) {
            'factor' => count($column['codes']),
            'invalid' => $column['len'],
            default => count($column['v']),
        };
    }

    private static function convert(string $name, mixed $value): array
    {
        if (! \is_array($value)) {
            throw self::invalid("column '$name' must be an array of values or a factor descriptor");
        }
        foreach ($value as $key => $_) {
            if (! \is_int($key)) {
                return self::descriptor($name, $value);
            }
        }
        $numbers = $bools = $strings = $nested = false;
        foreach ($value as $key => $v) {
            if (\ReflectionReference::fromArrayElement($value, $key) !== null) {
                // ext-php-rs classifies the raw slot, so a PHP reference counts as a nested value.
                $nested = true;

                continue;
            }
            if ($v === null) {
                continue;
            }
            if (\is_int($v) || \is_float($v)) {
                if (\is_float($v) && \is_infinite($v)) {
                    throw ErrorCode::exception(ErrorCode::NONFINITE, "column '$name' contains an infinite value");
                }
                $numbers = true;
            } elseif (\is_bool($v)) {
                $bools = true;
            } elseif (\is_string($v)) {
                $strings = true;
            } else {
                $nested = true;
            }
        }
        $len = count($value);
        if ($nested) {
            return ['t' => 'invalid', 'len' => $len, 'message' => 'list-columns (arrays or objects as values) cannot be used in a model'];
        }
        if ((int) $numbers + (int) $bools + (int) $strings > 1) {
            return ['t' => 'invalid', 'len' => $len, 'message' => 'values mix numbers, booleans and strings'];
        }
        $out = [];
        if ($strings) {
            foreach ($value as $v) {
                $out[] = self::utf8($v) ? $v : null;
            }

            return ['t' => 'text', 'v' => $out];
        }
        if ($bools) {
            foreach ($value as $v) {
                $out[] = \is_bool($v) ? $v : null;
            }

            return ['t' => 'lgl', 'v' => $out];
        }
        foreach ($value as $v) {
            $out[] = ($v === null || (\is_float($v) && \is_nan($v))) ? null : (float) $v;
        }

        return ['t' => 'num', 'v' => $out];
    }

    private static function descriptor(string $name, array $table): array
    {
        $malformed = fn (string $detail) => self::invalid(
            "column '$name' is not a valid factor descriptor (exactly levels: list<string>, codes: list<?int>, ordered: bool): $detail"
        );
        if (count($table) !== 3) {
            throw $malformed('unexpected keys');
        }
        if (! \is_array($table['levels'] ?? null)) {
            throw $malformed('levels must be an array');
        }
        $levels = [];
        foreach ($table['levels'] as $level) {
            if (! self::utf8($level)) {
                throw $malformed('levels must be strings');
            }
            $levels[] = $level;
        }
        if (! \is_array($table['codes'] ?? null)) {
            throw $malformed('codes must be an array');
        }
        $codes = [];
        foreach ($table['codes'] as $code) {
            if ($code === null) {
                $codes[] = null;
            } elseif (\is_int($code)) {
                if ($code < 0 || $code >= count($levels)) {
                    throw $malformed('codes must be 0-based indexes into levels');
                }
                $codes[] = $code;
            } else {
                throw $malformed('codes must be integers or null');
            }
        }
        if (! \is_bool($table['ordered'] ?? null)) {
            throw $malformed('ordered must be a boolean');
        }

        return ['t' => 'factor', 'levels' => $levels, 'codes' => $codes, 'ordered' => $table['ordered']];
    }
}
