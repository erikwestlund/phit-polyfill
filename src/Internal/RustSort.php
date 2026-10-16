<?php

declare(strict_types=1);

namespace PhpStats\Polyfill\Internal;

/**
 * Rust 1.99's `slice::sort_unstable_by` (core::slice::sort::unstable) for
 * Copy elements of 16 bytes, such as nalgebra's (f64, usize) pairs, so that
 * equal keys end up in the same order: insertion sort up to 20 elements,
 * otherwise ipnsort (an existing-run check, then quicksort with
 * pseudo-median pivots, branchless cyclic Lomuto partitions,
 * small_sort_general below 33 elements and a heapsort fallback).
 *
 * @internal
 */
final class RustSort
{
    private const int SMALL_SORT_THRESHOLD = 32;

    /** @var list<mixed> */
    private array $v;

    /** @var \Closure(mixed, mixed): bool */
    private \Closure $isLess;

    /**
     * @template T
     *
     * @param  list<T>  $values
     * @param  callable(T, T): bool  $isLess
     * @return list<T>
     */
    public static function sortUnstableBy(array $values, callable $isLess): array
    {
        $length = count($values);
        if ($length < 2) {
            return $values;
        }
        $sort = new self;
        $sort->v = $values;
        $sort->isLess = $isLess(...);
        if ($length <= 20) {
            $sort->insertionSort(0, $length, 1);
        } else {
            $sort->ipnsort($length);
        }

        return $sort->v;
    }

    private function less(mixed $a, mixed $b): bool
    {
        return ($this->isLess)($a, $b);
    }

    /** insertion_sort_shift_left on v[$base, $base + $length), v[..$offset] sorted. */
    private function insertionSort(int $base, int $length, int $offset): void
    {
        for ($tail = $base + $offset; $tail < $base + $length; $tail++) {
            self::insertTail($this->v, $base, $tail, $this->isLess);
        }
    }

    /** @param  list<mixed>  $v */
    private static function insertTail(array &$v, int $begin, int $tail, \Closure $isLess): void
    {
        $sift = $tail - 1;
        if ($isLess($v[$tail], $v[$sift]) === false) {
            return;
        }
        $tmp = $v[$tail];
        $gap = $tail;
        while (true) {
            $v[$gap] = $v[$sift];
            $gap = $sift;
            if ($sift === $begin) {
                break;
            }
            $sift--;
            if ($isLess($tmp, $v[$sift]) === false) {
                break;
            }
        }
        $v[$gap] = $tmp;
    }

    private function ipnsort(int $length): void
    {
        $descending = $this->less($this->v[1], $this->v[0]);
        $run = 2;
        if ($descending) {
            while ($run < $length && $this->less($this->v[$run], $this->v[$run - 1])) {
                $run++;
            }
        } else {
            while ($run < $length && $this->less($this->v[$run], $this->v[$run - 1]) === false) {
                $run++;
            }
        }
        if ($run === $length) {
            if ($descending) {
                $this->v = array_reverse($this->v);
            }

            return;
        }
        $limit = 2 * (63 - self::leadingZeros($length | 1)); // 2 * ilog2(len | 1)
        $this->quicksort(0, $length, null, $limit);
    }

    private static function leadingZeros(int $x): int
    {
        $n = 0;
        for ($bit = 63; $bit >= 0 && (($x >> $bit) & 1) === 0; $bit--) {
            $n++;
        }

        return $n;
    }

    /** @param  array{mixed}|null  $ancestor  the ancestor pivot's value, boxed */
    private function quicksort(int $base, int $length, ?array $ancestor, int $limit): void
    {
        while (true) {
            if ($length <= self::SMALL_SORT_THRESHOLD) {
                $this->smallSortGeneral($base, $length);

                return;
            }
            if ($limit === 0) {
                $this->heapsort($base, $length);

                return;
            }
            $limit--;
            $pivot = $this->choosePivot($base, $length);
            if ($ancestor !== null && $this->less($ancestor[0], $this->v[$base + $pivot]) === false) {
                $isLess = $this->isLess;
                $numLt = $this->partition($base, $length, $pivot, fn ($a, $b) => $isLess($b, $a) === false);
                $base += $numLt + 1;
                $length -= $numLt + 1;
                $ancestor = null;

                continue;
            }
            $numLt = $this->partition($base, $length, $pivot, $this->isLess);
            $this->quicksort($base, $numLt, $ancestor, $limit);
            $ancestor = [$this->v[$base + $numLt]];
            $base += $numLt + 1;
            $length -= $numLt + 1;
        }
    }

    private function choosePivot(int $base, int $length): int
    {
        $n8 = intdiv($length, 8);
        $a = $base;
        $b = $base + $n8 * 4;
        $c = $base + $n8 * 7;
        $index = $length < 64 ? $this->median3($a, $b, $c) : $this->median3Rec($a, $b, $c, $n8);

        return $index - $base;
    }

    private function median3Rec(int $a, int $b, int $c, int $n): int
    {
        if ($n * 8 >= 64) {
            $n8 = intdiv($n, 8);
            $a = $this->median3Rec($a, $a + $n8 * 4, $a + $n8 * 7, $n8);
            $b = $this->median3Rec($b, $b + $n8 * 4, $b + $n8 * 7, $n8);
            $c = $this->median3Rec($c, $c + $n8 * 4, $c + $n8 * 7, $n8);
        }

        return $this->median3($a, $b, $c);
    }

    private function median3(int $a, int $b, int $c): int
    {
        $x = $this->less($this->v[$a], $this->v[$b]);
        $y = $this->less($this->v[$a], $this->v[$c]);
        if ($x === $y) {
            $z = $this->less($this->v[$b], $this->v[$c]);

            return $z !== $x ? $c : $b;
        }

        return $a;
    }

    /** partition with partition_lomuto_branchless_cyclic; returns the pivot's final offset. */
    private function partition(int $base, int $length, int $pivot, \Closure $isLess): int
    {
        $v = &$this->v;
        [$v[$base], $v[$base + $pivot]] = [$v[$base + $pivot], $v[$base]];
        $pivotValue = $v[$base];
        $start = $base + 1;
        $count = $length - 1;
        $numLt = 0;
        if ($count > 0) {
            $gapValue = $v[$start];
            $gap = $start;
            for ($k = 1; $k <= $count; $k++) {
                $fromGap = $k === $count;
                $right = $fromGap ? $gapValue : $v[$start + $k];
                $rightIsLt = $isLess($right, $pivotValue);
                $left = $start + $numLt;
                $v[$gap] = $v[$left];
                $v[$left] = $right;
                $gap = $start + $k;
                $numLt += $rightIsLt ? 1 : 0;
            }
        }
        [$v[$base], $v[$base + $numLt]] = [$v[$base + $numLt], $v[$base]];

        return $numLt;
    }

    private function smallSortGeneral(int $base, int $length): void
    {
        if ($length < 2) {
            return;
        }
        $isLess = $this->isLess;
        $half = intdiv($length, 2);
        $scratch = [];
        if ($length >= 16) {
            $this->sort8Stable($base, $scratch, 0, $length, $isLess);
            $this->sort8Stable($base + $half, $scratch, $half, $length + 8, $isLess);
            $presorted = 8;
        } elseif ($length >= 8) {
            self::sort4Stable($this->v, $base, $scratch, 0, $isLess);
            self::sort4Stable($this->v, $base + $half, $scratch, $half, $isLess);
            $presorted = 4;
        } else {
            $scratch[0] = $this->v[$base];
            $scratch[$half] = $this->v[$base + $half];
            $presorted = 1;
        }
        foreach ([0, $half] as $offset) {
            $desired = $offset === 0 ? $half : $length - $half;
            for ($i = $presorted; $i < $desired; $i++) {
                $scratch[$offset + $i] = $this->v[$base + $offset + $i];
                self::insertTail($scratch, $offset, $offset + $i, $isLess);
            }
        }
        self::bidirectionalMerge($scratch, 0, $length, $this->v, $base, $isLess);
    }

    /** sort8_stable: v[$src, +8) into $dst[$dstAt, +8), using $dst[$tmpAt, +8) as scratch. */
    private function sort8Stable(int $src, array &$dst, int $dstAt, int $tmpAt, \Closure $isLess): void
    {
        self::sort4Stable($this->v, $src, $dst, $tmpAt, $isLess);
        self::sort4Stable($this->v, $src + 4, $dst, $tmpAt + 4, $isLess);
        $tmp = [];
        for ($i = 0; $i < 8; $i++) {
            $tmp[$i] = $dst[$tmpAt + $i];
        }
        self::bidirectionalMerge($tmp, 0, 8, $dst, $dstAt, $isLess);
    }

    private static function sort4Stable(array $src, int $at, array &$dst, int $dstAt, \Closure $isLess): void
    {
        $c1 = $isLess($src[$at + 1], $src[$at]);
        $c2 = $isLess($src[$at + 3], $src[$at + 2]);
        $a = $at + ($c1 ? 1 : 0);
        $b = $at + ($c1 ? 0 : 1);
        $c = $at + 2 + ($c2 ? 1 : 0);
        $d = $at + 2 + ($c2 ? 0 : 1);
        $c3 = $isLess($src[$c], $src[$a]);
        $c4 = $isLess($src[$d], $src[$b]);
        $min = $c3 ? $c : $a;
        $max = $c4 ? $b : $d;
        $unknownLeft = $c3 ? $a : ($c4 ? $c : $b);
        $unknownRight = $c4 ? $d : ($c3 ? $b : $c);
        $c5 = $isLess($src[$unknownRight], $src[$unknownLeft]);
        $lo = $c5 ? $unknownRight : $unknownLeft;
        $hi = $c5 ? $unknownLeft : $unknownRight;
        $dst[$dstAt] = $src[$min];
        $dst[$dstAt + 1] = $src[$lo];
        $dst[$dstAt + 2] = $src[$hi];
        $dst[$dstAt + 3] = $src[$max];
    }

    private static function bidirectionalMerge(array $src, int $at, int $length, array &$dst, int $dstAt, \Closure $isLess): void
    {
        $half = intdiv($length, 2);
        $left = $at;
        $right = $at + $half;
        $out = $dstAt;
        $leftRev = $at + $half - 1;
        $rightRev = $at + $length - 1;
        $outRev = $dstAt + $length - 1;
        for ($i = 0; $i < $half; $i++) {
            $takeLeft = $isLess($src[$right], $src[$left]) === false;
            $dst[$out++] = $takeLeft ? $src[$left] : $src[$right];
            $takeLeft ? $left++ : $right++;

            $takeLeft = $isLess($src[$rightRev], $src[$leftRev]) === false;
            $dst[$outRev--] = $takeLeft ? $src[$rightRev] : $src[$leftRev];
            $takeLeft ? $rightRev-- : $leftRev--;
        }
        if ($length % 2 === 1) {
            $dst[$out] = $left < $leftRev + 1 ? $src[$left] : $src[$right];
        }
    }

    private function heapsort(int $base, int $length): void
    {
        $v = &$this->v;
        for ($i = $length + intdiv($length, 2) - 1; $i >= 0; $i--) {
            if ($i >= $length) {
                $node = $i - $length;
            } else {
                [$v[$base], $v[$base + $i]] = [$v[$base + $i], $v[$base]];
                $node = 0;
            }
            $end = min($i, $length);
            while (true) {
                $child = 2 * $node + 1;
                if ($child >= $end) {
                    break;
                }
                if ($child + 1 < $end && $this->less($v[$base + $child], $v[$base + $child + 1])) {
                    $child++;
                }
                if ($this->less($v[$base + $node], $v[$base + $child]) === false) {
                    break;
                }
                [$v[$base + $node], $v[$base + $child]] = [$v[$base + $child], $v[$base + $node]];
                $node = $child;
            }
        }
    }
}
