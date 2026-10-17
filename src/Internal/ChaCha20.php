<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal;

/**
 * The extension's sampling generator: rand_chacha 0.9's ChaCha20Rng seeded by
 * rand_core 0.9's seed_from_u64, read as u64 words (BlockRng::next_u64) and
 * mapped to rand 0.9's Open01 doubles. The keystream is standard ChaCha20 with
 * a 64-bit block counter and a zero stream identifier.
 *
 * @internal
 */
final class ChaCha20
{
    /** Keystream bytes produced per refill: rand_chacha's four-block buffer. */
    private const int BUFFER_BYTES = 256;

    /** @var list<int> eight little-endian key words */
    private array $key;

    private int $counter = 0;

    /** @var list<int> buffered 32-bit words */
    private array $words = [];

    private int $index = 0;

    private bool $openssl;

    public function __construct(int $seed)
    {
        $this->key = self::seedWords($seed);
        $this->openssl = extension_loaded('openssl') && in_array('chacha20', openssl_get_cipher_methods(), true);
    }

    /**
     * rand_core's seed_from_u64: eight PCG32 outputs. The state is a u64 held
     * in a PHP int, so the LCG step multiplies in 32-bit halves.
     *
     * @return list<int>
     */
    private static function seedWords(int $state): array
    {
        $mulHi = 0x5851F42D;
        $mulLo = 0x4C957F2D;
        $incHi = 0xA17654E4;
        $incLo = 0x6FBE17F3;
        $hi = ($state >> 32) & 0xFFFFFFFF;
        $lo = $state & 0xFFFFFFFF;
        $words = [];
        for ($i = 0; $i < 8; $i++) {
            // state = state * MUL + INC (mod 2^64)
            [$productHi, $newLo] = self::multiply($lo, $mulLo);
            $newHi = ($productHi + self::multiply($hi, $mulLo)[1] + self::multiply($lo, $mulHi)[1]) & 0xFFFFFFFF;
            $sumLo = $newLo + $incLo;
            $lo = $sumLo & 0xFFFFFFFF;
            $hi = ($newHi + $incHi + ($sumLo >> 32)) & 0xFFFFFFFF;
            // xorshifted = (((state >> 18) ^ state) >> 27) as u32
            $shiftedHi = $hi >> 18;
            $shiftedLo = (($lo >> 18) | ($hi << 14)) & 0xFFFFFFFF;
            $xHi = $shiftedHi ^ $hi;
            $xLo = $shiftedLo ^ $lo;
            $xorshifted = (($xLo >> 27) | ($xHi << 5)) & 0xFFFFFFFF;
            $rot = $hi >> 27;
            $words[] = $rot === 0 ? $xorshifted : (($xorshifted >> $rot) | ($xorshifted << (32 - $rot))) & 0xFFFFFFFF;
        }

        return $words;
    }

    /**
     * The 64-bit product of two u32 values as [high word, low word], split so
     * that no intermediate leaves PHP's signed 64-bit range.
     *
     * @return array{0: int, 1: int}
     */
    private static function multiply(int $a, int $b): array
    {
        $low = ($a & 0xFFFF) * $b;
        $high = ($a >> 16) * $b;
        $t = $low + (($high & 0xFFFF) << 16);

        return [(($high >> 16) + ($t >> 32)) & 0xFFFFFFFF, $t & 0xFFFFFFFF];
    }

    /** The next Open01 double: (2f + 1) / 2^53 from the top 52 bits of a u64. */
    public function open01(): float
    {
        if ($this->index >= count($this->words)) {
            $this->refill();
        }
        $lo = $this->words[$this->index];
        $hi = $this->words[$this->index + 1];
        $this->index += 2;
        $fraction = ($hi << 20) | ($lo >> 12);

        return (float) (2 * $fraction + 1) * 1.1102230246251565E-16; // 2^-53
    }

    private function refill(): void
    {
        $bytes = $this->openssl ? $this->opensslBlocks() : $this->phpBlocks();
        $this->words = array_values(unpack('V*', $bytes));
        $this->index = 0;
        $this->counter += 4;
    }

    /** OpenSSL's chacha20 takes the 128-bit IV as state words 12 to 15. */
    private function opensslBlocks(): string
    {
        $key = pack('V*', ...$this->key);
        $iv = pack('VVVV', $this->counter & 0xFFFFFFFF, ($this->counter >> 32) & 0xFFFFFFFF, 0, 0);
        $out = openssl_encrypt(str_repeat("\0", self::BUFFER_BYTES), 'chacha20', $key, OPENSSL_RAW_DATA, $iv);
        if ($out === false) {
            return $this->phpBlocks();
        }

        return $out;
    }

    private function phpBlocks(): string
    {
        $out = '';
        for ($b = 0; $b < 4; $b++) {
            $counter = $this->counter + $b;
            $state = [
                0x61707865, 0x3320646E, 0x79622D32, 0x6B206574,
                ...$this->key,
                $counter & 0xFFFFFFFF, ($counter >> 32) & 0xFFFFFFFF, 0, 0,
            ];
            $x = $state;
            for ($round = 0; $round < 10; $round++) {
                self::quarter($x, 0, 4, 8, 12);
                self::quarter($x, 1, 5, 9, 13);
                self::quarter($x, 2, 6, 10, 14);
                self::quarter($x, 3, 7, 11, 15);
                self::quarter($x, 0, 5, 10, 15);
                self::quarter($x, 1, 6, 11, 12);
                self::quarter($x, 2, 7, 8, 13);
                self::quarter($x, 3, 4, 9, 14);
            }
            for ($i = 0; $i < 16; $i++) {
                $x[$i] = ($x[$i] + $state[$i]) & 0xFFFFFFFF;
            }
            $out .= pack('V*', ...$x);
        }

        return $out;
    }

    /** @param list<int> $x */
    private static function quarter(array &$x, int $a, int $b, int $c, int $d): void
    {
        $x[$a] = ($x[$a] + $x[$b]) & 0xFFFFFFFF;
        $v = $x[$d] ^ $x[$a];
        $x[$d] = (($v << 16) | ($v >> 16)) & 0xFFFFFFFF;
        $x[$c] = ($x[$c] + $x[$d]) & 0xFFFFFFFF;
        $v = $x[$b] ^ $x[$c];
        $x[$b] = (($v << 12) | ($v >> 20)) & 0xFFFFFFFF;
        $x[$a] = ($x[$a] + $x[$b]) & 0xFFFFFFFF;
        $v = $x[$d] ^ $x[$a];
        $x[$d] = (($v << 8) | ($v >> 24)) & 0xFFFFFFFF;
        $x[$c] = ($x[$c] + $x[$d]) & 0xFFFFFFFF;
        $v = $x[$b] ^ $x[$c];
        $x[$b] = (($v << 7) | ($v >> 25)) & 0xFFFFFFFF;
    }

    /** Force the pure-PHP keystream (for testing it against OpenSSL). */
    public static function withoutOpenssl(int $seed): self
    {
        $rng = new self($seed);
        $rng->openssl = false;

        return $rng;
    }
}
