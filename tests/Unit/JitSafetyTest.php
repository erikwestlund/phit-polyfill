<?php

declare(strict_types=1);

// PHP 8.5's function JIT (opcache.jit=1235) can fail to restore a value when
// compiled code is entered partway through a function (php-src #24214); see
// docs/notes/function-jit.md. Two call shapes in src/ are fixed by rule:
//
// - Internal functions are called fully qualified (\abs($x)). In a namespace an
//   unqualified call compiles to a by-name call, and the JIT adds an entry
//   point after each one.
// - count() is the exception and stays unqualified. \count() compiles to the
//   COUNT opcode, whose inlined result the JIT fails to restore at a loop's
//   on-stack-replacement entry (Dense::multiply looped forever).

/** @return list<string> "File.php:line: name()" for each call matching $offending */
function jitSafetyOffenders(callable $offending): array
{
    $offenders = [];
    $files = new RecursiveIteratorIterator(new RecursiveDirectoryIterator(dirname(__DIR__, 2).'/src'));
    foreach ($files as $file) {
        if ($file->getExtension() !== 'php') {
            continue;
        }
        $tokens = array_values(array_filter(
            PhpToken::tokenize(file_get_contents($file->getPathname())),
            fn (PhpToken $t) => ! $t->isIgnorable(),
        ));
        foreach ($tokens as $i => $token) {
            if (! $token->is([T_STRING, T_NAME_FULLY_QUALIFIED]) || ($tokens[$i + 1]->text ?? '') !== '(') {
                continue;
            }
            $previous = $tokens[$i - 1] ?? null;
            if ($previous !== null && $previous->is([T_FUNCTION, T_OBJECT_OPERATOR, T_NULLSAFE_OBJECT_OPERATOR, T_DOUBLE_COLON, T_NEW, T_CONST])) {
                continue;
            }
            $qualified = $token->is(T_NAME_FULLY_QUALIFIED);
            $name = strtolower(ltrim($token->text, '\\'));
            if (function_exists($name) && (new ReflectionFunction($name))->isInternal() && $offending($name, $qualified)) {
                $offenders[] = $file->getFilename().':'.$token->line.': '.$token->text.'()';
            }
        }
    }

    return $offenders;
}

it('fully qualifies calls to internal functions other than count()', function () {
    expect(jitSafetyOffenders(fn (string $name, bool $qualified) => ! $qualified && $name !== 'count'))->toBe([]);
});

it('leaves count() unqualified', function () {
    expect(jitSafetyOffenders(fn (string $name, bool $qualified) => $qualified && $name === 'count'))->toBe([]);
});
