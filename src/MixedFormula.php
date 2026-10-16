<?php

declare(strict_types=1);

namespace PhpStats\Polyfill;

use PhpStats\Polyfill\Formula\Design;
use PhpStats\Polyfill\Formula\Frontend;
use PhpStats\Polyfill\Internal\LinearMixed;
use PhpStats\Polyfill\Internal\LogisticRandomIntercept;
use PhpStats\Polyfill\Internal\ModelData;

/**
 * lme4-style formula interfaces over the mixed-model cores: stats_lmer and
 * stats_glmer (the extension's formula.rs).
 */
final class MixedFormula
{
    /** stats_lmer */
    public static function lmer(
        string $formula,
        mixed $data,
        bool $reml = true,
        string $missing = 'drop',
        mixed $weights = null,
        ?array $start = null,
        int $max_evaluations = 20000,
        float $tolerance = 1e-7,
        string $optimizer = 'bobyqa',
        float $objective_tolerance = 1e-12,
    ): array {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($missing, 'missing');
        Frontend::stringArgument($optimizer, 'optimizer');
        if ($optimizer !== 'bobyqa' && $optimizer !== 'nelder_mead') {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, "optimizer must be 'bobyqa' or 'nelder_mead'");
        }
        $startValues = null;
        if ($start !== null) {
            $startValues = [];
            foreach ($start as $value) {
                $value = ModelData::number($value);
                if ($value === null) {
                    throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'start cannot contain missing values');
                }
                $startValues[] = $value;
            }
        }
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        self::requireRandom($design, 'stats_lmer');
        $y = Frontend::numericResponse($design, 'stats_lmer');
        $terms = [];
        $termLabels = [];
        foreach ($design->random as $term) {
            [$indices, $labels] = ModelData::groupIndices($term['groups']);
            $termLabels[] = $labels;
            $terms[] = ['z' => $term['z'], 'groups' => $indices, 'correlated' => $term['correlated']];
        }
        if ($max_evaluations < 0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'max_evaluations must be positive');
        }
        $result = LinearMixed::fit($design->x, $y, $terms, [
            'reml' => $reml,
            'optimizer' => $optimizer,
            'objective_tolerance' => $objective_tolerance,
            'offset' => $design->offset ?? array_fill(0, count($y), 0.0),
            'weights' => $w,
            'start' => $startValues,
            'max_evaluations' => $max_evaluations,
            'tolerance' => $tolerance,
        ]);
        $output = MixedModels::lmmOutput($result, $reml, $termLabels);
        foreach ($design->random as $i => $term) {
            $output['random'][$i]['label'] = $term['label'];
            $output['random'][$i]['group'] = $term['group'];
            $output['random'][$i]['columns'] = $term['columns'];
            $output['random'][$i]['correlated'] = $term['correlated'];
        }
        $output = Frontend::addNames($output, $design);

        return Frontend::finish($output, $design, 'stats_lmer', $optimizer);
    }

    /** stats_glmer */
    public static function glmer(
        string $formula,
        mixed $data,
        string $family = 'binomial',
        string $link = '',
        string $missing = 'drop',
        mixed $weights = null,
        int $quadrature_points = 1,
        int $max_evaluations = 20000,
        float $tolerance = 1e-7,
    ): array {
        Frontend::stringArgument($formula, 'formula');
        Frontend::stringArgument($family, 'family');
        Frontend::stringArgument($link, 'link');
        Frontend::stringArgument($missing, 'missing');
        if ($family !== 'binomial' || ($link !== '' && $link !== 'logit')) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, "stats_glmer supports family 'binomial' with the logit link only (got '$family'/'$link')");
        }
        [$design, $w] = Frontend::prepare($formula, $data, $missing, $weights, true);
        self::requireRandom($design, 'stats_glmer');
        if (count($design->random) !== 1 || $design->random[0]['columns'] !== ['(Intercept)']) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'stats_glmer supports exactly one random intercept, (1 | g); random slopes and several terms are not implemented');
        }
        $y = Glm::binomialResponse($design, $w);
        [$indices, $labels] = ModelData::groupIndices($design->random[0]['groups']);
        if ($quadrature_points < 0) {
            throw ErrorCode::exception(ErrorCode::UNSUPPORTED, 'quadrature_points must be positive');
        }
        if ($max_evaluations < 0) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, 'max_evaluations must be positive');
        }
        $r = LogisticRandomIntercept::fit($design->x, $y, $indices, [
            'offset' => $design->offset ?? array_fill(0, count($y), 0.0),
            'weights' => $w,
            'quadrature_points' => $quadrature_points,
            'max_evaluations' => $max_evaluations,
            'tolerance' => $tolerance,
        ]);
        $output = MixedModels::logisticMixedOutput($r, $labels);
        $output['random'] = [[
            'label' => $design->random[0]['label'],
            'group' => $design->random[0]['group'],
            'columns' => $design->random[0]['columns'],
        ]];
        $output = Frontend::addNames($output, $design);

        return Frontend::finish($output, $design, 'stats_glmer', 'bobyqa');
    }

    private static function requireRandom(Design $design, string $function): void
    {
        if ($design->random === []) {
            throw ErrorCode::exception(ErrorCode::INVALID_INPUT, "$function needs at least one random-effect term such as (1 | g); use stats_lm or stats_glm for fixed effects only");
        }
    }
}
