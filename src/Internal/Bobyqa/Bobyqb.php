<?php

declare(strict_types=1);

namespace Phit\Polyfill\Internal\Bobyqa;

/**
 * The crate's bobyqb.rs: PRIMA's trust-region main loop with its helpers errbd,
 * shiftbase, redrat and redrho, plus the crate's restart and ftol extensions.
 *
 * @internal
 */
final class Bobyqb
{
    /** PRIMA bobyqb.f90 L188: convergence tolerance of the trust-region subproblem solver. */
    private const float TRTOL = 1.0e-2;

    public const float ETA1 = 0.1;

    public const float ETA2 = 0.7;

    public const float GAMMA1 = 0.5;

    public const float GAMMA2 = 2.0;

    /**
     * Run the solver from the revised x0 in $x, overwriting it with the best point.
     *
     * @param  callable(list<float>): (float|int)  $calfun
     * @param  array{float, float}|null  $ftol  [rel, abs], or null when the check is off
     * @param  list<float>  $x
     * @param  list<float>  $xl
     * @param  list<float>  $xu
     * @return array{float, int, int} [f, nf, info]
     */
    public static function solve(
        callable $calfun,
        int $maxfun,
        int $npt,
        float $ftarget,
        ?array $ftol,
        float $rhobeg,
        float $rhoend,
        array &$x,
        array $xl,
        array $xu,
        ?RestartState $restart,
        bool $fullModelDeviations,
    ): array {
        $eta1 = self::ETA1;
        $eta2 = self::ETA2;
        $n = count($x);
        $realmax = Linalg::REALMAX;

        $bmat = [];
        $zmat = [];
        $xpt = [];
        $hq = [];
        $fval = [];
        $pq = [];
        $gopt = array_fill(0, $n, 0.0);
        $sl = array_fill(0, $n, 0.0);
        $su = array_fill(0, $n, 0.0);
        $xbase = array_fill(0, $n, 0.0);
        $d = array_fill(0, $n, 0.0);
        $vlag = [];
        $den = [];
        $distsq = [];
        $ij = [];
        $dnormRec = [$realmax, $realmax];
        $moderrRec = [$realmax, $realmax];

        [$kopt, $nf, $subinfo] = Initialize::initxf(
            $calfun, $maxfun, $ftarget, $rhobeg, $xl, $xu, $x, $ij, $fval, $sl, $su, $xbase, $xpt, $n, $npt,
        );
        $x = Util::xinbd($xbase, Linalg::col($xpt, $n, $kopt), $xl, $xu, $sl, $su);
        $f = $fval[$kopt];

        if ($subinfo === Util::INFO_DFT) {
            Initialize::inith($ij, $xpt, $bmat, $zmat, $n, $npt);
            Initialize::initq($ij, $fval, $xpt, $gopt, $hq, $pq, $n, $npt);
            if (! self::modelIsFinite($gopt, $hq, $pq)) {
                $subinfo = Util::NAN_INF_MODEL;
            }
        }
        if ($subinfo !== Util::INFO_DFT) {
            return [$f, $nf, $subinfo];
        }

        $rho = $rhobeg;
        $delta = $rho;
        $ebound = 0.0;
        $rescued = false;
        $shortd = false;
        $ratio = -1.0;
        $knewTr = null;
        $itest = 0;

        $gamma3 = Linalg::max(1.0, Linalg::min(0.75 * self::GAMMA2, 1.5));

        $maxtr = $restart !== null ? $maxfun * 2 * ($restart->config->maxRestarts + 1) : $maxfun * 2;
        $info = Util::MAXTR_REACHED;
        $dnorm = 0.0;

        $bestF = $realmax;
        $bestX = array_fill(0, $n, 0.0);
        $fireRestart = false;
        $ftolStageF = null;

        for ($tr = 1; $tr <= $maxtr; $tr++) {
            if ($restart !== null) {
                if (! $fireRestart
                    && $restart->overCycleBudget($nf, $maxfun)
                    && $restart->shouldRestart($fval[$kopt], $nf, $maxfun)
                    && $maxfun - $nf > $npt) {
                    $fireRestart = true;
                }
                if ($fireRestart) {
                    $fireRestart = false;
                    $fopt = $fval[$kopt];
                    $rho = $rhobeg;
                    $delta = $rhobeg;
                    $dnormRec = [$realmax, $realmax];
                    $moderrRec = [$realmax, $realmax];
                    $restart->onRestart($fopt, $nf);

                    // hard_rebuild: fold the incumbent into the record, then reseed from it.
                    if ($fval[$kopt] < $bestF) {
                        $bestF = $fval[$kopt];
                        $bestX = Util::xinbd($xbase, Linalg::col($xpt, $n, $kopt), $xl, $xu, $sl, $su);
                    }
                    $x = Util::xinbd($xbase, Linalg::col($xpt, $n, $kopt), $xl, $xu, $sl, $su);
                    Util::reviseX0($x, $xl, $xu, $rhobeg);
                    [$koptNew, $nfInit, $subinfo] = Initialize::initxf(
                        $calfun, $maxfun - $nf, $ftarget, $rhobeg, $xl, $xu, $x, $ij, $fval, $sl, $su, $xbase, $xpt, $n, $npt,
                    );
                    $kopt = $koptNew;
                    $nf += $nfInit;
                    $x = Util::xinbd($xbase, Linalg::col($xpt, $n, $kopt), $xl, $xu, $sl, $su);
                    $f = $fval[$kopt];
                    if ($subinfo === Util::INFO_DFT) {
                        Initialize::inith($ij, $xpt, $bmat, $zmat, $n, $npt);
                        Initialize::initq($ij, $fval, $xpt, $gopt, $hq, $pq, $n, $npt);
                        if (! self::modelIsFinite($gopt, $hq, $pq)) {
                            $subinfo = Util::NAN_INF_MODEL;
                        }
                    }

                    $ebound = 0.0;
                    $rescued = false;
                    $shortd = false;
                    $ratio = -1.0;
                    $knewTr = null;
                    $itest = 0;
                    $ftolStageF = null;
                    if ($subinfo !== Util::INFO_DFT) {
                        $info = $subinfo;
                        break;
                    }

                    continue;
                }
            }

            $crvmin = TrustRegion::trsbox(
                $delta, $gopt, $hq, $pq, $sl, $su, self::TRTOL, Linalg::col($xpt, $n, $kopt), $xpt, $n, $npt, $d,
            );
            $dnorm = Linalg::min($delta, Linalg::norm($d));
            $shortd = $dnorm <= 0.5 * $rho;

            $qred = -PowAlg::quadinc($d, $xpt, $n, $npt, $gopt, $pq, $hq);
            $qredThreshold = Util::F32_1E_6 * ($rho * $rho);
            $trfail = is_nan($qred) || is_nan($qredThreshold) || $qred <= $qredThreshold;

            if ($shortd || $trfail) {
                $delta = 0.1 * $delta;
                if ($delta <= $gamma3 * $rho) {
                    $delta = $rho;
                }
                $ebound = self::errbd(
                    $crvmin, $d, $gopt, $hq, $moderrRec, $pq, $rho, $sl, $su, Linalg::col($xpt, $n, $kopt), $xpt, $n, $npt,
                );
            } else {
                $xnew = [];
                for ($i = 0; $i < $n; $i++) {
                    $xnew[] = $xpt[$i + $kopt * $n] + $d[$i];
                }
                $x = Util::xinbd($xbase, $xnew, $xl, $xu, $sl, $su);
                $f = Util::evaluate($calfun, $x);
                $nf++;
                $rescued = false;

                $subinfo = Util::checkexit($maxfun, $nf, $f, $ftarget, $x);
                if ($subinfo !== Util::INFO_DFT) {
                    $info = $subinfo;
                    break;
                }

                $dnormRec = [$dnormRec[1], $dnorm];
                $moderr = $f - $fval[$kopt] + $qred;
                $moderrRec = [$moderrRec[1], $moderr];

                $ratio = self::redrat($fval[$kopt] - $f, $qred, $eta1);

                $delta = TrustRegion::trrad($delta, $dnorm, $eta1, $eta2, self::GAMMA1, self::GAMMA2, $ratio);
                if ($delta <= $gamma3 * $rho) {
                    $delta = $rho;
                }

                $ximproved = $f < $fval[$kopt];

                $refBeta = PowAlg::calvlagAndDen($kopt, $bmat, $d, $xpt, $zmat, $n, $npt, $vlag, $den);
                $modelUnchanged = true;
                $vlagAbsSum = Linalg::sumAbs($vlag);
                $vmax = 0.0;
                for ($k = 0; $k < $npt; $k++) {
                    $sq = $vlag[$k] * $vlag[$k];
                    if ($k === 0 || $sq > $vmax) {
                        $vmax = $sq;
                    }
                }
                $denThreshold = $fullModelDeviations ? 0.5 * $vmax : $vmax;
                $denOk = false;
                foreach ($den as $value) {
                    if ($value > $denThreshold) {
                        $denOk = true;
                        break;
                    }
                }
                $toRescue = $ximproved && (! is_finite($vlagAbsSum) || ! $denOk);
                if ($toRescue) {
                    if ($rescued) {
                        $info = Util::DAMAGING_ROUNDING;
                        break;
                    }
                    $subinfo = Rescue::rescue(
                        $calfun, $maxfun, $delta, $ftarget, $xl, $xu, $kopt, $nf, $fval, $gopt, $hq, $pq,
                        $sl, $su, $xbase, $xpt, $bmat, $zmat, $n, $npt,
                    );
                    if ($subinfo !== Util::INFO_DFT) {
                        $info = $subinfo;
                        break;
                    }
                    $rescued = true;
                    $modelUnchanged = false;
                    $dnormRec = [$realmax, $realmax];
                    $moderrRec = [$realmax, $realmax];

                    for ($i = 0; $i < $n; $i++) {
                        $d[$i] = Linalg::max($sl[$i], Linalg::min($su[$i], $d[$i])) - $xpt[$i + $kopt * $n];
                    }
                    $moderr = $f - $fval[$kopt] - PowAlg::quadinc($d, $xpt, $n, $npt, $gopt, $pq, $hq);
                    $ximproved = $f < $fval[$kopt];
                }

                $knewTr = Geometry::setdropTr(
                    $kopt, $ximproved, $bmat, $d, $rho, $xpt, $zmat, $n, $npt, $modelUnchanged ? $den : null,
                );

                if ($knewTr !== null) {
                    $xdrop = Linalg::col($xpt, $n, $knewTr);
                    $xosav = Linalg::col($xpt, $n, $kopt);
                    Update::updateh($knewTr, $kopt, $d, $xpt, $bmat, $zmat, $n, $npt, $modelUnchanged ? [$vlag, $refBeta] : null);
                    $xnewClamped = [];
                    for ($i = 0; $i < $n; $i++) {
                        $xnewClamped[] = Linalg::max($sl[$i], Linalg::min($su[$i], $xosav[$i] + $d[$i]));
                    }
                    Update::updatexf($knewTr, $ximproved, $f, $xnewClamped, $kopt, $fval, $xpt, $n);
                    Update::updateq($knewTr, $ximproved, $bmat, $d, $moderr, $xdrop, $xosav, $xpt, $zmat, $gopt, $hq, $pq, $n, $npt);
                    $fvalShift = [];
                    for ($k = 0; $k < $npt; $k++) {
                        $fvalShift[] = $fval[$k] - $fval[$kopt];
                    }
                    Update::tryqalt(
                        $bmat, $fvalShift, $ratio, $sl, $su, Linalg::col($xpt, $n, $kopt), $xpt, $zmat, $itest, $gopt, $hq, $pq, $n, $npt,
                    );
                    if (! self::modelIsFinite($gopt, $hq, $pq)) {
                        $info = Util::NAN_INF_MODEL;
                        break;
                    }
                }
            }

            $accurateMod = abs($moderrRec[0]) <= $ebound && abs($moderrRec[1]) <= $ebound
                && $dnormRec[0] <= $rho && $dnormRec[1] <= $rho;
            $koff = $kopt * $n;
            $distsq = [];
            for ($k = 0; $k < $npt; $k++) {
                $offset = $k * $n;
                $sq = 0.0;
                for ($i = 0; $i < $n; $i++) {
                    $diff = $xpt[$offset + $i] - $xpt[$koff + $i];
                    $sq += $diff * $diff;
                }
                $distsq[] = $sq;
            }
            $closeBound = Linalg::max($delta * $delta, (10.0 * $rho) * (10.0 * $rho));
            $closeItpset = true;
            foreach ($distsq as $value) {
                if (is_nan($value) || is_nan($closeBound) || $value > $closeBound) {
                    $closeItpset = false;
                    break;
                }
            }
            $adequateGeo = ($shortd && $accurateMod) || $closeItpset;
            $smallTrrad = Linalg::max($delta, $dnorm) <= $rho;

            $badTrstep = $shortd || $trfail || $ratio <= $eta1 || $knewTr === null;
            $improveGeo = $badTrstep && ! $adequateGeo;
            $badTrstep = $shortd || $trfail || $ratio <= 0.0 || $knewTr === null;
            $reduceRho = $badTrstep && $adequateGeo && $smallTrrad;

            if ($improveGeo) {
                $knewGeo = 0;
                for ($k = 0; $k < $npt; $k++) {
                    if ($k === 0 || $distsq[$k] > $distsq[$knewGeo]) {
                        $knewGeo = $k;
                    }
                }
                $maxDistsq = 0.0;
                for ($k = 0; $k < $npt; $k++) {
                    if ($k === 0 || $distsq[$k] > $maxDistsq) {
                        $maxDistsq = $distsq[$k];
                    }
                }
                $delbar = Linalg::max(Linalg::min(0.1 * sqrt($maxDistsq), $delta), $rho);

                $d = Geometry::geostep($knewGeo, $kopt, $bmat, $delbar, $sl, $su, $xpt, $zmat, $n, $npt);

                $geoBeta = PowAlg::calvlagAndDen($kopt, $bmat, $d, $xpt, $zmat, $n, $npt, $vlag, $den);
                $vlagAbsSum = Linalg::sumAbs($vlag);
                $geoThreshold = 0.5 * ($vlag[$knewGeo] * $vlag[$knewGeo]);
                $toRescue = ! is_finite($vlagAbsSum) || is_nan($den[$knewGeo]) || is_nan($geoThreshold)
                    || $den[$knewGeo] <= $geoThreshold;
                if ($toRescue) {
                    if ($rescued) {
                        $info = Util::DAMAGING_ROUNDING;
                        break;
                    }
                    $subinfo = Rescue::rescue(
                        $calfun, $maxfun, $delta, $ftarget, $xl, $xu, $kopt, $nf, $fval, $gopt, $hq, $pq,
                        $sl, $su, $xbase, $xpt, $bmat, $zmat, $n, $npt,
                    );
                    if ($subinfo !== Util::INFO_DFT) {
                        $info = $subinfo;
                        break;
                    }
                    $rescued = true;
                    $dnormRec = [$realmax, $realmax];
                    $moderrRec = [$realmax, $realmax];
                } else {
                    $xnew = [];
                    for ($i = 0; $i < $n; $i++) {
                        $xnew[] = $xpt[$i + $kopt * $n] + $d[$i];
                    }
                    $x = Util::xinbd($xbase, $xnew, $xl, $xu, $sl, $su);
                    $f = Util::evaluate($calfun, $x);
                    $nf++;
                    $rescued = false;

                    $subinfo = Util::checkexit($maxfun, $nf, $f, $ftarget, $x);
                    if ($subinfo !== Util::INFO_DFT) {
                        $info = $subinfo;
                        break;
                    }

                    $dnorm = Linalg::min($delbar, Linalg::norm($d));
                    $dnormRec = [$dnormRec[1], $dnorm];
                    $moderr = $f - $fval[$kopt] - PowAlg::quadinc($d, $xpt, $n, $npt, $gopt, $pq, $hq);
                    $moderrRec = [$moderrRec[1], $moderr];

                    $ximproved = $f < $fval[$kopt];

                    $xdrop = Linalg::col($xpt, $n, $knewGeo);
                    $xosav = Linalg::col($xpt, $n, $kopt);
                    Update::updateh($knewGeo, $kopt, $d, $xpt, $bmat, $zmat, $n, $npt, [$vlag, $geoBeta]);
                    $xnewClamped = [];
                    for ($i = 0; $i < $n; $i++) {
                        $xnewClamped[] = Linalg::max($sl[$i], Linalg::min($su[$i], $xosav[$i] + $d[$i]));
                    }
                    Update::updatexf($knewGeo, $ximproved, $f, $xnewClamped, $kopt, $fval, $xpt, $n);
                    Update::updateq($knewGeo, $ximproved, $bmat, $d, $moderr, $xdrop, $xosav, $xpt, $zmat, $gopt, $hq, $pq, $n, $npt);
                    if (! self::modelIsFinite($gopt, $hq, $pq)) {
                        $info = Util::NAN_INF_MODEL;
                        break;
                    }
                }
            }

            if ($reduceRho) {
                if ($ftol !== null && $rho > $rhoend) {
                    $fopt = $fval[$kopt];
                    if ($ftolStageF !== null) {
                        $df = $ftolStageF - $fopt;
                        if ($df <= $ftol[0] * Linalg::max(abs($fopt), 1.0) + $ftol[1]) {
                            $info = Util::FTOL_REACHED;
                            break;
                        }
                    }
                    $ftolStageF = $fopt;
                }
                if ($restart !== null) {
                    $fopt = $fval[$kopt];
                    if ($rho <= $rhoend) {
                        $fire = ($restart->config->cycleBudgetFrac <= 0.0 || $restart->overCycleBudget($nf, $maxfun))
                            && $restart->shouldRestart($fopt, $nf, $maxfun);
                    } else {
                        $fire = $restart->noteReduction($fopt) && $restart->shouldRestart($fopt, $nf, $maxfun);
                    }
                    if ($fire && $maxfun - $nf > $npt) {
                        $fireRestart = true;

                        continue;
                    }
                }
                if ($rho <= $rhoend) {
                    $info = Util::SMALL_TR_RADIUS;
                    break;
                }
                $delta = Linalg::max(0.5 * $rho, self::redrho($rho, $rhoend));
                $rho = self::redrho($rho, $rhoend);
                $dnormRec = [$realmax, $realmax];
                $moderrRec = [$realmax, $realmax];
            }

            $xoptSq = 0.0;
            for ($i = 0; $i < $n; $i++) {
                $xoptSq += $xpt[$i + $kopt * $n] * $xpt[$i + $kopt * $n];
            }
            if ($xoptSq >= 1.0e3 * ($delta * $delta)) {
                for ($i = 0; $i < $n; $i++) {
                    $xi = $xpt[$i + $kopt * $n];
                    $sl[$i] = Linalg::min($sl[$i] - $xi, 0.0);
                    $su[$i] = Linalg::max($su[$i] - $xi, 0.0);
                }
                self::shiftbase($kopt, $xbase, $xpt, $zmat, $bmat, $pq, $hq, $n, $npt);
                for ($i = 0; $i < $n; $i++) {
                    $xbase[$i] = Linalg::max($xl[$i], Linalg::min($xu[$i], $xbase[$i]));
                }
            }
        }

        if ($info === Util::SMALL_TR_RADIUS && $shortd && $dnorm > 0.1 * $rhoend && $nf < $maxfun) {
            $xnew = [];
            for ($i = 0; $i < $n; $i++) {
                $xnew[] = $xpt[$i + $kopt * $n] + $d[$i];
            }
            $x = Util::xinbd($xbase, $xnew, $xl, $xu, $sl, $su);
            $f = Util::evaluate($calfun, $x);
            $nf++;
        }

        if ($fval[$kopt] < $f || is_nan($f)) {
            $x = Util::xinbd($xbase, Linalg::col($xpt, $n, $kopt), $xl, $xu, $sl, $su);
            $f = $fval[$kopt];
        }

        if ($bestF < $f || (is_nan($f) && $bestF < $realmax)) {
            $x = $bestX;
            $f = $bestF;
        }

        if ($restart !== null && $restart->cycleBoundaries !== []
            && in_array($info, [Util::DAMAGING_ROUNDING, Util::NAN_INF_MODEL, Util::NAN_INF_X, Util::NAN_INF_F], true)) {
            $info = Util::SMALL_TR_RADIUS;
        }

        return [$f, $nf, $info];
    }

    /** PRIMA ratio.f90 redrat: the reduction ratio of a trust-region step, Inf/NaN-safe. */
    private static function redrat(float $ared, float $pred, float $rshrink): float
    {
        if (is_nan($ared)) {
            return -Linalg::REALMAX;
        }
        if (is_nan($pred) || $pred <= 0.0) {
            return $ared > 0.0 ? 0.5 * $rshrink : -Linalg::REALMAX;
        }
        if (is_infinite($pred) && $pred > 0.0 && is_infinite($ared)) {
            return $ared > 0.0 ? 1.0 : -Linalg::REALMAX;
        }

        return $ared / $pred;
    }

    /** PRIMA redrho.f90: the reduced RHO. */
    private static function redrho(float $rho, float $rhoend): float
    {
        $rhoRatio = $rho / $rhoend;
        if ($rhoRatio > 250.0) {
            return 0.1 * $rho;
        }
        if ($rhoRatio <= 16.0) {
            return $rhoend;
        }

        return sqrt($rhoRatio) * $rhoend;
    }

    /**
     * @param  list<float>  $gopt
     * @param  list<float>  $hq
     * @param  list<float>  $pq
     */
    private static function modelIsFinite(array $gopt, array $hq, array $pq): bool
    {
        return Linalg::allFinite($gopt) && Linalg::allFinite($hq) && Linalg::allFinite($pq);
    }

    /**
     * PRIMA shiftbase_lfqint: move XBASE to XBASE + XOPT, updating BMAT and HQ.
     *
     * @param  list<float>  $xbase
     * @param  list<float>  $xpt
     * @param  list<float>  $zmat
     * @param  list<float>  $bmat
     * @param  list<float>  $pq
     * @param  list<float>  $hq
     */
    private static function shiftbase(int $kopt, array &$xbase, array &$xpt, array $zmat, array &$bmat, array $pq, array &$hq, int $n, int $npt): void
    {
        $nz = $npt - $n - 1;
        $xopt = Linalg::col($xpt, $n, $kopt);
        $xoptsq = Linalg::inprod($xopt, $xopt);

        $sxpt = Linalg::matprod12($xopt, $xpt, $n, $npt);
        for ($k = 0; $k < $npt; $k++) {
            $sxpt[$k] -= 0.5 * $xoptsq;
        }

        $qxoptq = 0.25 * $xoptsq;
        $ymat = [];
        for ($k = 0; $k < $npt; $k++) {
            for ($i = 0; $i < $n; $i++) {
                $xptxav = $xpt[$i + $k * $n] - 0.5 * $xopt[$i];
                $ymat[] = $sxpt[$k] * $xptxav + $qxoptq * $xopt[$i];
            }
        }

        $bmatNpt = array_slice($bmat, 0, $n * $npt);
        $ymatT = [];
        for ($i = 0; $i < $n; $i++) {
            for ($k = 0; $k < $npt; $k++) {
                $ymatT[] = $ymat[$i + $k * $n];
            }
        }
        $bymat = Linalg::matprod22($bmatNpt, $n, $npt, $ymatT, $n);
        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $n; $j++) {
                $bmat[$i + ($npt + $j) * $n] += $bymat[$i + $j * $n] + $bymat[$j + $i * $n];
            }
        }

        $yzmat = Linalg::matprod22($ymat, $n, $npt, $zmat, $nz);
        $yzmatT = [];
        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $nz; $j++) {
                $yzmatT[] = $yzmat[$i + $j * $n];
            }
        }
        $yzyz = Linalg::matprod22($yzmat, $n, $nz, $yzmatT, $n);
        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $n; $j++) {
                $bmat[$i + ($npt + $j) * $n] += $yzyz[$i + $j * $n];
            }
        }

        $zmatT = [];
        for ($k = 0; $k < $npt; $k++) {
            for ($j = 0; $j < $nz; $j++) {
                $zmatT[] = $zmat[$k + $j * $npt];
            }
        }
        $yzZt = Linalg::matprod22($yzmat, $n, $nz, $zmatT, $npt);
        $limit = $n * $npt;
        for ($index = 0; $index < $limit; $index++) {
            $bmat[$index] += $yzZt[$index];
        }

        $v = Linalg::matprod21($xpt, $n, $npt, $pq);
        $pqSum = Linalg::sum($pq);
        for ($i = 0; $i < $n; $i++) {
            $v[$i] -= (0.5 * $pqSum) * $xopt[$i];
        }

        for ($i = 0; $i < $n; $i++) {
            for ($j = 0; $j < $n; $j++) {
                $hq[$i + $j * $n] = ($v[$i] * $xopt[$j] + $v[$j] * $xopt[$i]) + $hq[$i + $j * $n];
            }
        }

        for ($i = 0; $i < $n; $i++) {
            $xbase[$i] += $xopt[$i];
        }
        for ($k = 0; $k < $npt; $k++) {
            for ($i = 0; $i < $n; $i++) {
                $xpt[$i + $k * $n] -= $xopt[$i];
            }
        }
        for ($i = 0; $i < $n; $i++) {
            $xpt[$i + $kopt * $n] = 0.0;
        }
    }

    /**
     * PRIMA bobyqb.f90 errbd: the bound used to test whether recent model errors are small.
     *
     * @param  list<float>  $d
     * @param  list<float>  $gopt
     * @param  list<float>  $hq
     * @param  array{float, float}  $moderrRec
     * @param  list<float>  $pq
     * @param  list<float>  $sl
     * @param  list<float>  $su
     * @param  list<float>  $xopt
     * @param  list<float>  $xpt
     */
    private static function errbd(
        float $crvmin,
        array $d,
        array $gopt,
        array $hq,
        array $moderrRec,
        array $pq,
        float $rho,
        array $sl,
        array $su,
        array $xopt,
        array $xpt,
        int $n,
        int $npt,
    ): float {
        $hm = PowAlg::hessMul($d, $xpt, $n, $npt, $pq, $hq);
        $maxAbsModerr = 0.0;
        foreach ($moderrRec as $value) {
            $maxAbsModerr = Linalg::max($maxAbsModerr, abs($value));
        }

        $bfirst = [];
        for ($i = 0; $i < $n; $i++) {
            $xnew = $xopt[$i] + $d[$i];
            $gnew = $gopt[$i] + $hm[$i];
            $b = $maxAbsModerr;
            if ($xnew <= $sl[$i]) {
                $b = $gnew * $rho;
            }
            if ($xnew >= $su[$i]) {
                $b = -$gnew * $rho;
            }
            $bfirst[] = $b;
        }

        $rhosq = $rho * $rho;
        $v = array_fill(0, $n, 0.0);
        for ($k = 0; $k < $npt; $k++) {
            $pqk = $pq[$k];
            $offset = $k * $n;
            for ($i = 0; $i < $n; $i++) {
                $xk = $xpt[$offset + $i];
                $v[$i] += ($xk * $xk) * $pqk;
            }
        }

        $ebound = INF;
        for ($i = 0; $i < $n; $i++) {
            $bsecond = Linalg::midpoint($hq[$i + $i * $n], $v[$i]) * $rhosq;
            $value = Linalg::max($bfirst[$i], $bfirst[$i] + $bsecond);
            if ($value < $ebound) {
                $ebound = $value;
            }
        }
        if ($crvmin > 0.0) {
            $ebound = Linalg::min($ebound, 0.125 * $crvmin * $rhosq);
        }

        return $ebound;
    }
}
