//! Reference harness for the PHP port of the `bobyqa` crate (=0.3.2).
//!
//! `cargo run --release --offline -- fixture` prints the fixed problems used by
//! tests/Unit/BobyqaTest.php; `-- sweep <count> <seed>` prints a randomized sweep. Every
//! float is written as the 16 hex digits of its IEEE-754 bits (big-endian), and each result
//! carries the CRC-32 of the full evaluation trace (every point the solver evaluated).

use bobyqa::{Bobyqa, Config, RestartConfig};
use std::fmt::Write as _;

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Rosenbrock,
    Quad,
    Quartic,
    L1,
    Step,
    AllInf,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Rosenbrock => "rosenbrock",
            Kind::Quad => "quad",
            Kind::Quartic => "quartic",
            Kind::L1 => "l1",
            Kind::Step => "step",
            Kind::AllInf => "allinf",
        }
    }
}

#[derive(Clone)]
struct Problem {
    name: String,
    n: usize,
    kind: Kind,
    c: f64,
    g: Vec<f64>,
    h: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    // 0: none, 1: +inf when sum(x) > thr, 2: NaN when sum(x) > thr
    barrier: u8,
    thr: f64,
    x0: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    npt: usize,
    rho_begin: f64,
    rho_end: f64,
    max_fun: usize,
    f_target: f64,
    ftol_rel: Option<f64>,
    ftol_abs: Option<f64>,
    restart: Option<(usize, f64, usize, f64)>,
    prima_parity: bool,
}

/// The objective, written with the exact operation order the PHP test mirrors.
fn objective(p: &Problem, x: &[f64]) -> f64 {
    let n = p.n;
    if p.barrier != 0 {
        let mut s = 0.0;
        for i in 0..n {
            s += x[i];
        }
        if s > p.thr {
            return if p.barrier == 1 { f64::INFINITY } else { f64::NAN };
        }
    }
    let quad = |x: &[f64]| {
        let mut f = p.c;
        for i in 0..n {
            f += p.g[i] * x[i];
        }
        for i in 0..n {
            for j in 0..n {
                f += 0.5 * p.h[i * n + j] * x[i] * x[j];
            }
        }
        f
    };
    match p.kind {
        Kind::Rosenbrock => {
            if n == 1 {
                let b = 1.0 - x[0];
                return b * b;
            }
            let mut f = 0.0;
            for i in 0..n - 1 {
                let a = x[i + 1] - x[i] * x[i];
                let b = 1.0 - x[i];
                f += 100.0 * a * a + b * b;
            }
            f
        }
        Kind::Quad => quad(x),
        Kind::Quartic => {
            let mut f = quad(x);
            for i in 0..n {
                let t = x[i] - p.b[i];
                let t2 = t * t;
                f += p.a[i] * t2 * t2;
            }
            f
        }
        Kind::L1 => {
            let mut f = quad(x);
            for i in 0..n {
                f += p.a[i] * (x[i] - p.b[i]).abs();
            }
            f
        }
        Kind::Step => {
            // Quantized landscape: flat steps of height 1/a[0].
            let q = quad(x) * p.a[0];
            q.floor() / p.a[0]
        }
        Kind::AllInf => f64::INFINITY,
    }
}

fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (0xEDB8_8320 & mask);
        }
    }
    !crc
}

fn hex(v: f64) -> String {
    format!("{:016x}", v.to_bits())
}

fn hex_list(v: &[f64]) -> String {
    let parts: Vec<String> = v.iter().map(|&x| format!("\"{}\"", hex(x))).collect();
    format!("[{}]", parts.join(","))
}

fn opt_hex(v: Option<f64>) -> String {
    v.map_or("null".to_string(), |x| format!("\"{}\"", hex(x)))
}

fn run(p: &Problem, with_trace: bool) -> String {
    let mut config = Config::new(p.n);
    config.npt = p.npt;
    config.rho_begin = p.rho_begin;
    config.rho_end = p.rho_end;
    config.max_fun = p.max_fun;
    config.f_target = p.f_target;
    config.ftol_rel = p.ftol_rel;
    config.ftol_abs = p.ftol_abs;
    config.prima_parity = p.prima_parity;
    if let Some((max_restarts, improve, stall, frac)) = p.restart {
        let mut r = RestartConfig::new();
        r.max_restarts = max_restarts;
        r.improve_rel_tol = improve;
        r.stall_reductions = stall;
        r.cycle_budget_frac = frac;
        config.restart = Some(r);
    }

    let mut spec = String::new();
    write!(
        spec,
        "\"name\":\"{}\",\"n\":{},\"kind\":\"{}\",\"c\":\"{}\",\"g\":{},\"h\":{},\"a\":{},\"b\":{},\
         \"barrier\":{},\"thr\":\"{}\",\"x0\":{},\"lower\":{},\"upper\":{},\"npt\":{},\
         \"rho_begin\":\"{}\",\"rho_end\":\"{}\",\"max_fun\":{},\"f_target\":\"{}\",\
         \"ftol_rel\":{},\"ftol_abs\":{},\"restart\":{},\"prima_parity\":{}",
        p.name,
        p.n,
        p.kind.name(),
        hex(p.c),
        hex_list(&p.g),
        hex_list(&p.h),
        hex_list(&p.a),
        hex_list(&p.b),
        p.barrier,
        hex(p.thr),
        hex_list(&p.x0),
        hex_list(&p.lower),
        hex_list(&p.upper),
        p.npt,
        hex(p.rho_begin),
        hex(p.rho_end),
        p.max_fun,
        hex(p.f_target),
        opt_hex(p.ftol_rel),
        opt_hex(p.ftol_abs),
        p.restart.map_or("null".to_string(), |(m, i, s, f)| format!(
            "{{\"max_restarts\":{m},\"improve_rel_tol\":\"{}\",\"stall_reductions\":{s},\"cycle_budget_frac\":\"{}\"}}",
            hex(i),
            hex(f)
        )),
        p.prima_parity,
    )
    .unwrap();

    let result = match Bobyqa::new(p.n, config) {
        Err(e) => format!("\"error\":\"{e}\""),
        Ok(mut solver) => {
            let mut trace = String::new();
            let mut x = p.x0.clone();
            let outcome = solver.minimize(
                |v: &[f64]| {
                    for (i, &xi) in v.iter().enumerate() {
                        if i > 0 {
                            trace.push(',');
                        }
                        trace.push_str(&hex(xi));
                    }
                    trace.push(';');
                    objective(p, v)
                },
                &mut x,
                &p.lower,
                &p.upper,
            );
            let mut r = format!(
                "\"status\":\"{:?}\",\"display\":\"{}\",\"n_eval\":{},\"f\":\"{}\",\"found_finite\":{},\
                 \"x\":{},\"restarts\":{},\"trace_crc\":{}",
                outcome.status,
                outcome.status,
                outcome.n_eval,
                hex(outcome.f),
                outcome.found_finite(),
                hex_list(&x),
                solver.last_restart_count(),
                crc32(trace.as_bytes()),
            );
            if with_trace {
                write!(r, ",\"trace\":\"{trace}\"").unwrap();
            }
            r
        }
    };
    format!("{{{spec},{result}}}")
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }
}

fn base(name: &str, n: usize, kind: Kind) -> Problem {
    Problem {
        name: name.to_string(),
        n,
        kind,
        c: 0.0,
        g: vec![0.0; n],
        h: vec![0.0; n * n],
        a: vec![0.0; n],
        b: vec![0.0; n],
        barrier: 0,
        thr: 0.0,
        x0: vec![0.0; n],
        lower: vec![-10.0; n],
        upper: vec![10.0; n],
        npt: 2 * n + 1,
        rho_begin: 1.0,
        rho_end: 1e-6,
        max_fun: 500 * n,
        f_target: f64::NEG_INFINITY,
        ftol_rel: None,
        ftol_abs: None,
        restart: None,
        prima_parity: true,
    }
}

/// Positive definite (or, with `indefinite`, possibly indefinite) symmetric matrix.
fn random_h(rng: &mut Rng, n: usize, indefinite: bool) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for v in m.iter_mut() {
        *v = rng.range(-1.0, 1.0);
    }
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..n {
                s += m[i * n + k] * m[j * n + k];
            }
            h[i * n + j] = s;
        }
        h[i * n + i] += rng.range(0.05, 2.0);
    }
    if indefinite {
        let i = rng.below(n);
        h[i * n + i] -= rng.range(0.0, 4.0);
    }
    h
}

fn random_problem(rng: &mut Rng, index: usize) -> Problem {
    let n = 1 + rng.below(8);
    let kind = match rng.below(10) {
        0 | 1 => Kind::Rosenbrock,
        2..=4 => Kind::Quad,
        5..=7 => Kind::Quartic,
        8 => Kind::L1,
        _ => {
            if rng.chance(0.4) {
                Kind::Step
            } else if rng.chance(0.3) {
                Kind::AllInf
            } else {
                Kind::Quad
            }
        }
    };
    let mut p = base(&format!("sweep-{index}"), n, kind);
    let scale = if rng.chance(0.08) { 1e13 } else if rng.chance(0.08) { 1e-4 } else if rng.chance(0.03) { 1e150 } else if rng.chance(0.03) { 1e-150 } else { 1.0 };
    p.c = rng.range(-5.0, 5.0) * scale;
    for i in 0..n {
        p.g[i] = rng.range(-3.0, 3.0) * scale;
        p.a[i] = rng.range(0.01, 3.0) * scale;
        p.b[i] = rng.range(-2.0, 2.0);
    }
    let indefinite = kind == Kind::Quartic && rng.chance(0.4);
    p.h = random_h(rng, n, indefinite);
    for v in p.h.iter_mut() {
        *v *= scale;
    }
    if kind == Kind::L1 {
        for v in p.h.iter_mut() {
            *v *= 0.05;
        }
    }
    if kind == Kind::Step {
        p.a[0] = [1.0, 10.0, 1000.0, 1e6][rng.below(4)] / scale;
    }

    // Configuration: the extension's shape most of the time, crate defaults otherwise.
    let ext_like = rng.chance(0.6);
    p.rho_begin = if ext_like {
        0.1
    } else {
        [1.0, 0.5, 0.25, 2.0][rng.below(4)]
    };
    let rho_end: f64 = [1e-4, 1e-6, 1e-8, 1e-10, 1e-3][rng.below(5)];
    p.rho_end = rho_end.min(p.rho_begin);
    let max_npt = (n + 1) * (n + 2) / 2;
    p.npt = match rng.below(4) {
        0 => n + 2,
        1 => max_npt,
        2 => n + 2 + rng.below(max_npt - n - 1),
        _ => (2 * n + 1).min(max_npt),
    };
    p.max_fun = match rng.below(5) {
        0 => p.npt + 1 + rng.below(3 * n + 3),
        1 => p.npt + 1 + rng.below(60),
        2 => 100 * n,
        _ => 500 * n,
    };
    p.max_fun = p.max_fun.max(p.npt + 1);

    // Bounds: width at least 2 * rho_begin, sometimes infinite, sometimes tight near x0.
    let rb = p.rho_begin;
    for i in 0..n {
        let center = rng.range(-2.0, 2.0);
        let lo = match rng.below(6) {
            0 => f64::NEG_INFINITY,
            1 => center - 2.0 * rb - rng.range(0.0, 0.5),
            _ => center - rng.range(2.0 * rb, 5.0),
        };
        let hi = match rng.below(6) {
            0 => f64::INFINITY,
            1 => center + rng.range(0.0, 0.5),
            _ => center + rng.range(0.1, 5.0),
        };
        let (lo, hi) = if hi - lo < 2.0 * rb {
            (lo, lo + 2.0 * rb + rng.range(0.0, 1.0))
        } else {
            (lo, hi)
        };
        p.lower[i] = lo;
        p.upper[i] = hi;
        let flo = if lo.is_finite() { lo } else { center - 6.0 };
        let fhi = if hi.is_finite() { hi } else { center + 6.0 };
        p.x0[i] = match rng.below(8) {
            0 => flo,
            1 => fhi,
            2 => flo + rng.range(0.0, rb),
            3 => fhi + 3.0, // outside the box: clamped
            _ => rng.range(flo, fhi),
        };
    }
    if kind != Kind::AllInf && rng.chance(0.15) {
        p.barrier = if rng.chance(0.7) { 1 } else { 2 };
        let mut s = 0.0;
        for i in 0..n {
            s += p.x0[i].max(p.lower[i]).min(p.upper[i]);
        }
        p.thr = s + rng.range(0.0, 1.5);
    }
    if rng.chance(0.05) {
        p.f_target = p.c - rng.range(0.0, 2.0) * scale.max(1.0);
    }
    if rng.chance(0.08) {
        p.ftol_rel = Some([0.0, 1e-8, 1e-4][rng.below(3)]);
        if rng.chance(0.5) {
            p.ftol_abs = Some(1e-6);
        }
    }
    if rng.chance(0.1) {
        p.restart = Some((
            1 + rng.below(3),
            [1e-6, 0.0, 1e-3][rng.below(3)],
            rng.below(3),
            [0.125, 0.0, 0.5][rng.below(3)],
        ));
    }
    if rng.chance(0.1) {
        p.prima_parity = false;
    }
    p
}

fn fixed_problems() -> Vec<Problem> {
    let mut out = Vec::new();

    let mut p = base("sphere-default", 2, Kind::Quad);
    p.h = vec![2.0, 0.0, 0.0, 2.0];
    p.x0 = vec![1.0, 2.0];
    p.lower = vec![-5.0, -5.0];
    p.upper = vec![5.0, 5.0];
    out.push(p);

    let mut p = base("rosenbrock-2-default", 2, Kind::Rosenbrock);
    p.x0 = vec![-1.2, 1.0];
    p.lower = vec![-5.0, -5.0];
    p.upper = vec![5.0, 5.0];
    out.push(p);

    let mut p = base("rosenbrock-4-extension", 4, Kind::Rosenbrock);
    p.x0 = vec![0.5, 0.5, 0.5, 0.5];
    p.lower = vec![f64::NEG_INFINITY; 4];
    p.upper = vec![f64::INFINITY; 4];
    p.rho_begin = 0.1;
    p.rho_end = 1e-8;
    p.max_fun = 10_000;
    out.push(p);

    // Variance-component shape: theta >= 0 with the optimum on the bound for theta[1].
    let mut p = base("quad-3-bound-active", 3, Kind::Quad);
    p.g = vec![-1.0, 2.0, -0.5];
    p.h = vec![2.0, 0.3, 0.1, 0.3, 1.5, -0.2, 0.1, -0.2, 1.0];
    p.x0 = vec![1.0, 1.0, 0.0];
    p.lower = vec![0.0, 0.0, f64::NEG_INFINITY];
    p.upper = vec![f64::INFINITY; 3];
    p.rho_begin = 0.1;
    p.rho_end = 1e-10;
    p.max_fun = 10_000;
    out.push(p);

    let mut p = base("quartic-5-infinite-region", 5, Kind::Quartic);
    p.g = vec![0.1, -0.2, 0.3, -0.1, 0.05];
    p.h = vec![
        1.0, 0.1, 0.0, 0.0, 0.0, 0.1, 1.0, 0.1, 0.0, 0.0, 0.0, 0.1, 1.0, 0.1, 0.0, 0.0, 0.0, 0.1,
        1.0, 0.1, 0.0, 0.0, 0.0, 0.1, 1.0,
    ];
    p.a = vec![0.5, 1.0, 0.25, 2.0, 0.75];
    p.b = vec![1.0, -0.5, 0.25, 0.5, -1.0];
    p.barrier = 1;
    p.thr = 0.7;
    p.x0 = vec![0.0; 5];
    p.lower = vec![-3.0; 5];
    p.upper = vec![3.0; 5];
    p.rho_begin = 0.1;
    p.rho_end = 1e-6;
    p.max_fun = 2000;
    out.push(p);

    let mut p = base("rosenbrock-4-max-fun", 4, Kind::Rosenbrock);
    p.x0 = vec![-1.0, 1.0, -1.0, 1.0];
    p.lower = vec![-4.0; 4];
    p.upper = vec![4.0; 4];
    p.rho_begin = 0.1;
    p.rho_end = 1e-8;
    p.max_fun = 40;
    out.push(p);

    let mut p = base("all-infinite", 2, Kind::AllInf);
    p.x0 = vec![1.0, 2.0];
    p.lower = vec![-5.0, -5.0];
    p.upper = vec![5.0, 5.0];
    out.push(p);

    let mut p = base("nan-region-target", 3, Kind::Quad);
    p.c = 3.0;
    p.g = vec![-1.0, -1.0, -1.0];
    p.h = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    p.barrier = 2;
    p.thr = 2.5;
    p.x0 = vec![0.2, 0.2, 0.2];
    p.lower = vec![-2.0; 3];
    p.upper = vec![2.0; 3];
    p.rho_begin = 0.1;
    p.rho_end = 1e-8;
    p.f_target = 1.6;
    out.push(p);

    let mut p = base("one-dimensional", 1, Kind::Quartic);
    p.g = vec![0.3];
    p.h = vec![0.5];
    p.a = vec![1.0];
    p.b = vec![0.7];
    p.npt = 3;
    p.x0 = vec![-0.4];
    p.lower = vec![-2.0];
    p.upper = vec![2.0];
    p.rho_begin = 0.1;
    p.rho_end = 1e-8;
    p.max_fun = 400;
    out.push(p);

    let mut p = base("l1-8-full-npt", 8, Kind::L1);
    let mut rng = Rng(42);
    p.h = random_h(&mut rng, 8, false);
    for v in p.h.iter_mut() {
        *v *= 0.05;
    }
    for i in 0..8 {
        p.g[i] = rng.range(-1.0, 1.0);
        p.a[i] = rng.range(0.1, 1.0);
        p.b[i] = rng.range(-1.0, 1.0);
        p.x0[i] = rng.range(-1.0, 1.0);
    }
    p.lower = vec![-3.0; 8];
    p.upper = vec![3.0; 8];
    p.npt = 45;
    p.rho_begin = 0.1;
    p.rho_end = 1e-6;
    p.max_fun = 3000;
    out.push(p);

    let mut p = base("quartic-6-restart-ftol", 6, Kind::Quartic);
    let mut rng = Rng(7);
    p.h = random_h(&mut rng, 6, true);
    for i in 0..6 {
        p.g[i] = rng.range(-1.0, 1.0);
        p.a[i] = rng.range(0.1, 1.0);
        p.b[i] = rng.range(-1.0, 1.0);
    }
    p.lower = vec![-2.0; 6];
    p.upper = vec![2.0; 6];
    p.npt = 10;
    p.rho_begin = 0.5;
    p.rho_end = 1e-8;
    p.max_fun = 3000;
    p.restart = Some((2, 1e-6, 2, 0.125));
    p.ftol_abs = Some(1e-12);
    out.push(p);

    let mut p = base("quad-2-full-model-no-parity", 2, Kind::Quad);
    p.g = vec![0.5, -0.25];
    p.h = vec![1.0, 0.2, 0.2, 3.0];
    p.npt = 6;
    p.lower = vec![-3.0, -3.0];
    p.upper = vec![3.0, 3.0];
    p.prima_parity = false;
    p.rho_end = 1e-9;
    out.push(p);

    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map_or("fixture", String::as_str);
    let mut rows = Vec::new();
    match mode {
        "fixture" => {
            for p in fixed_problems() {
                rows.push(run(&p, false));
            }
        }
        "sweep" => {
            let count: usize = args.get(2).map_or(500, |s| s.parse().unwrap());
            let seed: u64 = args.get(3).map_or(1, |s| s.parse().unwrap());
            let mut rng = Rng(seed);
            for i in 0..count {
                let p = random_problem(&mut rng, i);
                rows.push(run(&p, true));
            }
        }
        "probe" => {
            // Semantics the port depends on: f64::max/min on signed zeros, midpoint, empty sum.
            let z = std::hint::black_box(-0.0f64);
            let pz = std::hint::black_box(0.0f64);
            let empty: [f64; 0] = std::hint::black_box([]);
            println!(
                "max(-0,0)={:x} max(0,-0)={:x} min(-0,0)={:x} min(0,-0)={:x} sum([])={:x} sum([-0])={:x}",
                z.max(pz).to_bits(),
                pz.max(z).to_bits(),
                z.min(pz).to_bits(),
                pz.min(z).to_bits(),
                empty.iter().sum::<f64>().to_bits(),
                std::hint::black_box([-0.0f64]).iter().sum::<f64>().to_bits(),
            );
            let big = std::hint::black_box(f64::MAX);
            let tiny = std::hint::black_box(5e-324f64);
            println!(
                "mid(MAX,MAX)={:x} mid(tiny,MAX)={:x} mid(1,3)={:x} mid(tiny,tiny)={:x}",
                big.midpoint(big).to_bits(),
                tiny.midpoint(big).to_bits(),
                std::hint::black_box(1.0f64).midpoint(3.0).to_bits(),
                tiny.midpoint(tiny).to_bits(),
            );
            return;
        }
        other => panic!("unknown mode {other}"),
    }
    println!("[\n{}\n]", rows.join(",\n"));
}
