//! Reference harness for the PHP port of argmin 0.11.0's Nelder-Mead solver
//! (src/Internal/NelderMead.php). Every float is written as the 16 hex digits
//! of its IEEE-754 bit pattern so the comparison is exact.
//!
//!   cargo run --release --offline -- fixture > ../../fixtures/neldermead-reference.json
//!   cargo run --release --offline -- sweep <count> <seed> > sweep.json

use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::neldermead::NelderMead;
use std::cell::{Cell, RefCell};
use std::fmt::Write as _;
use std::sync::Mutex;

static PANIC: Mutex<String> = Mutex::new(String::new());

#[derive(Clone)]
struct Problem {
    name: String,
    kind: &'static str, // "quad" | "rosen"
    centre: Vec<f64>,
    weights: Vec<f64>,
    simplex: Vec<Vec<f64>>,
    sd_tolerance: f64,
    max_iters: u64,
    budget: Option<usize>,
    lower: Option<Vec<f64>>,
    nan_above: Option<f64>,
    neginf_below: Option<f64>,
    plateau: Option<f64>,
    error_message: String,
}

struct Objective<'a> {
    p: &'a Problem,
    evaluations: &'a Cell<usize>,
    trace: &'a RefCell<Vec<(Vec<f64>, Option<f64>)>>,
}

fn base_cost(p: &Problem, x: &[f64]) -> f64 {
    let x: Vec<f64> = match p.plateau {
        Some(s) => x.iter().map(|v| (v / s).floor() * s).collect(),
        None => x.to_vec(),
    };
    let n = x.len();
    let mut acc = 0.0;
    if p.kind == "quad" {
        for i in 0..n {
            let d = x[i] - p.centre[i];
            acc += p.weights[i] * d * d;
        }
    } else {
        for i in 0..n - 1 {
            let a = x[i + 1] - x[i] * x[i];
            let b = 1.0 - x[i];
            acc += 100.0 * a * a + b * b;
        }
        let b = p.centre[0] - x[n - 1];
        acc += b * b;
    }
    acc
}

impl CostFunction for Objective<'_> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, x: &Vec<f64>) -> Result<f64, Error> {
        let p = self.p;
        if let Some(budget) = p.budget {
            if self.evaluations.get() >= budget {
                self.trace.borrow_mut().push((x.clone(), None));
                return Err(Error::msg(p.error_message.clone()));
            }
        }
        self.evaluations.set(self.evaluations.get() + 1);
        let c = if p
            .lower
            .as_ref()
            .is_some_and(|l| x.iter().zip(l).any(|(v, l)| !v.is_finite() || *v < *l || *v > 1e4))
        {
            f64::INFINITY
        } else if p.nan_above.is_some_and(|t| x[0] > t) {
            f64::NAN
        } else if p.neginf_below.is_some_and(|t| x[0] < t) {
            f64::NEG_INFINITY
        } else {
            base_cost(p, x)
        };
        self.trace.borrow_mut().push((x.clone(), Some(c)));
        Ok(c)
    }
}

fn hex(v: f64) -> String {
    format!("\"{:016x}\"", v.to_bits())
}

fn hexes(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|x| hex(*x)).collect::<Vec<_>>().join(","))
}

fn opt_hex(v: Option<f64>) -> String {
    v.map_or("null".into(), hex)
}

fn json_str(s: &str) -> String {
    let mut out = String::from("\"");
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            c if (c as u32) < 0x20 => write!(out, "\\u{:04x}", c as u32).unwrap(),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn run(p: &Problem, trace_limit: usize) -> String {
    let evaluations = Cell::new(0);
    let trace = RefCell::new(Vec::new());
    let mut spec = String::new();
    write!(
        spec,
        "{{\"kind\":\"{}\",\"centre\":{},\"weights\":{},\"simplex\":[{}],\"sd_tolerance\":{},\"max_iters\":{},\"budget\":{},\"lower\":{},\"nan_above\":{},\"neginf_below\":{},\"plateau\":{},\"error_message\":{}}}",
        p.kind,
        hexes(&p.centre),
        hexes(&p.weights),
        p.simplex.iter().map(|v| hexes(v)).collect::<Vec<_>>().join(","),
        hex(p.sd_tolerance),
        if p.max_iters == u64::MAX { "null".to_string() } else { p.max_iters.to_string() },
        p.budget.map_or("null".into(), |b| b.to_string()),
        p.lower.as_ref().map_or("null".into(), |l| hexes(l)),
        opt_hex(p.nan_above),
        opt_hex(p.neginf_below),
        opt_hex(p.plateau),
        json_str(&p.error_message),
    )
    .unwrap();

    let result = match NelderMead::new(p.simplex.clone()).with_sd_tolerance(p.sd_tolerance) {
        Err(e) => format!("{{\"status\":\"invalid\",\"message\":{}}}", json_str(&e.to_string())),
        Ok(solver) => {
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                Executor::new(Objective { p, evaluations: &evaluations, trace: &trace }, solver)
                    .configure(|state| state.max_iters(p.max_iters))
                    .run()
            }));
            match outcome {
                Err(_) => format!(
                    "{{\"status\":\"panic\",\"message\":{}}}",
                    json_str(&PANIC.lock().unwrap())
                ),
                Ok(Err(e)) => format!("{{\"status\":\"error\",\"message\":{}}}", json_str(&e.to_string())),
                Ok(Ok(r)) => {
                    let s = &r.state;
                    format!(
                        "{{\"status\":\"ok\",\"termination\":{},\"iterations\":{},\"best_param\":{},\"best_cost\":{}}}",
                        s.get_termination_reason().map_or("null".into(), |t| json_str(&t.to_string())),
                        s.get_iter(),
                        s.get_best_param().map_or("null".into(), |b| hexes(b)),
                        hex(s.get_best_cost()),
                    )
                }
            }
        }
    };
    let trace = trace
        .borrow()
        .iter()
        .take(trace_limit)
        .map(|(x, c)| format!("[{},{}]", hexes(x), opt_hex(*c)))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"name\":{},\"problem\":{},\"result\":{},\"evaluations\":{},\"trace\":[{}]}}",
        json_str(&p.name),
        spec,
        result,
        evaluations.get(),
        trace
    )
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.unit()
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }
}

fn lmm_simplex(start: &[f64]) -> Vec<Vec<f64>> {
    let mut simplex = vec![start.to_vec()];
    for j in 0..start.len() {
        let mut v = start.to_vec();
        v[j] = if v[j] + 0.1 > 1e4 { v[j] - 0.1 } else { v[j] + 0.1 };
        simplex.push(v);
    }
    simplex
}

fn base(name: &str, kind: &'static str, n: usize) -> Problem {
    Problem {
        name: name.into(),
        kind,
        centre: (0..n).map(|i| 0.5 + i as f64 * 0.25).collect(),
        weights: (0..n).map(|i| 1.0 + i as f64).collect(),
        simplex: lmm_simplex(&vec![1.0; n]),
        sd_tolerance: 1e-8,
        max_iters: 2000,
        budget: None,
        lower: None,
        nan_above: None,
        neginf_below: None,
        plateau: None,
        error_message: "Nelder-Mead evaluation budget exhausted".into(),
    }
}

fn fixture() -> Vec<Problem> {
    let mut v = Vec::new();
    for n in 1..=6 {
        v.push(base(&format!("quad n={n}"), "quad", n));
        let mut r = base(&format!("rosen n={n}"), "rosen", n);
        r.simplex = lmm_simplex(&vec![-0.5; n]);
        r.centre = vec![1.0; n];
        v.push(r);
    }
    let mut p = base("lmm-like bounds and budget", "quad", 3);
    p.lower = Some(vec![0.0, f64::NEG_INFINITY, 0.0]);
    p.centre = vec![-0.2, 0.3, 0.05];
    p.simplex = lmm_simplex(&[0.0, 0.0, 1.0]);
    p.budget = Some(400);
    v.push(p);
    let mut p = base("budget error mid-run", "rosen", 2);
    p.budget = Some(25);
    v.push(p);
    let mut p = base("budget error inside shrink", "quad", 2);
    p.plateau = Some(1.0);
    p.sd_tolerance = 0.0;
    p.budget = Some(6);
    v.push(p);
    let mut p = base("budget shorter than the simplex panics", "quad", 3);
    p.budget = Some(2);
    v.push(p);
    let mut p = base("max_iters 0", "quad", 2);
    p.max_iters = 0;
    v.push(p);
    let mut p = base("max_iters 7", "rosen", 3);
    p.max_iters = 7;
    v.push(p);
    let mut p = base("infinite start", "quad", 2);
    p.lower = Some(vec![2.0, 2.0]);
    p.max_iters = 50;
    v.push(p);
    let mut p = base("all infinite", "quad", 2);
    p.lower = Some(vec![1e5, 1e5]);
    p.max_iters = 20;
    v.push(p);
    let mut p = base("NaN reflection is an unreachable point", "quad", 2);
    p.nan_above = Some(1.05);
    v.push(p);
    let mut p = base("negative infinity reaches the target cost", "quad", 2);
    p.neginf_below = Some(0.8);
    v.push(p);
    let mut p = base("plateau ties", "quad", 3);
    p.plateau = Some(0.25);
    p.sd_tolerance = 0.0;
    p.max_iters = 60;
    v.push(p);
    let mut p = base("near the 1e4 bound", "quad", 2);
    p.simplex = lmm_simplex(&[9999.95, 3.0]);
    p.centre = vec![9990.0, 0.0];
    v.push(p);
    let mut p = base("sd_tolerance NaN never converges", "quad", 2);
    p.sd_tolerance = f64::NAN;
    p.max_iters = 40;
    v.push(p);
    let mut p = base("sd_tolerance -0 is accepted", "quad", 1);
    p.sd_tolerance = -0.0;
    p.max_iters = 300;
    v.push(p);
    let mut p = base("sd_tolerance negative is rejected", "quad", 1);
    p.sd_tolerance = -1e-300;
    v.push(p);
    let mut p = base("sd_tolerance -inf is rejected", "quad", 1);
    p.sd_tolerance = f64::NEG_INFINITY;
    v.push(p);
    let mut p = base("unbounded iterations", "quad", 2);
    p.max_iters = u64::MAX;
    v.push(p);
    let mut p = base("cost error text is passed through", "quad", 1);
    p.budget = Some(3);
    p.error_message = "custom \"quoted\" failure".into();
    v.push(p);
    let mut p = base("mismatched vertex lengths panic", "quad", 3);
    p.simplex = vec![vec![1.0, 1.0, 1.0], vec![1.1, 1.0], vec![1.0, 1.1, 1.0], vec![1.0, 1.0, 1.1]];
    v.push(p);
    let mut p = base("zero-dimensional vertices panic", "quad", 1);
    p.sd_tolerance = 0.0;
    p.simplex = vec![vec![], vec![]];
    v.push(p);
    let mut p = base("empty simplex panics", "quad", 1);
    p.simplex = vec![];
    v.push(p);
    let mut p = base("single vertex panics", "quad", 1);
    p.simplex = vec![vec![1.0]];
    v.push(p);
    let mut p = base("single vertex with max_iters 0", "quad", 1);
    p.simplex = vec![vec![1.0]];
    p.max_iters = 0;
    v.push(p);
    let mut rng = Rng(2);
    let mut p = (0..=261).map(|i| random_problem(&mut rng, i)).last().unwrap();
    p.name = "NaN best vertex leaves no best param (sweep 2 #261)".into();
    v.push(p);
    let mut rng = Rng(7);
    v.extend((0..40).map(|i| random_problem(&mut rng, i)).filter(|p| p.simplex.len() <= 7));
    v
}

fn random_problem(rng: &mut Rng, id: usize) -> Problem {
    let n = if rng.chance(0.05) { 18 + rng.below(6) as usize } else { 1 + rng.below(6) as usize };
    let kind = if rng.chance(0.6) { "quad" } else { "rosen" };
    let mut p = base(&format!("random #{id}"), kind, n);
    p.centre = (0..n).map(|_| rng.range(-2.0, 2.0)).collect();
    p.weights = (0..n).map(|_| rng.range(0.1, 10.0)).collect();
    let start: Vec<f64> = (0..n).map(|_| rng.range(-3.0, 3.0)).collect();
    p.simplex = match rng.below(4) {
        0 | 1 => lmm_simplex(&start),
        2 => (0..=n).map(|_| (0..n).map(|_| rng.range(-3.0, 3.0)).collect()).collect(),
        _ => {
            // degenerate: a repeated vertex
            let mut s = lmm_simplex(&start);
            s[n] = s[0].clone();
            s
        }
    };
    p.sd_tolerance = [0.0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-3, f64::EPSILON, 1e-8][rng.below(8) as usize];
    p.max_iters = [0, 1, 2, 3, 10, 50, 200, 1000, 5000][rng.below(9) as usize];
    if rng.chance(0.25) {
        p.budget = Some(rng.below(400) as usize);
    }
    if rng.chance(0.2) {
        p.lower = Some((0..n).map(|i| if rng.chance(0.5) { f64::NEG_INFINITY } else { start[i] - rng.range(0.0, 1.0) }).collect());
    }
    if rng.chance(0.1) {
        p.nan_above = Some(start[0] + rng.range(-0.05, 0.5));
    }
    if rng.chance(0.05) {
        p.neginf_below = Some(start[0] - rng.range(0.0, 1.0));
    }
    if rng.chance(0.15) {
        p.plateau = Some([0.5, 0.25, 0.1, 1.0][rng.below(4) as usize]);
    }
    if rng.chance(0.03) {
        p.sd_tolerance = f64::NAN;
    }
    if rng.chance(0.02) {
        p.sd_tolerance = -1e-9;
    }
    p
}

fn main() {
    std::panic::set_hook(Box::new(|info| {
        let msg = if let Some(s) = info.payload().downcast_ref::<&str>() {
            s.to_string()
        } else if let Some(s) = info.payload().downcast_ref::<String>() {
            s.clone()
        } else {
            String::from("Box<dyn Any>")
        };
        *PANIC.lock().unwrap() = msg;
    }));
    let args: Vec<String> = std::env::args().collect();
    let (problems, trace_limit) = match args.get(1).map(String::as_str) {
        Some("sweep") => {
            let count: usize = args[2].parse().unwrap();
            let mut rng = Rng(args[3].parse().unwrap());
            ((0..count).map(|i| random_problem(&mut rng, i)).collect::<Vec<_>>(), usize::MAX)
        }
        // The committed fixture keeps the first 120 evaluated points of each run.
        _ => (fixture(), 120),
    };
    let body: Vec<String> = problems.iter().map(|p| run(p, trace_limit)).collect();
    println!(
        "{{\"source\":\"argmin 0.11.0 NelderMead + Executor, argmin-math 0.5.1 (vec), tests/reference/neldermead\",\"problems\":[\n{}\n]}}",
        body.join(",\n")
    );
}
