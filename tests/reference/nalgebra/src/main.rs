//! Reference outputs from nalgebra 0.35.0 (matrixmultiply 0.3.11) for the
//! polyfill's Dense, Svd and SymmetricEigen ports. Every float is printed as
//! its IEEE-754 bit pattern in hex; matrices are column-major.
//!
//! Usage: nalgebra-ref fixture | sweep <seed>

use nalgebra::{DMatrix, DVector, SymmetricEigen};
use std::fmt::Write;

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn unif(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn sym(&mut self) -> f64 {
        2.0 * self.unif() - 1.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    fn normal(&mut self) -> f64 {
        let u = self.unif().max(1e-300);
        let v = self.unif();
        (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos()
    }
}

const KINDS: usize = 14;

fn generate(rng: &mut Rng, r: usize, c: usize, kind: usize) -> DMatrix<f64> {
    match kind {
        0 => DMatrix::from_fn(r, c, |_, _| rng.sym()),
        1 => DMatrix::from_fn(r, c, |_, _| rng.normal()),
        2 => {
            // badly scaled columns
            let scales: Vec<f64> = (0..c).map(|_| 10f64.powi(rng.below(17) as i32 - 8)).collect();
            DMatrix::from_fn(r, c, |_, j| rng.normal() * scales[j])
        }
        3 => DMatrix::from_fn(r, c, |_, _| if rng.unif() < 0.4 { 0.0 } else { rng.normal() }),
        4 => {
            // rank deficient: product of r x k and k x c with k < min
            let k = rng.below(r.min(c).max(1));
            let a = DMatrix::from_fn(r, k, |_, _| rng.normal());
            let b = DMatrix::from_fn(k, c, |_, _| rng.normal());
            if k == 0 { DMatrix::zeros(r, c) } else { a * b }
        }
        5 => DMatrix::from_fn(r, c, |_, _| (rng.below(7) as f64) - 3.0),
        6 => {
            // duplicate columns and rows
            let mut m = DMatrix::from_fn(r, c, |_, _| rng.normal());
            if c > 1 {
                let j = rng.below(c - 1) + 1;
                let col = m.column(0).into_owned();
                m.set_column(j, &col);
            }
            if r > 1 {
                let i = rng.below(r - 1) + 1;
                let row = m.row(0).into_owned();
                m.set_row(i, &row);
            }
            m
        }
        7 => DMatrix::zeros(r, c),
        8 => {
            // extreme magnitudes
            let e = if rng.unif() < 0.5 { 1e-200 } else { 1e150 };
            DMatrix::from_fn(r, c, |_, _| rng.normal() * e)
        }
        9 => DMatrix::from_fn(r, c, |_, _| match rng.below(4) {
            0 => -0.0,
            1 => 0.0,
            _ => rng.normal(),
        }),
        10 => {
            // nearly rank deficient
            let k = r.min(c).saturating_sub(1).max(1);
            let a = DMatrix::from_fn(r, k, |_, _| rng.normal());
            let b = DMatrix::from_fn(k, c, |_, _| rng.normal());
            let noise = DMatrix::from_fn(r, c, |_, _| rng.normal() * 1e-13);
            a * b + noise
        }
        11 => {
            // upper bidiagonal-ish / diagonal structure
            DMatrix::from_fn(r, c, |i, j| if i == j || i + 1 == j { rng.normal() } else { 0.0 })
        }
        12 => {
            // a single nonzero row or column
            let row = rng.unif() < 0.5;
            let k = if row { rng.below(r) } else { rng.below(c) };
            DMatrix::from_fn(r, c, |i, j| if (row && i == k) || (!row && j == k) { rng.normal() } else { 0.0 })
        }
        _ => {
            // repeated singular values: scaled orthogonal-ish integer pattern
            DMatrix::from_fn(r, c, |i, j| if (i + j) % 2 == 0 { 1.0 } else { -1.0 })
        }
    }
}

fn h(x: f64) -> String {
    format!("\"{:016x}\"", x.to_bits())
}

fn arr(xs: &[f64]) -> String {
    let parts: Vec<String> = xs.iter().map(|x| h(*x)).collect();
    format!("[{}]", parts.join(","))
}

struct Out {
    svd: Vec<String>,
    solve: Vec<String>,
    eigen: Vec<String>,
    mul: Vec<String>,
    norm: Vec<String>,
    hypot: Vec<String>,
    fma: Vec<String>,
}

fn svd_case(out: &mut Out, a: &DMatrix<f64>, cu: bool, cv: bool) {
    let svd = a.clone().svd(cu, cv);
    let mut s = String::new();
    write!(
        s,
        "{{\"rows\":{},\"cols\":{},\"u\":{},\"v\":{},\"a\":{},\"s\":{},\"U\":{},\"Vt\":{}}}",
        a.nrows(),
        a.ncols(),
        cu,
        cv,
        arr(a.as_slice()),
        arr(svd.singular_values.as_slice()),
        svd.u.as_ref().map_or("null".to_string(), |m| arr(m.as_slice())),
        svd.v_t.as_ref().map_or("null".to_string(), |m| arr(m.as_slice())),
    )
    .unwrap();
    out.svd.push(s);
}

fn solve_case(out: &mut Out, a: &DMatrix<f64>, b: &DMatrix<f64>, eps: f64) {
    let svd = a.clone().svd(true, true);
    let x = svd.solve(b, eps).unwrap();
    let mut s = String::new();
    write!(
        s,
        "{{\"rows\":{},\"cols\":{},\"bcols\":{},\"a\":{},\"b\":{},\"eps\":{},\"x\":{}}}",
        a.nrows(),
        a.ncols(),
        b.ncols(),
        arr(a.as_slice()),
        arr(b.as_slice()),
        h(eps),
        arr(x.as_slice()),
    )
    .unwrap();
    out.solve.push(s);
}

fn eigen_case(out: &mut Out, a: &DMatrix<f64>) {
    let eig = SymmetricEigen::new(a.clone());
    let mut s = String::new();
    write!(
        s,
        "{{\"n\":{},\"a\":{},\"values\":{},\"vectors\":{}}}",
        a.nrows(),
        arr(a.as_slice()),
        arr(eig.eigenvalues.as_slice()),
        arr(eig.eigenvectors.as_slice()),
    )
    .unwrap();
    out.eigen.push(s);
}

fn mul_case(out: &mut Out, op: &str, a: &DMatrix<f64>, b: &DMatrix<f64>) {
    let c = match op {
        "mul" => a * b,
        "mul_bt" => a * b.transpose(),
        "tr_mul" => a.tr_mul(b),
        "mul_vec" => {
            let x = DVector::from_column_slice(b.as_slice());
            DMatrix::from_column_slice(a.nrows(), 1, (a * x).as_slice())
        }
        "tr_mul_vec" => {
            let x = DVector::from_column_slice(b.as_slice());
            DMatrix::from_column_slice(a.ncols(), 1, a.tr_mul(&x).as_slice())
        }
        _ => unreachable!(),
    };
    let mut s = String::new();
    write!(
        s,
        "{{\"op\":\"{}\",\"ar\":{},\"ac\":{},\"br\":{},\"bc\":{},\"a\":{},\"b\":{},\"c\":{}}}",
        op,
        a.nrows(),
        a.ncols(),
        b.nrows(),
        b.ncols(),
        arr(a.as_slice()),
        arr(b.as_slice()),
        arr(c.as_slice()),
    )
    .unwrap();
    out.mul.push(s);
}

fn norm_case(out: &mut Out, a: &DMatrix<f64>) {
    let mut s = String::new();
    let v = DVector::from_column_slice(a.as_slice());
    write!(
        s,
        "{{\"rows\":{},\"cols\":{},\"a\":{},\"norm_squared\":{},\"norm\":{},\"vector_norm_squared\":{},\"vector_norm\":{}}}",
        a.nrows(),
        a.ncols(),
        arr(a.as_slice()),
        h(a.norm_squared()),
        h(a.norm()),
        h(v.norm_squared()),
        h(v.norm()),
    )
    .unwrap();
    out.norm.push(s);
}

fn jacobi(points: usize) -> DMatrix<f64> {
    // logistic_mixed.rs normal_quadrature
    let mut m = DMatrix::zeros(points, points);
    for j in 1..points {
        m[(j - 1, j)] = (j as f64).sqrt();
        m[(j, j - 1)] = (j as f64).sqrt();
    }
    m
}

fn symmetric(rng: &mut Rng, n: usize, kind: usize) -> DMatrix<f64> {
    let a = generate(rng, n, n, kind);
    let mut s = &a + a.transpose();
    if kind == 13 {
        s = DMatrix::identity(n, n) * 2.0;
    }
    s
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("fixture");
    let seed: u64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(20251016);
    let fixture = mode == "fixture";
    let mut rng = Rng(seed);
    let mut out = Out { svd: vec![], solve: vec![], eigen: vec![], mul: vec![], norm: vec![], hypot: vec![], fma: vec![] };

    // SVD over shapes and kinds.
    let max_dim = if fixture { 6 } else { 12 };
    for r in 1..=max_dim {
        for c in 1..=max_dim {
            let kinds: Vec<usize> = if fixture {
                vec![rng.below(KINDS)]
            } else {
                (0..KINDS).collect()
            };
            for kind in kinds {
                let a = generate(&mut rng, r, c, kind);
                let flags = if fixture { rng.below(4) } else { (r * 7 + c * 3 + kind) % 4 };
                svd_case(&mut out, &a, flags & 1 == 1 || !fixture && flags == 0, flags & 2 == 2 || !fixture && flags == 0);
                if !fixture {
                    svd_case(&mut out, &a, false, false);
                }
            }
        }
    }
    if fixture {
        // More than 32 singular values with exact zero ties: Rust's quicksort orders them.
        for kind in [4, 12] {
            let a = generate(&mut rng, 36, 34, kind);
            svd_case(&mut out, &a, true, true);
        }
    }
    let tall: &[(usize, usize)] = if fixture { &[(40, 6), (6, 30)] } else { &[(200, 6), (6, 200), (100, 3), (50, 8), (64, 20), (30, 25), (25, 30), (45, 35), (35, 45), (24, 22), (40, 40)] };
    for &(r, c) in tall {
        let kinds: Vec<usize> = if fixture { vec![0, 2] } else { (0..KINDS).collect() };
        for kind in kinds {
            let a = generate(&mut rng, r, c, kind);
            svd_case(&mut out, &a, true, true);
            if !fixture {
                svd_case(&mut out, &a, false, false);
            }
        }
    }

    // SVD solve.
    let solves = if fixture { 24 } else { 600 };
    for t in 0..solves {
        let r = 1 + rng.below(if t % 10 == 0 { 60 } else { 10 });
        let c = 1 + rng.below(10);
        let kind = rng.below(KINDS);
        let a = generate(&mut rng, r, c, kind);
        let bcols = if rng.unif() < 0.7 { 1 } else { 1 + rng.below(4) };
        let b = DMatrix::from_fn(r, bcols, |_, _| rng.normal());
        let svd = a.clone().svd(false, false);
        let smax = svd.singular_values[0];
        let eps = match rng.below(5) {
            0 => 0.0,
            1 => f64::EPSILON * r.max(c) as f64 * smax,
            2 => 1e-12,
            3 => smax * rng.unif(),
            _ => 1e300,
        };
        solve_case(&mut out, &a, &b, eps);
    }

    // Symmetric eigen: Gauss-Hermite Jacobi matrices and random symmetric matrices.
    let points: Vec<usize> = if fixture { vec![1, 2, 3, 5, 9, 15, 20, 31] } else { (1..=64).collect() };
    for p in points {
        eigen_case(&mut out, &jacobi(p));
    }
    let eigens = if fixture { 16 } else { 800 };
    for _ in 0..eigens {
        let n = 1 + rng.below(if fixture { 8 } else { 14 });
        let kind = rng.below(KINDS);
        let a = symmetric(&mut rng, n, kind);
        eigen_case(&mut out, &a);
    }
    if !fixture {
        for n in [20, 33, 48] {
            for kind in 0..KINDS {
                let a = symmetric(&mut rng, n, kind);
                eigen_case(&mut out, &a);
            }
        }
    }

    // Products.
    let dims: Vec<usize> = if fixture { vec![1, 5, 6, 9] } else { vec![1, 2, 5, 6, 7, 9, 13] };
    let ks: Vec<usize> = if fixture { vec![3, 6] } else { vec![1, 5, 6, 7, 31, 256, 257, 513] };
    for &m in &dims {
        for &n in &dims {
            for &k in &ks {
                let kind = rng.below(KINDS);
                let a = generate(&mut rng, m, k, kind);
                let b = { let kk = rng.below(KINDS); generate(&mut rng, k, n, kk) };
                mul_case(&mut out, "mul", &a, &b);
                if (!fixture && (k < 300 || m == n)) || (m == 6 && n == 9) {
                    let bt = b.transpose();
                    mul_case(&mut out, "mul_bt", &a, &bt);
                    let at = a.transpose();
                    mul_case(&mut out, "tr_mul", &at, &b);
                }
            }
            let a = { let kk = rng.below(KINDS); generate(&mut rng, m, n, kk) };
            let x = { let kk = rng.below(KINDS); generate(&mut rng, n, 1, kk) };
            mul_case(&mut out, "mul_vec", &a, &x);
            let y = { let kk = rng.below(KINDS); generate(&mut rng, m, 1, kk) };
            mul_case(&mut out, "tr_mul_vec", &a, &y);
        }
    }
    // Deep products cross matrixmultiply's 256-deep blocks; kind 8 underflows every product.
    let deep: &[(usize, usize, usize, usize)] = if fixture {
        &[(6, 257, 6, 1), (6, 256, 7, 8), (6, 300, 6, 3)]
    } else {
        &[(70, 300, 70, 1), (65, 40, 1030, 1), (200, 6, 6, 1), (6, 200, 6, 1), (7, 1000, 7, 1), (6, 256, 9, 8), (9, 600, 6, 8)]
    };
    for &(m, k, n, kind) in deep {
        let a = generate(&mut rng, m, k, kind);
        let b = generate(&mut rng, k, n, kind);
        mul_case(&mut out, "mul", &a, &b);
    }

    // Norms.
    let norms = if fixture { 12 } else { 400 };
    for _ in 0..norms {
        let r = 1 + rng.below(30);
        let c = 1 + rng.below(5);
        let a = { let kk = rng.below(KINDS); generate(&mut rng, r, c, kk) };
        norm_case(&mut out, &a);
    }

    // hypot as Rust computes it (the system libm on macOS).
    let hypots = if fixture { 40 } else { 20000 };
    for _ in 0..hypots {
        let scale = |rng: &mut Rng| 10f64.powi(rng.below(41) as i32 * 15 - 300);
        let x = rng.normal() * scale(&mut rng);
        let y = rng.normal() * scale(&mut rng);
        out.hypot.push(format!("[{},{},{}]", h(x), h(y), h(x.hypot(y))));
    }

    // Fused multiply-add as the NEON kernel performs it, stressed near underflow.
    let fmas = if fixture { 60 } else { 200000 };
    for _ in 0..fmas {
        let wide = rng.below(2) == 0;
        let e = |rng: &mut Rng| 2f64.powi(-(rng.below(1100) as i32));
        let w = |rng: &mut Rng| 2f64.powi(rng.below(2046) as i32 - 1022) * 2f64.powi(-(rng.below(53) as i32));
        let a = if wide { rng.normal() * w(&mut rng) } else { rng.normal() * e(&mut rng).sqrt() };
        let b = if wide { rng.normal() * w(&mut rng) } else { rng.normal() * e(&mut rng).sqrt() };
        let c = match rng.below(4) {
            0 => 0.0,
            1 => -0.0,
            2 => -(a * b) * (1.0 + rng.normal() * 1e-10),
            _ => rng.normal() * e(&mut rng),
        };
        out.fma.push(format!("[{},{},{},{}]", h(a), h(b), h(c), h(a.mul_add(b, c))));
    }

    println!("{{");
    println!("\"generator\":\"tests/reference/nalgebra (nalgebra 0.35.0, matrixmultiply 0.3.11, simba 0.10.2) mode={} seed={}\",", mode, seed);
    let sections = [
        ("svd", &out.svd),
        ("solve", &out.solve),
        ("eigen", &out.eigen),
        ("mul", &out.mul),
        ("norm", &out.norm),
        ("hypot", &out.hypot),
        ("fma", &out.fma),
    ];
    for (i, (name, items)) in sections.iter().enumerate() {
        println!("\"{}\":[", name);
        println!("{}", items.join(",\n"));
        println!("]{}", if i + 1 == sections.len() { "" } else { "," });
    }
    println!("}}");
}
