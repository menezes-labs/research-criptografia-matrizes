use matrixcalculus1000::{generate_matrix, Matrix, MatrixKind, MAX_ORDER};
use rayon::ThreadPoolBuilder;
use std::env;
use std::process;
use std::time::Instant;

const SEED_STEP: u64 = 0x9E3779B97F4A7C15;

#[derive(Debug, Clone)]
struct Config {
    size: usize,
    kind: MatrixKind,
    seed: u64,
    repeat: usize,
    tolerance: Option<f64>,
    threads: Option<usize>,
    print_matrix: bool,
    csv: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            size: 5,
            kind: MatrixKind::Singular,
            seed: 42,
            repeat: 1,
            tolerance: None,
            threads: None,
            print_matrix: false,
            csv: false,
        }
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        eprintln!("use --help to see the available options");
        process::exit(2);
    }
}

fn run() -> Result<(), String> {
    let Some(config) = parse_config()? else {
        return Ok(());
    };

    if let Some(threads) = config.threads {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|error| format!("could not configure Rayon thread pool: {error}"))?;
    }

    if config.csv && config.print_matrix {
        return Err("--csv and --print cannot be used together".into());
    }

    if config.csv {
        println!(
            "iteration,size,kind,seed,rank,nullity,singular,tolerance,generation_ms,analysis_ms,row_swaps,det_sign,log_abs_det"
        );
    }

    let mut analysis_times_ms = Vec::with_capacity(config.repeat);

    for iteration in 0..config.repeat {
        let seed = config
            .seed
            .wrapping_add((iteration as u64).wrapping_mul(SEED_STEP));

        let generation_started = Instant::now();
        let matrix = generate_matrix(config.size, config.kind, seed)
            .map_err(|error| error.to_string())?;
        let generation_ms = generation_started.elapsed().as_secs_f64() * 1_000.0;

        if config.print_matrix && iteration == 0 {
            print_dense_matrix(&matrix);
        }

        let analysis_started = Instant::now();
        let analysis = matrix
            .analyze(config.tolerance)
            .map_err(|error| error.to_string())?;
        let analysis_ms = analysis_started.elapsed().as_secs_f64() * 1_000.0;
        analysis_times_ms.push(analysis_ms);

        verify_expected_kind(config.kind, analysis.singular)?;

        if config.csv {
            let log_abs_det = analysis
                .log_abs_determinant
                .map(|value| value.to_string())
                .unwrap_or_default();
            println!(
                "{},{},{},{},{},{},{},{:.17e},{:.6},{:.6},{},{},{}",
                iteration + 1,
                analysis.order,
                config.kind,
                seed,
                analysis.rank,
                analysis.nullity,
                analysis.singular,
                analysis.tolerance,
                generation_ms,
                analysis_ms,
                analysis.row_swaps,
                analysis.determinant_sign,
                log_abs_det
            );
        } else {
            println!("\n=== run {}/{} ===", iteration + 1, config.repeat);
            println!("matrix       : {}x{} ({})", config.size, config.size, config.kind);
            println!("seed         : {seed}");
            println!("rank         : {}", analysis.rank);
            println!("nullity      : {}", analysis.nullity);
            println!("singular     : {}", analysis.singular);
            println!("tolerance    : {:.6e}", analysis.tolerance);
            println!("row swaps    : {}", analysis.row_swaps);
            println!("generation   : {:.3} ms", generation_ms);
            println!("analysis     : {:.3} ms", analysis_ms);
            println!(
                "matrix memory : {:.2} MiB",
                matrix.memory_bytes() as f64 / 1_048_576.0
            );
            println!(
                "working memory: ~{:.2} MiB + pivot buffer",
                (matrix.memory_bytes() * 2) as f64 / 1_048_576.0
            );

            if let Some(log_abs_det) = analysis.log_abs_determinant {
                println!("det sign     : {}", analysis.determinant_sign);
                println!("ln|det(A)|   : {:.6}", log_abs_det);
            } else {
                println!("det(A)       : numerically zero");
            }
        }
    }

    if !config.csv && analysis_times_ms.len() > 1 {
        let min = analysis_times_ms.iter().copied().fold(f64::INFINITY, f64::min);
        let max = analysis_times_ms
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let avg = analysis_times_ms.iter().sum::<f64>() / analysis_times_ms.len() as f64;

        println!("\n=== benchmark summary ===");
        println!("runs : {}", analysis_times_ms.len());
        println!("min  : {:.3} ms", min);
        println!("avg  : {:.3} ms", avg);
        println!("max  : {:.3} ms", max);
    }

    Ok(())
}

fn verify_expected_kind(kind: MatrixKind, singular: bool) -> Result<(), String> {
    match (kind, singular) {
        (MatrixKind::Singular, false) => Err(
            "generated singular matrix was classified as invertible; review the tolerance".into(),
        ),
        (MatrixKind::Invertible, true) => Err(
            "strictly diagonally dominant matrix was classified as singular; review the tolerance"
                .into(),
        ),
        _ => Ok(()),
    }
}

fn print_dense_matrix(matrix: &Matrix) {
    if matrix.order() > 20 {
        eprintln!(
            "--print ignored: printing a {}x{} dense matrix would be excessive (limit: 20x20)",
            matrix.order(),
            matrix.order()
        );
        return;
    }

    println!("matrix:");
    for row in 0..matrix.order() {
        for col in 0..matrix.order() {
            print!("{:>11.5} ", matrix.get(row, col));
        }
        println!();
    }
}

fn parse_config() -> Result<Option<Config>, String> {
    let mut config = Config::default();
    let mut args = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                print_help();
                return Ok(None);
            }
            "--size" => {
                config.size = next_value(&mut args, "--size")?
                    .parse::<usize>()
                    .map_err(|_| "--size expects an integer".to_string())?;
            }
            "--kind" => {
                config.kind = next_value(&mut args, "--kind")?.parse::<MatrixKind>()?;
            }
            "--seed" => {
                config.seed = next_value(&mut args, "--seed")?
                    .parse::<u64>()
                    .map_err(|_| "--seed expects an unsigned 64-bit integer".to_string())?;
            }
            "--repeat" => {
                config.repeat = next_value(&mut args, "--repeat")?
                    .parse::<usize>()
                    .map_err(|_| "--repeat expects a positive integer".to_string())?;
            }
            "--tolerance" => {
                let value = next_value(&mut args, "--tolerance")?;
                config.tolerance = if value.eq_ignore_ascii_case("auto") {
                    None
                } else {
                    Some(
                        value
                            .parse::<f64>()
                            .map_err(|_| "--tolerance expects a non-negative number or 'auto'".to_string())?,
                    )
                };
            }
            "--threads" => {
                config.threads = Some(
                    next_value(&mut args, "--threads")?
                        .parse::<usize>()
                        .map_err(|_| "--threads expects a positive integer".to_string())?,
                );
            }
            "--print" => config.print_matrix = true,
            "--csv" => config.csv = true,
            unknown => return Err(format!("unknown option '{unknown}'")),
        }
    }

    if !(1..=MAX_ORDER).contains(&config.size) {
        return Err(format!("--size must be between 1 and {MAX_ORDER}"));
    }
    if config.repeat == 0 {
        return Err("--repeat must be at least 1".into());
    }
    if config.threads == Some(0) {
        return Err("--threads must be at least 1".into());
    }
    if let Some(tolerance) = config.tolerance {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err("--tolerance must be finite and non-negative".into());
        }
    }

    Ok(Some(config))
}

fn next_value(
    args: &mut impl Iterator<Item = String>,
    flag: &str,
) -> Result<String, String> {
    args.next()
        .ok_or_else(|| format!("missing value after {flag}"))
}

fn print_help() {
    println!(
        r#"matrixcalculus1000

Dense real-matrix singularity analysis up to 1000x1000.

USAGE:
    matrixcalculus1000 [OPTIONS]

OPTIONS:
    --size N              Matrix order, 1..=1000 [default: 5]
    --kind KIND           singular | invertible | random [default: singular]
    --seed N              Deterministic u64 seed [default: 42]
    --repeat N            Repeat generation + analysis N times [default: 1]
    --tolerance VALUE     Numerical pivot tolerance or 'auto' [default: auto]
    --threads N           Rayon worker threads [default: logical CPU count]
    --print               Print the matrix when N <= 20
    --csv                 Emit machine-readable benchmark rows
    -h, --help            Show this help

EXAMPLES:
    cargo run --release -- --size 1000 --kind singular
    cargo run --release -- --size 1000 --kind invertible --repeat 5
    cargo run --release -- --size 256 --kind random --seed 2026 --threads 8
    cargo run --release -- --size 500 --kind singular --repeat 10 --csv
"#
    );
}
