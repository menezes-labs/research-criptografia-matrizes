use crate::{Matrix, MatrixError};
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixKind {
    Singular,
    Invertible,
    Random,
}

impl fmt::Display for MatrixKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Singular => write!(f, "singular"),
            Self::Invertible => write!(f, "invertible"),
            Self::Random => write!(f, "random"),
        }
    }
}

impl FromStr for MatrixKind {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "singular" => Ok(Self::Singular),
            "invertible" => Ok(Self::Invertible),
            "random" => Ok(Self::Random),
            other => Err(format!(
                "invalid matrix kind '{other}'; expected singular, invertible, or random"
            )),
        }
    }
}

pub fn generate_matrix(
    order: usize,
    kind: MatrixKind,
    seed: u64,
) -> Result<Matrix, MatrixError> {
    match kind {
        MatrixKind::Singular => generate_singular(order, seed),
        MatrixKind::Invertible => generate_invertible(order, seed),
        MatrixKind::Random => generate_random(order, seed),
    }
}

fn generate_singular(order: usize, seed: u64) -> Result<Matrix, MatrixError> {
    let mut rng = SplitMix64::new(seed);
    let mut matrix = Matrix::zeros(order)?;

    for row in 0..order {
        for col in 0..order {
            matrix.set(row, col, rng.signed_unit());
        }
    }

    if order == 1 {
        matrix.set(0, 0, 0.0);
    } else if order == 2 {
        for col in 0..order {
            matrix.set(1, col, 2.0 * matrix.get(0, col));
        }
    } else {
        let last = order - 1;
        for col in 0..order {
            let dependent_value = 1.25 * matrix.get(0, col) - 0.75 * matrix.get(1, col);
            matrix.set(last, col, dependent_value);
        }
    }

    Ok(matrix)
}

fn generate_invertible(order: usize, seed: u64) -> Result<Matrix, MatrixError> {
    let mut rng = SplitMix64::new(seed);
    let mut matrix = Matrix::zeros(order)?;

    for row in 0..order {
        let mut off_diagonal_sum = 0.0;

        for col in 0..order {
            if row == col {
                continue;
            }

            let value = 0.25 * rng.signed_unit();
            matrix.set(row, col, value);
            off_diagonal_sum += value.abs();
        }

        // Strict row diagonal dominance guarantees nonsingularity.
        let margin = 1.0 + rng.unit_interval();
        matrix.set(row, row, off_diagonal_sum + margin);
    }

    Ok(matrix)
}

fn generate_random(order: usize, seed: u64) -> Result<Matrix, MatrixError> {
    let mut rng = SplitMix64::new(seed);
    let mut matrix = Matrix::zeros(order)?;

    for row in 0..order {
        for col in 0..order {
            matrix.set(row, col, rng.signed_unit());
        }
    }

    Ok(matrix)
}

#[derive(Debug, Clone, Copy)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    fn unit_interval(&mut self) -> f64 {
        const DENOMINATOR: f64 = (1u64 << 53) as f64;
        ((self.next_u64() >> 11) as f64) / DENOMINATOR
    }

    fn signed_unit(&mut self) -> f64 {
        2.0 * self.unit_interval() - 1.0
    }
}
