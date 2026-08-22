use rayon::prelude::*;
use std::error::Error;
use std::fmt;

pub const MAX_ORDER: usize = 1000;
const PARALLEL_THRESHOLD: usize = 192;
const AUTO_TOLERANCE_FACTOR: f64 = 64.0;

#[derive(Debug, Clone, PartialEq)]
pub enum MatrixError {
    InvalidOrder { order: usize },
    InvalidDataLength { expected: usize, actual: usize },
    NonFiniteValue { index: usize },
    InvalidTolerance { value: f64 },
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOrder { order } => write!(
                f,
                "matrix order must be between 1 and {MAX_ORDER}; received {order}"
            ),
            Self::InvalidDataLength { expected, actual } => write!(
                f,
                "invalid data length: expected {expected} elements, received {actual}"
            ),
            Self::NonFiniteValue { index } => {
                write!(f, "matrix contains NaN or infinity at flat index {index}")
            }
            Self::InvalidTolerance { value } => {
                write!(f, "tolerance must be finite and non-negative; received {value}")
            }
        }
    }
}

impl Error for MatrixError {}

#[derive(Debug, Clone)]
pub struct Analysis {
    pub order: usize,
    pub rank: usize,
    pub nullity: usize,
    pub singular: bool,
    pub tolerance: f64,
    pub row_swaps: usize,
    pub determinant_sign: i8,
    pub log_abs_determinant: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct Matrix {
    order: usize,
    data: Vec<f64>,
}

impl Matrix {
    pub fn from_vec(order: usize, data: Vec<f64>) -> Result<Self, MatrixError> {
        validate_order(order)?;
        let expected = order * order;

        if data.len() != expected {
            return Err(MatrixError::InvalidDataLength {
                expected,
                actual: data.len(),
            });
        }

        if let Some((index, _)) = data.iter().enumerate().find(|(_, value)| !value.is_finite()) {
            return Err(MatrixError::NonFiniteValue { index });
        }

        Ok(Self { order, data })
    }

    pub fn zeros(order: usize) -> Result<Self, MatrixError> {
        validate_order(order)?;
        Ok(Self {
            order,
            data: vec![0.0; order * order],
        })
    }

    pub fn identity(order: usize) -> Result<Self, MatrixError> {
        let mut matrix = Self::zeros(order)?;
        for i in 0..order {
            matrix.set(i, i, 1.0);
        }
        Ok(matrix)
    }

    #[inline]
    pub fn order(&self) -> usize {
        self.order
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.order + col]
    }

    #[inline]
    pub(crate) fn set(&mut self, row: usize, col: usize, value: f64) {
        let index = row * self.order + col;
        self.data[index] = value;
    }

    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    pub fn memory_bytes(&self) -> usize {
        self.data.len() * std::mem::size_of::<f64>()
    }

    pub fn infinity_norm(&self) -> f64 {
        self.data
            .chunks(self.order)
            .map(|row| row.iter().map(|value| value.abs()).sum::<f64>())
            .fold(0.0, f64::max)
    }

    pub fn auto_tolerance(&self) -> f64 {
        let scale = self.infinity_norm().max(1.0);
        f64::EPSILON * self.order as f64 * scale * AUTO_TOLERANCE_FACTOR
    }

    pub fn is_singular(&self, tolerance: Option<f64>) -> Result<bool, MatrixError> {
        Ok(self.analyze(tolerance)?.singular)
    }

    /// Computes numerical rank through Gaussian elimination with partial pivoting.
    ///
    /// The matrix is copied before elimination, so the original data remains unchanged.
    /// For sufficiently large matrices, elimination of rows below each pivot is parallelized
    /// with Rayon. Complexity is O(n^3) time and O(n^2) additional memory.
    pub fn analyze(&self, tolerance: Option<f64>) -> Result<Analysis, MatrixError> {
        let tolerance = match tolerance {
            Some(value) if !value.is_finite() || value < 0.0 => {
                return Err(MatrixError::InvalidTolerance { value });
            }
            Some(value) => value,
            None => self.auto_tolerance(),
        };

        let n = self.order;
        let mut work = self.clone();
        let mut pivot_row = 0usize;
        let mut row_swaps = 0usize;
        let mut determinant_sign = 1i8;
        let mut log_abs_determinant = 0.0f64;

        for pivot_col in 0..n {
            if pivot_row == n {
                break;
            }

            let mut best_row = pivot_row;
            let mut best_abs = work.get(best_row, pivot_col).abs();

            for row in (pivot_row + 1)..n {
                let candidate = work.get(row, pivot_col).abs();
                if candidate > best_abs {
                    best_abs = candidate;
                    best_row = row;
                }
            }

            if best_abs <= tolerance {
                continue;
            }

            if best_row != pivot_row {
                work.swap_rows(pivot_row, best_row);
                row_swaps += 1;
                determinant_sign = -determinant_sign;
            }

            let pivot = work.get(pivot_row, pivot_col);
            if pivot < 0.0 {
                determinant_sign = -determinant_sign;
            }
            log_abs_determinant += pivot.abs().ln();

            let first_row_below = pivot_row + 1;
            if first_row_below < n {
                let pivot_snapshot = work.data[pivot_row * n..(pivot_row + 1) * n].to_vec();
                let tail = &mut work.data[first_row_below * n..];

                let eliminate = |target: &mut [f64]| {
                    let value = target[pivot_col];
                    if value.abs() <= tolerance {
                        target[pivot_col] = 0.0;
                        return;
                    }

                    let factor = value / pivot;
                    target[pivot_col] = 0.0;

                    for col in (pivot_col + 1)..n {
                        target[col] -= factor * pivot_snapshot[col];
                    }
                };

                let rows_below = n - first_row_below;
                if n >= PARALLEL_THRESHOLD && rows_below >= 8 {
                    tail.par_chunks_mut(n).for_each(eliminate);
                } else {
                    tail.chunks_mut(n).for_each(eliminate);
                }
            }

            pivot_row += 1;
        }

        let rank = pivot_row;
        let singular = rank < n;

        Ok(Analysis {
            order: n,
            rank,
            nullity: n - rank,
            singular,
            tolerance,
            row_swaps,
            determinant_sign: if singular { 0 } else { determinant_sign },
            log_abs_determinant: if singular {
                None
            } else {
                Some(log_abs_determinant)
            },
        })
    }

    fn swap_rows(&mut self, row_a: usize, row_b: usize) {
        if row_a == row_b {
            return;
        }

        for col in 0..self.order {
            let a = row_a * self.order + col;
            let b = row_b * self.order + col;
            self.data.swap(a, b);
        }
    }
}

fn validate_order(order: usize) -> Result<(), MatrixError> {
    if !(1..=MAX_ORDER).contains(&order) {
        return Err(MatrixError::InvalidOrder { order });
    }
    Ok(())
}
