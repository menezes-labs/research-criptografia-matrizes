pub mod generator;
pub mod matrix;

pub use generator::{generate_matrix, MatrixKind};
pub use matrix::{Analysis, Matrix, MatrixError, MAX_ORDER};
