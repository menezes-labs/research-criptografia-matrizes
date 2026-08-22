use matrixcalculus1000::{generate_matrix, Matrix, MatrixError, MatrixKind, MAX_ORDER};

#[test]
fn detects_known_singular_matrix() {
    let matrix = Matrix::from_vec(
        3,
        vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            5.0, 7.0, 9.0,
        ],
    )
    .unwrap();

    let analysis = matrix.analyze(None).unwrap();
    assert_eq!(analysis.rank, 2);
    assert_eq!(analysis.nullity, 1);
    assert!(analysis.singular);
    assert_eq!(analysis.determinant_sign, 0);
    assert!(analysis.log_abs_determinant.is_none());
}

#[test]
fn detects_identity_as_invertible() {
    let matrix = Matrix::identity(32).unwrap();
    let analysis = matrix.analyze(None).unwrap();

    assert_eq!(analysis.rank, 32);
    assert_eq!(analysis.nullity, 0);
    assert!(!analysis.singular);
    assert_eq!(analysis.determinant_sign, 1);
    assert_eq!(analysis.log_abs_determinant, Some(0.0));
}

#[test]
fn generated_singular_matrix_has_deficient_rank() {
    let matrix = generate_matrix(64, MatrixKind::Singular, 12345).unwrap();
    let analysis = matrix.analyze(None).unwrap();

    assert!(analysis.singular);
    assert!(analysis.rank < 64);
    assert!(analysis.nullity >= 1);
}

#[test]
fn generated_diagonally_dominant_matrix_is_invertible() {
    let matrix = generate_matrix(64, MatrixKind::Invertible, 12345).unwrap();
    let analysis = matrix.analyze(None).unwrap();

    assert!(!analysis.singular);
    assert_eq!(analysis.rank, 64);
    assert_eq!(analysis.nullity, 0);
}

#[test]
fn supports_structural_allocation_at_order_1000() {
    let matrix = Matrix::identity(MAX_ORDER).unwrap();

    assert_eq!(matrix.order(), 1000);
    assert_eq!(matrix.len(), 1_000_000);
    assert_eq!(matrix.memory_bytes(), 8_000_000);
    assert_eq!(matrix.get(999, 999), 1.0);
}

#[test]
fn rejects_order_above_limit_before_allocation() {
    let error = Matrix::zeros(MAX_ORDER + 1).unwrap_err();
    assert!(matches!(
        error,
        MatrixError::InvalidOrder { order } if order == MAX_ORDER + 1
    ));
}

#[test]
fn rejects_non_finite_input() {
    let error = Matrix::from_vec(2, vec![1.0, 0.0, f64::NAN, 1.0]).unwrap_err();
    assert!(matches!(error, MatrixError::NonFiniteValue { index: 2 }));
}

#[test]
fn rejects_negative_tolerance() {
    let matrix = Matrix::identity(2).unwrap();
    let error = matrix.analyze(Some(-1.0)).unwrap_err();

    assert!(matches!(error, MatrixError::InvalidTolerance { .. }));
}
