row_vector col_sums(matrix X) {
  row_vector[cols(X)] s;
  for (j in 1:cols(X)) s[j] = sum(col(X, j));
  return s;
}

vector row_sums(matrix X) {
  vector[rows(X)] s;
  for (i in 1:rows(X)) s[i] = sum(X[i]);
  return s;
}


matrix log_mix_par(matrix theta, matrix A, matrix B) {
  int ncol = cols(A);
  int nrow = rows(A);
  matrix[ncol, nrow] M;
  for (j in 1:nrow) {
    for (i in 1:ncol) {
      M[i, j] = log_mix(theta[i, j], A[i, j], B[i, j]);
    }
  }
  return M;
}
