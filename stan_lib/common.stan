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
