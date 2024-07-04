void InverseMatrix(double mu[3][3]) {
  double tmp[3][3];
  double delta;
  
  delta = mu[0][0] * mu[1][1] * mu[2][2] + mu[1][0] * mu[2][1] * mu[0][2]
    + mu[2][0] * mu[0][1] * mu[1][2] - mu[0][0] * mu[1][2] * mu[2][1]
    - mu[0][1] * mu[1][0] * mu[2][2] - mu[0][2] * mu[1][1] * mu[2][0];
  
  tmp[0][0] = (mu[1][1] * mu[2][2] - mu[1][2] * mu[2][1]) / delta;
  tmp[1][0] = -(mu[1][0] * mu[2][2] - mu[1][2] * mu[2][0]) / delta;
  tmp[2][0] = (mu[1][0] * mu[2][1] - mu[1][1] * mu[2][0]) / delta;
  
  tmp[0][1] = -(mu[0][1] * mu[2][2] - mu[0][2] * mu[2][1]) / delta;
  tmp[1][1] = (mu[0][0] * mu[2][2] - mu[0][2] * mu[2][0]) / delta;
  tmp[2][1] = -(mu[0][0] * mu[2][1] - mu[0][1] * mu[2][0]) / delta;
  
  tmp[0][2] = (mu[0][1] * mu[1][2] - mu[0][2] * mu[1][1]) / delta;
  tmp[1][2] = -(mu[0][0] * mu[1][2] - mu[0][2] * mu[1][0]) / delta;
  tmp[2][2] = (mu[0][0] * mu[1][1] - mu[0][1] * mu[1][0]) / delta;
  
  for (long long int i = 0; i < 3; i++) {
    for (long long int j = 0; j < 3; j++) {
      mu[i][j] = tmp[i][j];
    }
  }
  return;
}
