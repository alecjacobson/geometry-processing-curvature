#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  int n = l_sqr.rows();
  A.resize(n, 3);

  for(int i = 0; i < n; i++){
    for (int j = 0; j < 3; j++){
      // opposite side to angle
      double a = l_sqr(i, (j+1)% 3);
      double b = l_sqr(i, (j+2)% 3);
      double c = l_sqr(i, j);
      A(i,j) = acos((a + b - c) / (2 * sqrt(a) * sqrt(b)));
    }
  }
}