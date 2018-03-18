#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A.resize(l_sqr.rows(), l_sqr.cols());

  for(int i = 0; i < A.rows(); i++) {
    for(int j = 0; j < 3; j++) {
      double a2 = l_sqr(i, (j + 1) % 3);
      double b2 = l_sqr(i, (j + 2) % 3);
      double c2 = l_sqr(i, (j + 0) % 3);
      double a = sqrt(a2);
      double b = sqrt(b2);

      A(i, j) = acos((a2 + b2 - c2) / (2.0 * a * b));
    }
  }
}
