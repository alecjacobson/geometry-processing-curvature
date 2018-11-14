#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  A = Eigen::MatrixXd::Zero(l_sqr.rows(), l_sqr.cols());
  for (int m = 0; m < A.rows(); m++) {
    for (int i = 0; i < 3; i++) {
      double l1 = l_sqr(m, i);
      double l2 = l_sqr(m, (i + 1) % 3);
      double l3 = l_sqr(m, (i + 2) % 3);
      A(m, i) = acos((l3 + l2 - l1) / (2.0 * sqrt(l3 * l2)));
    }
  }
}
