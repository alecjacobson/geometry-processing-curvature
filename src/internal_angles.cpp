#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A = Eigen::MatrixXd(l_sqr.rows(), 3);
  for (int m = 0; m < A.rows(); m++)
    for (int i = 0; i < 3; i++)
    {
      int j = (i + 1) % 3, k = (i + 2) % 3;
      double cosine = (l_sqr(m, i) - l_sqr(m, j) - l_sqr(m, k)) /
        (-2 * sqrt(l_sqr(m, j)*l_sqr(m, k)));
      A(m, i) = acos(cosine);
    }
}
