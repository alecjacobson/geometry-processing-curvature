#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    A.resize(l_sqr.rows(), l_sqr.cols());
    for (int i = 0; i < l_sqr.rows(); i++) {
        A(i, 0) = acos((l_sqr(i, 1) + l_sqr(i, 2) - l_sqr(i, 0)) / (2.0 * sqrt(l_sqr(i, 1) * l_sqr(i, 2))));
        A(i, 1) = acos((l_sqr(i, 2) + l_sqr(i, 0) - l_sqr(i, 1)) / (2.0 * sqrt(l_sqr(i, 2) * l_sqr(i, 0))));
        A(i, 2) = acos((l_sqr(i, 0) + l_sqr(i, 1) - l_sqr(i, 2)) / (2.0 * sqrt(l_sqr(i, 0) * l_sqr(i, 1))));
    }
}
