#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A.resize(l_sqr.rows(), 3);
  // Law of Cosines:
  for (int i = 0; i < l_sqr.rows(); i ++) {
    double a2 = l_sqr(i, 0);
    double b2 = l_sqr(i, 1);
    double c2 = l_sqr(i, 2);

    double cos_C = (a2 + b2 - c2) / (2.0 * sqrt(a2) * sqrt(b2));
    double angle_C = std::acos(cos_C);
    
    double cos_A = (b2 + c2 - a2) / (2.0 * sqrt(b2) * sqrt(c2));
    double angle_A = std::acos(cos_A);

    double angle_B = 180 - angle_C - angle_A;

    A(i, 0) = angle_A;
    A(i, 1) = angle_B;
    A(i, 2) = angle_C;
  }
}
