#include "../include/internal_angles.h"
#include <math.h>

// Return the angle corresponding to side C
double get_angle(double side_a, double side_b, double side_c) {
  double num = side_a + side_b - side_c;
  double denom = 2 * sqrt(side_a) * sqrt(side_b);
  return acos(num / denom);
}

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A.resize(l_sqr.rows(), 3);
  double side_a, side_b, side_c;
  for (int i = 0; i < l_sqr.rows(); i++) {
    side_a = l_sqr(i, 0);
    side_b = l_sqr(i, 1);
    side_c = l_sqr(i, 2);
    A(i, 0) = get_angle(side_b, side_c, side_a);
    A(i, 1) = get_angle(side_c, side_a, side_b);
    A(i, 2) = get_angle(side_a, side_b, side_c);
  }
}
