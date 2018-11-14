#include "../include/internal_angles.h"
#include <cmath>
#include <math.h>
#include <iostream>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // use law of cosines...
  // A = arccos((-a^2 + b^2 +c^2)/2bc)
  A.resize(l_sqr.rows(), l_sqr.cols());
  double a2, b2, c2, a, b, c, B, C;
  for (int i = 0; i < l_sqr.rows(); i++) {
    a2 = l_sqr(i, 0);
    b2 = l_sqr(i, 1);
    c2 = l_sqr(i, 2);
    a = std::sqrt(a2);
    b = std::sqrt(b2);
    c = std::sqrt(c2);
    B = std::acos((-b2 + a2 + c2) / (2 * a * c));
    C = std::acos((-c2 + b2 + a2) / (2 * a * b));
    A(i, 0) = M_PI - B - C; // saves us one call to acos lol
    A(i, 1) = B;
    A(i, 2) = C;
  }
}
