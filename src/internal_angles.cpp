#include "../include/internal_angles.h"
#include <math.h> 

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  int num = l_sqr.rows();
  A.resize(num, 3);

  for (int i = 0; i < num; i++) {
    for (int j = 0; j < 3; j++) {
      // side: b, c   angle to: a
      double a = l_sqr(i, j % 3);
      double b = l_sqr(i, (j + 1) % 3);
      double c = l_sqr(i, (j + 2) % 3);
      
      double cosAngle = ((b + c - a) * 1.0 / (2 * sqrt(b * c)));
      A(i, j) = acos(cosAngle);
    }
  }

}
