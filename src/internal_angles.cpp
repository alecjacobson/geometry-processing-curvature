#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    A.resize(l_sqr.rows(), 3);
  // Add with your code
  for (int i = 0; i < l_sqr.rows(); i++) {
      for (int j = 0; j < 3; j++) {
          double a = l_sqr(i, (j+1)%3);
          double b = l_sqr(i, (j+2)%3);
          double c = l_sqr(i, j);
          // Using law of cosines to calculate internal angle
          A(i, j) = acos( (a+b-c) / (2.*sqrt(a)*sqrt(b)) );
      }
  }
}
