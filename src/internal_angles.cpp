#include "../include/internal_angles.h"
#include <cmath>

void internal_angles(
    const Eigen::MatrixXd &l_sqr,
    Eigen::MatrixXd &A)
{
  A.resize(l_sqr.rows(), 3);
  // Add with your code
  for (int i = 0; i < l_sqr.rows(); i++)
  {
    //get edge squares
    double a = l_sqr(i, 0); //1,2
    double b = l_sqr(i, 1); //2,0
    double c = l_sqr(i, 2); //0,1

    //by cosine law
    A(i, 0) = acos((b + c - a) / (2.00 * sqrt(b) * sqrt(c)));
    A(i, 1) = acos((a + c - b) / (2.00 * sqrt(a) * sqrt(c)));
    A(i, 2) = acos((a + b - c) / (2.00 * sqrt(a) * sqrt(b)));
  }
}
