#include "../include/internal_angles.h"
#include <cmath>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  A.resizeLike(l_sqr);
  int row = l_sqr.rows();
  int col = l_sqr.cols();
  for (int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      double c2 = l_sqr(i,j);
      double b2 = l_sqr(i, (j+1)%col);
      double a2 = l_sqr(i, (j+2)%col);
      double cosc = (a2 + b2 - c2) / (2 * sqrt(a2) * sqrt(b2));
      A(i,j) = acos(cosc);
    }
  }
}
