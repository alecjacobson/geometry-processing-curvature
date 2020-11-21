#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  A.resize(l_sqr.rows(), l_sqr.cols());

  for (int i = 0; i < l_sqr.rows(); ++i) {
    const Eigen::RowVector3d& l = l_sqr.row(i);   
    for (int j = 0; j < 3; ++j) {
      int a = j;
      int b = (j + 1)%3;
      int c = (j + 2)%3;

      A(i, a) = std::acos(
          (l(a) - l(b) - l(c)) / (-2.0*std::sqrt(l(b)*l(c)))
        );
    } // end loop j
  } // end loop i
}
