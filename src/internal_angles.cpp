#include "../include/internal_angles.h"
#include <cmath>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  A.resizeLike(l_sqr);

    for(int i = 0; i < l_sqr.rows(); ++i) {
        auto&& el = l_sqr.row(i);
        for(int j = 0; j < 3; ++j) {
            double a = el((j+0)%3);
            double b = el((j+1)%3);
            double c = el((j+2)%3);

            double ct = -(a - b - c) / (2 * std::sqrt(b*c));
            A(i,j) = std::acos(ct);




        }
    }
}
