#include "../include/internal_angles.h"

#include <cmath>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    A.resizeLike(l_sqr);

    // return angle opposite of edge `c`
    auto angle = [](double a2, double b2, double c2) {
        return std::acos((a2 + b2 - c2) / (2 * sqrt(a2) * sqrt(b2)));
    };

    double a2, b2, c2;
    for (int i = 0; i < l_sqr.rows(); ++i) 
    {
        // row = [A, B, C] where `A` opposite of edge `a`

        a2 = l_sqr(i, 0);
        b2 = l_sqr(i, 1);
        c2 = l_sqr(i, 2);

        A(i, 0) = angle(b2, c2, a2);
        A(i, 1) = angle(c2, a2, b2);
        A(i, 2) = angle(a2, b2, c2);
    }
}