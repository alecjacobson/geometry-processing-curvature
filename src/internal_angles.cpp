#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A.resize(l_sqr.rows(), 3);
  for (int i = 0; i < l_sqr.rows(); i++) {
    double asq = l_sqr(i, 0);
    double bsq = l_sqr(i, 1);
    double csq = l_sqr(i, 2);
    A(i, 0) = acos((bsq + csq - asq) / (2*sqrt(bsq*csq)));
    A(i, 1) = acos((asq + csq - bsq) / (2*sqrt(asq*csq)));
    A(i, 2) = acos((asq + bsq - csq) / (2*sqrt(asq*bsq)));
  }
}
