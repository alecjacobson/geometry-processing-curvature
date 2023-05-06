#include "../include/internal_angles.h"

// Returns internal angle opposite A
static double inline cosine_law(double a_sqr, double b_sqr, double c_sqr) {
  return acos((b_sqr + c_sqr - a_sqr) / (2*sqrt(b_sqr)*sqrt(c_sqr)));
}

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A.resize(l_sqr.rows(), l_sqr.cols());
  for (int i = 0; i < A.rows(); i++) {
    A(i,0) = cosine_law(l_sqr(i,0),l_sqr(i,1),l_sqr(i,2));
    A(i,1) = cosine_law(l_sqr(i,1),l_sqr(i,0),l_sqr(i,2));
    A(i,2) = cosine_law(l_sqr(i,2),l_sqr(i,0),l_sqr(i,1));
  }
}
