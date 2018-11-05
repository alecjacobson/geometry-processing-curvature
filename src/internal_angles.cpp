#include "../include/internal_angles.h"
#include <cmath>

// returns theta_a in radians
double cosinerule(double la_sq, double lb_sq, double lc_sq); 


void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  A.resize(l_sqr.rows(), l_sqr.cols());

  // iterate
  for (int i = 0; i < l_sqr.rows(); i++) {
    A(i, 0) = cosinerule(l_sqr(i, 0), l_sqr(i, 1), l_sqr(i, 2));
    A(i, 1) = cosinerule(l_sqr(i, 1), l_sqr(i, 0), l_sqr(i, 2));
    A(i, 2) = cosinerule(l_sqr(i, 2), l_sqr(i, 1), l_sqr(i, 0));
  }

}

// returns the angle in radians formed by vertex a
double cosinerule(double la_sq, double lb_sq, double lc_sq){
  double cosinea = (lc_sq + lb_sq - la_sq)/(2.0*sqrt(lc_sq*lb_sq));
  double anglea = acos(cosinea);
  return anglea;
}
