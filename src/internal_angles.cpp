#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A.resize(l_sqr.rows(), 3);

  // Apply the cosine law to compute each internal angle
  for (int i = 0; i < l_sqr.rows(); i++) {
	A(i, 0) = (l_sqr(i, 0) - l_sqr(i,1) - l_sqr(i,2)) / (-2.0 * std::sqrt(l_sqr(i,1) * l_sqr(i,2))); 
	A(i, 1) = (l_sqr(i, 1) - l_sqr(i,2) - l_sqr(i,0)) / (-2.0 * std::sqrt(l_sqr(i,0) * l_sqr(i,2))); 
	A(i, 2) = (l_sqr(i, 2) - l_sqr(i,0) - l_sqr(i,1)) / (-2.0 * std::sqrt(l_sqr(i,0) * l_sqr(i,1)));
	A(i, 0) = std::acos(A(i,0));
 	A(i, 1) = std::acos(A(i,1));
	A(i, 2) = std::acos(A(i,2));
  }
}
