#include "../include/internal_angles.h"
#include <math.h> 

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
	A = Eigen::MatrixXd::Zero(l_sqr.rows(), l_sqr.cols());
  	for (int i = 0; i < l_sqr.rows(); i++) {
  		for (int j = 0; j < 3; j++) {
  			double oppo = l_sqr(i, j);
  			double side1 = l_sqr(i, (j+1) % 3);
  			double side2 = l_sqr(i, (j+2) % 3);
  			A(i, j) = std::acos((side1 * side1 + side2 * side2 - oppo * oppo) / (2 * side1 * side2));
  		} 
  	}
}
