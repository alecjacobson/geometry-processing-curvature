#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
   	const Eigen::MatrixXd & l_sqr,
  	Eigen::MatrixXd & A)
{
	A.resizeLike(l_sqr);
  	for (int i = 0; i < l_sqr.rows(); i++) {
  		for (int j = 0; j < 3; j++) {
	  		double a = l_sqr(i, j % 3);
	  		double b = l_sqr(i, (j + 1) % 3);
	  		double c = l_sqr(i, (j + 2) % 3);

	  		// A in units radians
	  		A(i, j) = acos((b + c - a) / (2.0 * sqrt(b * c)));
  		}
  	}
}
