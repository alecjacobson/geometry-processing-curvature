#include "../include/internal_angles.h"
#include <iostream>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code

	// For any triangle, one of its angles can be expressed in terms of its lengths via cosines:
	// a^2 = b^2 + c^2 - 2*b*c*cos(theta)
	// from which we can find theta.

	// We're going to be cycling through vertices in the order 0, 1, 2.

	A.resize(l_sqr.rows(), 3);

	// Loop through each set of squared lengths
	for (int ii = 0; ii < l_sqr.rows(); ii++)
	{
		// Extract this set of edge lengths
		Eigen::VectorXd this_l_sqr = l_sqr.row(ii);
		double a, b, c;

		// Compute the angle corresponding to the first vertex
		a = this_l_sqr(0);
		b = this_l_sqr(1);
		c = this_l_sqr(2);

		A(ii, 0) = std::acos((-a + b + c)/std::sqrt(b*c)/2.0);

		// Compute the angle corresponding to the second vertex
		a = this_l_sqr(1);
		b = this_l_sqr(2);
		c = this_l_sqr(0);

		A(ii, 1) = std::acos((-a + b + c)/std::sqrt(b*c)/2.0);

		// Compute the angle corresponding to the third vertex
		a = this_l_sqr(2);
		b = this_l_sqr(0);
		c = this_l_sqr(1);

		A(ii, 2) = std::acos((-a + b + c)/std::sqrt(b*c)/2.0);

	}

	return;

}
