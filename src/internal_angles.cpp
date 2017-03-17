#include "../include/internal_angles.h"
#include <cmath>

void internal_angles(
	const Eigen::MatrixXd & l_sqr,
	Eigen::MatrixXd & A)
{
	// Add with your code
	auto num_faces = l_sqr.rows();

	const Eigen::MatrixXd l = l_sqr.cwiseSqrt();

	A.resize(num_faces, 3);
	for (int f = 0; f < num_faces; f++) {
		for (int v = 0; v < 3; v++) {
			double numerator = l_sqr(f, (v + 1) % 3) + l_sqr(f, (v + 2) % 3) - l_sqr(f, v);
			double denominator = 2 * l(f, (v + 1) % 3)*l(f, (v + 2) % 3);
			A(f, v) = acos(numerator / denominator);
		}
	}

}
