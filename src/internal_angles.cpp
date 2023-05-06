#include "../include/internal_angles.h"

#include <Eigen/Dense>

void internal_angles(const Eigen::MatrixXd & l_sqr, Eigen::MatrixXd & A) {
	// Add with your code

	int fnum = l_sqr.rows();
	A.resize(fnum, 3);
	for (int i = 0; i < fnum; i++) {
		// c^2 = a^2 + b^2 - 2abcos
		// length: 12, 23, 31
		// angle order, c, a, b
		double f12 = l_sqr(i, 0);
		double f23 = l_sqr(i, 1);
		double f31 = l_sqr(i, 2);

		double a = (f12 * f12 + f31 * f31 - f23 * f23) / 2 / f12 / f31;
		A(i, 0) = std::acos(a);

		double b = (f12 * f12 + f23 * f23 - f31 * f31) / 2 / f12 / f23;
		A(i, 1) = std::acos(b);

		double c = (f23 * f23 + f31 * f31 - f12 * f12) / 2 / f23 / f31;
		A(i, 2) = std::acos(c);
	}

	return;
}
