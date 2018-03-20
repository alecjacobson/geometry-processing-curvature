#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
	int F_count = l_sqr.rows();
	A.resize(F_count, l_sqr.cols());
	for (int i = 0; i < F_count; i++) {
		double a2 = l_sqr(i, 0);
		double b2 = l_sqr(i, 1);
		double c2 = l_sqr(i, 2);

		double a = sqrt(a2);
		double b = sqrt(b2);
		double c = sqrt(c2);

		double cosa = (b2 + c2 - a2) / (2.0 * b * c);
		double cosb = (a2 + c2 - b2) / (2.0 * a * c);
		double cosc = (a2 + b2 - c2) / (2.0 * a * b);

		A(i, 0) = acos(cosa);
		A(i, 1) = acos(cosb);
		A(i, 2) = acos(cosc);
	}
}
