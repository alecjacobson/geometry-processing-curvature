#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
	int n = l_sqr.rows();
	A = Eigen::MatrixXd(n, 3);
	for (int i = 0; i < n; ++i) 
	{
		double a = l_sqr(i, 0);
		double b = l_sqr(i, 1);
		double c = l_sqr(i, 2);

		double cos_A = (a - b - c) / (-2 * sqrt(b * c));
		double cos_B = (b - c - a) / (-2 * sqrt(c * a));
		double cos_C = (c - a - b) / (-2 * sqrt(a * b));

		double A_angle = acos(cos_A);
		double B_angle = acos(cos_B);
		double C_angle = acos(cos_C);

		A.row(i) = Eigen::RowVector3d(A_angle, B_angle, C_angle);
	}
}
