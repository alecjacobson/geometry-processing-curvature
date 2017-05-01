#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
	auto computeAngle = [](auto a, auto b, auto c) { return acos((c - a - b ) / (2*sqrt(a)*sqrt(b)));};
	
	A.resize(l_sqr.rows(), 3);

	for (int i = 0; i < l_sqr.rows(); ++i) 
	{
		auto v = l_sqr.row(i);
		double a = v[0], b = v[1], c = v[2];
		A(i, 0) = computeAngle(b, c, a);
		A(i, 1) = computeAngle(a, c, b);
		A(i, 2) = computeAngle(a, b, c);
	}
}
