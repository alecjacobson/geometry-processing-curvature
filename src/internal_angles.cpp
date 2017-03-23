#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
	A.resize(l_sqr.rows(), 3);

	//Use law-of-cosines to calculate the internal angles given 
	//the edge lengths of the face.
	//https://www.mathsisfun.com/algebra/trig-solving-sss-triangles.html

	for (int i = 0; i < l_sqr.rows();i++) {
		double a = l_sqr(i, 0);
		double b = l_sqr(i, 1);
		double c = l_sqr(i, 2);

		double _A = acos((b + c - a) / (2 * sqrt(b*c)));
		double _B = acos((c + a - b) / (2 * sqrt(c*a)));
		double _C = acos((a + b - c) / (2 * sqrt(a*b)));

		A(i, 0) = _A; A(i, 1) = _B; A(i, 2) = _C;
	}
}
