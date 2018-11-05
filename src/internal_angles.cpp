#include "../include/internal_angles.h"

/**
* @brief Solves with the cosine law to produce the cosine of the angle
* opposite to the edge c.
*/
double cosineLaw(double a, double b, double c) {
	double top = a * a + b * b - c * c;
	double bottom = 2 * a * b;
	return top / bottom;
}

void internal_angles(const Eigen::MatrixXd & l_sqr, Eigen::MatrixXd & A) {

	A.resizeLike(l_sqr);
	int faceCount = l_sqr.rows();
	for (int i = 0; i < faceCount; i++) {
		double a = sqrt(l_sqr(i, 0));
		double b = sqrt(l_sqr(i, 1));
		double c = sqrt(l_sqr(i, 2));

		double cosA = cosineLaw(b, c, a);
		double cosB = cosineLaw(a, c, b);
		double cosC = cosineLaw(a, b, c);

		double angleA = acos(cosA);
		double angleB = acos(cosB);
		double angleC = acos(cosC);

		A(i, 0) = angleA;
		A(i, 1) = angleB;
		A(i, 2) = angleC;
	}
}
