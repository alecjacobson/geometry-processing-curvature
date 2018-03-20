#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
	A.resize(l_sqr.rows(),3);

	for (int i = 0; i < l_sqr.rows(); i++)
	{
		auto lengths = l_sqr.row(i);
		for (int l = 0; l < 3; l++)
		{
			float a_sq = lengths(l % 3); // 1,2
			float b_sq = lengths((l + 1) % 3); // 2,0
			float c_sq = lengths((l + 2) % 3); // 0,1
			double angle = acos((-a_sq + b_sq + c_sq) / (2 * sqrt(b_sq) * sqrt(c_sq)));
			A(i, l) = angle;
		}
	}
}
