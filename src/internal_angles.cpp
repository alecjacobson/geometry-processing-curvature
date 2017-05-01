#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
	int f = l_sqr.rows();
	A.resize(f, 3);

	for (int i = 0; i < f; ++i)
	{
		for (int v = 0; v < 3; ++v)
		{
			auto c2 = l_sqr(i, v);
			auto b2 = l_sqr(i, (v + 1) % 3);
			auto a2 = l_sqr(i, (v + 2) % 3);

			A(i, v) = std::acos((a2 + b2 - c2) / (2 * std::sqrt(a2*b2)));
		}
	}
}
