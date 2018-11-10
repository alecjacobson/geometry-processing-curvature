#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	double PI = 3.14159265358979;
	Eigen::MatrixXd l_sqr;
	igl::squared_edge_lengths(V, F, l_sqr);
	Eigen::MatrixXd A;
	internal_angles(l_sqr, A);
	D = Eigen::VectorXd::Constant(V.rows(), 2 * PI);
	for (int i = 0; i < A.rows(); i++) {
		for (int j = 0; j < 3; j++) {
			D(F(i, j)) -= A(i, j);
		}
	}
}
