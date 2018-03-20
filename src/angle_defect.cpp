#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <math.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	// Calculate internal angles
	Eigen::MatrixXd l_sqr;
	igl::squared_edge_lengths(V, F, l_sqr);
	Eigen::MatrixXd A;
	internal_angles(l_sqr, A);

	// Calculate angle defects
	D.resize(V.rows());
	D.setConstant(2.0 * M_PI);
	for (int i = 0; i < F.rows(); i++) {
		D(F(i, 0)) -= A(i, 0);
		D(F(i, 1)) -= A(i, 1);
		D(F(i, 2)) -= A(i, 2);
	}
}
