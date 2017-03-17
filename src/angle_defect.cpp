#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>

void angle_defect(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::VectorXd & D)
{
	Eigen::MatrixXd l_sqr;
	igl::squared_edge_lengths(V, F, l_sqr);
	Eigen::MatrixXd angles;
	internal_angles(l_sqr, angles);
	D = Eigen::VectorXd::Constant(V.rows(), 2 * M_PI);

	for (int f = 0; f < F.rows(); f++) {
		D(F(f, 0)) -= angles(f, 0);
		D(F(f, 1)) -= angles(f, 1);
		D(F(f, 2)) -= angles(f, 2);
	}
}
