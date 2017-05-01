#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	//Set the angle defect initially to 2PI
	D.resize(V.rows());
	D.setConstant(2 * M_PI);

	//Calculate the internal angles per vertex/corner
	Eigen::MatrixXd sq_length, int_angle;
	igl::squared_edge_lengths(V, F, sq_length);
	internal_angles(sq_length, int_angle);

	//Subtract the internal angle for every incident face on the ith vertex
	for (int f = 0; f < F.rows(); f++) {
		for (int i = 0; i < 3; i++) {
			D(F(f, i)) -= int_angle(f, i);
		}
	}

}
