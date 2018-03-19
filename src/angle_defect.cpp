#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <igl/PI.h>

using namespace std;

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	Eigen::MatrixXd edgeLengthSquare;
	igl::squared_edge_lengths(V, F, edgeLengthSquare);

	Eigen::MatrixXd internalAngles;
	internal_angles(edgeLengthSquare, internalAngles);

	D = Eigen::VectorXd::Zero(V.rows());
	D.array() += 2 * igl::PI;
	for(int ii = 0; ii < F.rows(); ii++){
		D(F(ii,0)) -= internalAngles(ii,0);
		D(F(ii,1)) -= internalAngles(ii,1);
		D(F(ii,2)) -= internalAngles(ii,2);
  }
}
