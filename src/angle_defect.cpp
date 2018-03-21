#include "../include/angle_defect.h"
#include <igl/squared_edge_lengths.h>
#include "../include/internal_angles.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  // D = Eigen::VectorXd::Zero(V.rows());

	// First, we need the internal angles, for which we need the squared lengths via igl.
	Eigen::MatrixXd l_sqr;
	igl::squared_edge_lengths(V, F, l_sqr);

	Eigen::MatrixXd A;
	internal_angles(l_sqr, A);

	// Initialize the angle defects to an original value of 2*PI
	D.resize(V.rows());
	for (int ii = 0; ii < V.rows(); ii++)
		D(ii) = 2.0*M_PI;

	// Now the angle defect is just 2*pi minus the sum of internal angles of each triangle incident on a vertex
	for (int ii = 0; ii < A.rows(); ii++)
	{
		// Loop through the faces incident on each vertex
		for (int jj = 0; jj < 3; jj++)
			D(F(ii, jj)) -= A(ii, jj);

	}

  	return;

}
