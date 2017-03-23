#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include "igl/squared_edge_lengths.h"
#include "igl/adjacency_matrix.h"


using namespace Eigen;
using namespace igl;

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	D = VectorXd::Zero(V.rows());
	
	MatrixXd lsqd, A;
	squared_edge_lengths(V, F, lsqd);
	internal_angles(lsqd, A);

	//iterate over every vertex and add up angles
	
	for (int i = 0; i < F.rows(); ++i) 
	{
		for (int k = 0; k < 3; ++k) 
		{
			auto f = F(i, k);
			auto l = A(i, k);

			D[f] += l;
		}
	}

	VectorXd GG = 2 * M_PI *  VectorXd::Ones(V.rows());

	D = GG - D;

}	
