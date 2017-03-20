#include "../include/angle_defect.h"
#include <igl/squared_edge_lengths.h>
#include "../include/internal_angles.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	int n_F = F.rows();
	int n_V = V.rows();
	
	Eigen::MatrixXd l_sqr;
	igl::squared_edge_lengths(V, F, l_sqr);

	Eigen::MatrixXd IA;
	internal_angles(l_sqr, IA);

	D = Eigen::VectorXd(n_V);
	D.setConstant(2 * M_PI);
	for(int i = 0; i < n_F; ++i)
	{
		// IA and F are in same order
		for (int k = 0; k < 3; ++k) 
		{
			D(F(i, k)) -= IA(i, k);
		}
	}
}
