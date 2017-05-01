#include "../include/angle_defect.h"
#include <igl/massmatrix.h>
#include <igl/squared_edge_lengths.h>
#include "internal_angles.h"

using namespace Eigen;

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	VectorXd area = M.diagonal();

	MatrixXd l_sqr;
	igl::squared_edge_lengths(V, F, l_sqr);

	MatrixXd angle;
	internal_angles(l_sqr, angle);

	int f = F.rows(), n = V.rows();
	D = 2*M_PI*VectorXd::Ones(n);

	for (int i = 0; i < f; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			auto v = F(i, j);
			D(v) -= angle(i, j);
		}
	}

	D.cwiseQuotient(area);
}
