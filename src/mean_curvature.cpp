#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

using namespace Eigen;
using namespace igl;

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
	SparseMatrix<double> M, L, Minv;
	massmatrix(V, F, MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
	cotmatrix(V, F, L);
	invert_diag(M, Minv);
	
	MatrixXd Hnorm = Minv * L * V;

	H = Hnorm.rowwise().norm(); 

	Eigen::MatrixXd VN;
	per_vertex_normals(V, F, VN);
	
	for (int i = 0; i < V.rows(); ++i)
	{
		if (Hnorm.row(i).dot(VN.row(i)) > 0) // smth wrong with sign, but seems to work
		{
			H.row(i) = -1 * H.row(i);
		}
	}
}
