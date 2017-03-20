#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
	int n_V = V.rows();
	
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);

	Eigen::SparseMatrix<double> M_inv;
	igl::invert_diag(M, M_inv);

	Eigen::MatrixXd A = M_inv*L*V;

	H = A.rowwise().norm(); // get magnitude

	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	// get sign
	for(int i = 0; i < n_V; ++i)
	{
		int sign = A.row(i).dot(N.row(i)) > 0 ? 1 : -1;
		H.row(i) = sign*H.row(i);
	}

}
