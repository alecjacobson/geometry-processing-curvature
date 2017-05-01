#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

using namespace Eigen;

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
	SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	for (int i = 0; i < M.outerSize(); ++i)
		for (SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
			it.valueRef() = 1. / it.value();

	MatrixXd H_vec(M*(-L)*V);

	MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	H.resize(H_vec.rows());

	for (int i = 0; i < H_vec.rows(); ++i)
		H(i) = H_vec.row(i).dot(N.row(i));
}
