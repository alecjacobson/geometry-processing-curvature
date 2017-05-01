#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::VectorXd & H)
{
	
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	
	Eigen::SparseMatrix<double> M, M_inv;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, M_inv);

	
	auto h_vecs = (M_inv*L*V).eval();
	
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	H.resize(V.rows());
	for (int i = 0; i < h_vecs.rows(); i++) {
		H(i) = (h_vecs.row(i).dot(N.row(i)) < 0) ? h_vecs.row(i).norm() : -h_vecs.row(i).norm();
	}
}
