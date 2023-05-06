#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseCholesky>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  	const Eigen::MatrixXd & V,
  	const Eigen::MatrixXi & F,
  	Eigen::VectorXd & H)
{
	// compute mass matrix and cotangent matrix
	Eigen::SparseMatrix<double> M, M_I, L;
	igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
	igl::cotmatrix(V, F, L);

	// compute mean curvature normal
	Eigen::MatrixXd HN;
	// Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> cholesky(M);
	// HN = cholesky.solve(L * V);
	igl::invert_diag(M, M_I);
	HN = M_I * L * V;

	// compute normal
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	// add sign to curvature
	H = Eigen::VectorXd::Zero(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		if (N.row(i).dot(HN.row(i)) > 0) {
			H(i) = HN.row(i).norm();
		} else if (N.row(i).dot(HN.row(i)) < 0) {
			H(i) = (-1.0) * HN.row(i).norm();
		}
	}
}
