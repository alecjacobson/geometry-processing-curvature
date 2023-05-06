#include "../include/mean_curvature.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
	// Calculate mass matrix M and cotangent matrix L
	Eigen::SparseMatrix<double> M, L;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::cotmatrix(V, F, L);
	
	// Caluclate M^-1
	Eigen::SparseMatrix<double> M_inv;
	igl::invert_diag(M, M_inv);
	
	// Calculate BH
	Eigen::MatrixXd BH(V.rows(), V.cols());
	BH = M_inv * L* V;
	
	// Calculate per vertex normal to get sign of mean curvature
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);
	
	// Calculate H
	H.resize(V.rows());
	H = (-BH * N.transpose()).diagonal();
}
