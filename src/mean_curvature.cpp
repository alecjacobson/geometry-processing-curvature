#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"
#include "igl/massmatrix.h"
#include "igl/squared_edge_lengths.h"
#include "igl/per_vertex_normals.h"

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Construct Laplacian and mass matrices
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> M_inv;
   
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, M);
  igl::cotmatrix(V, F, L);
  igl::invert_diag(M, M_inv);

  // Compute mean curvature vectors as 1/2 (M^(-1) L V)
  Eigen::MatrixXd H_vec(V.rows(), 3);
  for (int i = 0; i <= 2; i++) {
	H_vec.col(i) = 0.5 * M_inv * L * V.col(i); 
  }

  // Produce mean curvature by comparing normals
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);
 
  H = Eigen::VectorXd::Zero(V.rows());
  for (int i = 0; i < V.rows(); i++) {
	if (N(i, 2) > 0) {
		H(i) = H_vec.row(i).norm();
	}
	else if (N(i, 2) < 0) {
		H(i) = -H_vec.row(i).norm();
	}
	else {
		H(i) = 0;
	}
  }
}
