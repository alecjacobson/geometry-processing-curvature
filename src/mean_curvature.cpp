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
  // Replace with your code
  H = Eigen::VectorXd::Zero(V.rows());

  // Compute L
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // Compute M
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  // Invert M (M is diagonal)
  Eigen::SparseMatrix<double> M_inv;
  igl::invert_diag(M, M_inv);

  // Compute normals
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  // Compute H = M^-1 L V
  Eigen::MatrixXd H_vec = M_inv*L*V;

  for (int i = 0; i < V.rows(); ++i) {
    H(i) = H_vec.row(i).norm();

    // Check sign
    if (H_vec.row(i).dot(N.row(i)) > 0)
      H(i) *= -1.0;

  } // end loop i

}
