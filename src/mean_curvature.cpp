#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  H.resize(V.rows());
  // Compute cot laplacian and mass matrix
  Eigen::SparseMatrix<double> L, M, inverse_M;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, inverse_M);

  // HN = M-1 LV
  Eigen::MatrixXd HN = inverse_M * L * V;

  // Compute normal vectors
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  // Compute H and its sign
  for(int i = 0; i < HN.rows(); i ++) {
    H(i) = HN.row(i).norm();
    if (HN.row(i).dot( N.row(i) ) < 0) {
      // Non-consistency
      H(i) *= -1;
    }
  }
}
