#include "../include/mean_curvature.h"
#include <Eigen/Sparse>
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
  Eigen::SparseMatrix<double> L, M, Mi;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, Mi);

  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  Eigen::MatrixXd MLV = Mi * L * V;

  H.resize(MLV.rows());
  for (int i = 0; i < H.size(); i++) {
    int sign = (MLV.row(i).dot(N.row(i)) >= 0) ? 1 : -1;
    H(i) = sign * MLV.row(i).norm();
  }
}
