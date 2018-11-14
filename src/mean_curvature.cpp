#include "../include/mean_curvature.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/invert_diag.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // H =~ (M.inv() * L * V) . N
  Eigen::SparseMatrix<double> M, M_inv, L;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, M_inv);
  
  Eigen::MatrixXd H_mat = M_inv * L * V;
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);
  H = Eigen::VectorXd::Zero(V.rows());
  double dot;
  for (int i = 0; i < H_mat.rows(); i++) {
    dot = H_mat.row(i).dot(N.row(i));
    if (dot < 0) {
      H(i) = H_mat.row(i).norm();
    } else {
      H(i) = (-1.0) * H_mat.row(i).norm();
    }
  }
}
