#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <Eigen/Sparse>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  
  Eigen::MatrixXd Hn = M.cwiseInverse() * L * V;
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);
  H.resize(V.rows());
  
  for (int i = 0; i < Hn.rows(); i++) {
    H[i] = Hn.row(i).norm();
    // make sure sign is correct
    if (Hn.row(i).dot(N.row(i)) < 0) {
      H.row(i) *= -1;
    }
  }
}
