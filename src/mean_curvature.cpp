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
  Eigen::SparseMatrix<double> L, M, M_1;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT,M);
  igl::invert_diag(M,M_1);
  // curvature normal
  Eigen::MatrixXd N_hat = M_1 * L * V;
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V,F,N);
  for (int i=0; i<V.rows(); i++) {
    Eigen::RowVectorXd n = N.row(i);
    Eigen::RowVectorXd n_hat = N_hat.row(i);
    double nn = n.dot(n_hat);
    if (nn > 0) {
      H[i] = N_hat.row(i).norm();
    }
    else {
      H[i] = -N_hat.row(i).norm();
    }
  }
}
