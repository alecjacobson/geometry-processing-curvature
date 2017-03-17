#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <iostream>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  int num_v = V.rows();
  Eigen::SparseMatrix<double> L(num_v, num_v);
  Eigen::SparseMatrix<double> M(num_v, num_v);
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  
  Eigen::SparseMatrix<double> M_inv;
  igl::invert_diag(M, M_inv);

  Eigen::MatrixXd Hn;
  Hn = M_inv * L * V;
  
  Eigen::MatrixXd N(num_v, 3);
  igl::per_vertex_normals(V, F, N);
  
  H.resize(num_v);
  double d;
  for (int i = 0; i < num_v; i++) {
    d = ( abs(Hn(i,0))/N(i,0) + abs(Hn(i,1))/N(i,1) + abs(Hn(i,2))/N(i,2) ) / 3;
    H[i] = 0.5 * d;
  }
}
