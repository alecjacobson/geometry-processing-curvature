#include "../include/mean_curvature.h"
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <iostream>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  H = Eigen::VectorXd::Zero(V.rows());
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Y;
  Eigen::MatrixXd N;

  igl::cotmatrix(V, F, L); 
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, Y);
  igl::per_vertex_normals(V, F, N);
  std::cout << Y.cols() << std::endl;
  std::cout << L.rows() << std::endl;
  std::cout << V.rows() << std::endl;

  Eigen::MatrixXd temp = Y * L * V * N.transpose();
  for (int i = 0; i < V.rows(); i++)
    H(i) = -temp(i,i);
}
