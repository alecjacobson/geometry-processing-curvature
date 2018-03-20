#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>
#include <iostream>
#include <cmath>
#include <functional>
double sign_func(double x)
{
    if (x > 0)
        return +1.0;
    else if (x == 0)
        return 0.0;
    else
        return -1.0;
}

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  H = Eigen::VectorXd::Zero(V.rows());
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  Eigen::SparseMatrix<double> invertM;
  igl::invert_diag(M, invertM);
  auto mean_normal  = invertM * L * V;

  Eigen::MatrixXd N;
  H.resize(V.rows(), 1);
  igl::per_vertex_normals(V, F, N);
  auto mag = (mean_normal.array() * mean_normal.array()).rowwise().sum().sqrt();
  auto sign =  (mean_normal.array() * N.array()).rowwise().sum().unaryExpr(std::ptr_fun(sign_func));
  H = (mag.array() * sign.array()).matrix();
  
  // for (int i = 0; i < mean_normal.rows(); i++){
  // 	auto row = mean_normal.row(i);
  // 	auto dotprod = N.row(i).dot(row);
  // 	int sign = -1;
  // 	if (dotprod > 0) sign = 1; 
  // 	H(i) = sign * row.norm();
  // }
}
