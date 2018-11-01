#include "../include/mean_curvature.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  Eigen::MatrixXd N;
  Eigen::SparseMatrix<double> L, M;

  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::per_vertex_normals(V, F, N);

  Eigen::MatrixXd Hn = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>(M).solve(L * V);
  H = -(Hn * N.transpose()).diagonal();
}
