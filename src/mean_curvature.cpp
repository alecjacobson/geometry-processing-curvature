#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

using namespace Eigen;

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{

  H = Eigen::VectorXd::Zero(V.rows());

  // cotangent matrix
  SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

  // mass matrix
  SparseMatrix<double> M, M_inv;
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, M_inv);

  // normals: to determine the sign
  MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  // curvature normals
  MatrixXd HN = M_inv * L * V;

  // check and update H
  for (int i = 0; i < V.rows(); i++) {
    double cur = HN.row(i).dot(N.row(i));
    int norm = HN.row(i).norm();
    if (cur > 0) {
      H(i) = -norm;
    } 
    else {
      H(i) = norm;
    }
  }


}
