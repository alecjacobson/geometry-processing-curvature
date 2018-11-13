#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include <igl/invert_diag.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include "igl/per_vertex_normals.h"
#include <iostream>
#include <math.h>

using namespace std;

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  H = Eigen::VectorXd::Zero(V.rows());
  Eigen::SparseMatrix<double> M(V.rows(), V.cols());
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SparseMatrix<double> L(V.rows(), V.cols());
  igl::cotmatrix(V, F, L);

  // Calculate M^(-1) * L * V
  Eigen::SparseMatrix<double> Minv;
  igl::invert_diag(M, Minv);
  Eigen::MatrixXd Q = L*V;
  Eigen::MatrixXd Hn = Minv * Q;
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);
  for (int i = 0; i < Hn.rows(); i++) {
        // magnitude
        double mag = sqrt(Hn(i, 0)*Hn(i, 0) + Hn(i, 1)*Hn(i, 1) + Hn(i, 2)*Hn(i, 2));

        // dot product
        double dir = Hn(i, 0)*N(i, 0) + Hn(i, 1)*N(i, 1) + Hn(i, 2)*N(i, 2);
        if (dir >= 0) {
            H(i) = -mag;
        } else {
            H(i) = mag;
        }
    }
}
