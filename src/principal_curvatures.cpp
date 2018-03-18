#include "../include/principal_curvatures.h"

#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/pinv.h>
#include <iostream>
#include <set>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  K1.resize(V.rows());
  K2.resize(V.rows());
  D1.resize(V.rows(), 3);
  D2.resize(V.rows(), 3);

  Eigen::MatrixXd N;
  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F, A);
  igl::per_vertex_normals(V, F, N);

  // There's probably a more elegant way to do this whole thing but I'm tired
  for(int i = 0; i < V.rows(); i++) {
    // Get two-ring of v
    std::set<int> two_ring;
    for(Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it) {
      int one_ring_i = it.row();
      two_ring.insert(one_ring_i);
      for(Eigen::SparseMatrix<double>::InnerIterator it2(A,one_ring_i); it2; ++it2) {
        int two_ring_i = it2.row();
        two_ring.insert(two_ring_i);
      }
    }

    int j = 0;
    Eigen::MatrixXd P(two_ring.size(), 3);
    for(auto index : two_ring) {
      P.row(j++) = V.row(index) - V.row(i);
    }

    // Compute tangent plane
    Eigen::JacobiSVD<Eigen::MatrixXd> PTPsvd(P.transpose() * P, Eigen::ComputeFullU);
    Eigen::MatrixXd R = PTPsvd.matrixU();
    Eigen::MatrixXd PR = P * R;
    Eigen::MatrixXd S = PR.leftCols(2);
    Eigen::VectorXd B = PR.rightCols(1);

    // Fit quadratic surface
    Eigen::MatrixXd quadratic(S.rows(), 5);
    quadratic << S, S.col(0).array().square(), S.col(0).cwiseProduct(S.col(1)), S.col(1).array().square();

    Eigen::MatrixXd X;
    igl::pinv(quadratic, X);
    Eigen::VectorXd A = X * B; // Coeffs a1..a5

    // Construct shape operator and solve for curvature
    double e = 2.0 * A[2] / sqrt(A[0] * A[0] + 1.0 + A[1] * A[1]);
    double f = A[3] / sqrt(A[0] * A[0] + 1.0 + A[1] * A[1]);
    double g = 2.0 * A[4] / sqrt(A[0] * A[0] + 1.0 + A[1] * A[1]);
    double E = 1.0 + A[0] * A[0];
    double F = A[0] * A[1];
    double G = 1.0 + A[1] * A[1];
    Eigen::MatrixXd S_right_inv(2, 2);
    S_right_inv << E, F,
                   F, G;
    Eigen::MatrixXd S_left(2, 2);
    S_left << e, f,
              f, g;

    Eigen::MatrixXd shape_op = -S_right_inv.transpose().lu().solve(S_left.transpose()).transpose();
    Eigen::EigenSolver<Eigen::MatrixXd> shape_op_eigen(shape_op);
    Eigen::MatrixXd principal_directions = (R.leftCols(2) * shape_op_eigen.eigenvectors().real()).transpose();

    D1.row(i) = principal_directions.row(0);
    D2.row(i) = principal_directions.row(1);
    K1[i] = shape_op_eigen.eigenvalues()[0].real();
    K2[i] = shape_op_eigen.eigenvalues()[1].real();
  }
}
