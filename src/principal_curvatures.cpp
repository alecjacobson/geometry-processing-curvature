#include "../include/principal_curvatures.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <igl/adjacency_matrix.h>
#include <set>
#include <igl/pinv.h>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2) {
  // Replace with your code
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(), 3);
  D2 = Eigen::MatrixXd::Zero(V.rows(), 3);

  Eigen::SparseMatrix<int> adj;
  igl::adjacency_matrix(F, adj);

  for (int i = 0; i < V.rows(); i++) {
    // Gather two ring vertices
    std::set<int> twoRing;
    for(Eigen::SparseMatrix<int>::InnerIterator it(adj, i); it; ++it) {
      int oneRing = it.row();
      twoRing.insert(oneRing);
      for(Eigen::SparseMatrix<int>::InnerIterator it1(adj, oneRing); it1; ++it1) {
        int twoRingi = it1.row();
        twoRing.insert(twoRingi);
      }
    }
    // Construct P
    Eigen::MatrixXd P(twoRing.size(), 3);
    int j = 0;
    for (auto idx : twoRing) {
      P.row(j++) = V.row(idx) - V.row(i);
    }

    // Compute plane passing through V, eigen decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(P.transpose() * P);
    Eigen::VectorXd u = eigenSolver.eigenvectors().col(2);
    Eigen::VectorXd v = eigenSolver.eigenvectors().col(1);
    Eigen::VectorXd w = eigenSolver.eigenvectors().col(0);
    Eigen::VectorXd B = P * w;
    Eigen::MatrixXd S(twoRing.size(), 2);
    S.col(0) = P * u;
    S.col(1) = P * v;

    // Solve for known coefficients
    Eigen::MatrixXd ax(P.rows(), 5);
    ax << S, S.col(0).cwiseProduct(S.col(0)), S.col(0).cwiseProduct(S.col(1)), S.col(1).cwiseProduct(S.col(1));
    Eigen::MatrixXd X;
    igl::pinv(ax, X);
    Eigen::VectorXd A = X * B;

    // Shape operator
    double e, f, g, E, Fd, G;
    E = 1 + A(0) * A(0);
    Fd = A(0) * A(1);
    G = 1 + A(1) * A(1);
    e = (2 * A[2]) / sqrt(A[0] * A[0] + 1 + A[1] * A[1]);
    f = (1 * A[3]) / sqrt(A[0] * A[0] + 1 + A[1] * A[1]);
    g = (2 * A[4]) / sqrt(A[0] * A[0] + 1 + A[1] * A[1]);
    Eigen::MatrixXd left(2, 2);
    left << e, f,
            f, g;
    Eigen::MatrixXd right(2, 2);
    right << E, Fd,
             Fd, G;
    Eigen::MatrixXd shapeOp(2, 2);
    shapeOp = - left * right.inverse();

    // Eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolverSO(shapeOp);
    K1(i) = eigenSolverSO.eigenvalues()(0);
    K2(i) = eigenSolverSO.eigenvalues()(1);
    D1.row(i) = eigenSolverSO.eigenvectors()(1, 1) * u + eigenSolverSO.eigenvectors()(1, 0) * v;
    D2.row(i) = eigenSolverSO.eigenvectors()(0, 1) * u + eigenSolverSO.eigenvectors()(0, 0) * v;
  }

}