#include "../include/principal_curvatures.h"
#include <igl/pinv.h>
#include <igl/slice.h>
#include <igl/adjacency_matrix.h>
#include <set>
#include <Eigen/Eigen>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  // Replace with your code
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);

  // Grab adjacent vertex
  Eigen::SparseMatrix<int> adj_matrix;
  igl::adjacency_matrix(F, adj_matrix);

  for (int i = 0; i < V.rows(); i ++) {
    // Find adjacent vertices
    std::set<int> adj_ver;
    // Idea inspired by Eigen Sparse Matrix documentation
    for (Eigen::SparseMatrix<int>::InnerIterator it1(adj_matrix, i); it1; ++it1) {
      int level_1 = it1.row();
      adj_ver.insert(level_1);
      for (Eigen::SparseMatrix<int>::InnerIterator it2(adj_matrix, i); it2; ++it2) {
        int level_2 = it2.row();
        if (level_2 != i ) {
          adj_ver.insert(level_2);
        }
      }
    }
    
    // Find P: (vi - v)
    Eigen::MatrixXd P(adj_ver.size(), 3);
    int j = 0;
    for(auto cur_vi : adj_ver) {
      P.row(j) = V.row(cur_vi) - V.row(i);
      j ++;
    }

    // PCA: by issue #18 and documentation
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(P.transpose()*P);
    // Since R3, two most:
    Eigen::VectorXd u = es.eigenvectors().col(1);
    Eigen::VectorXd v = es.eigenvectors().col(2);
    // one leasT:
    Eigen::VectorXd w = P * es.eigenvectors().col(0);

    // Construct coefficients
    Eigen::MatrixXd coeff(P.rows(), 5);
    coeff.col(0) = P * u;
    coeff.col(1) = P * v;
    coeff.col(2) = (P * u).array().square();
    coeff.col(3) = (P * u).cwiseProduct(P * v);
    coeff.col(4) = (P * v).array().square();

    // Compute As
    Eigen::MatrixXd inverse_coeff;
    igl::pinv(coeff, inverse_coeff);
    Eigen::VectorXd A = inverse_coeff * w;

    // Compute variables
    double E = 1 + A(0) * A(0);
    double F = A(0) * A(1);
    double G = 1 + A(1) * A(1);
    double common = std::sqrt(E + G - 1);
    double e = (2 * A(2)) / common;
    double f = A(3) / common;
    double g = (2 * A(4)) / common;

    // Construct S
    Eigen::Matrix2d left, Right, S;
    left << e, f,
            f, g;
    Right << E, F,
             F, G;
    S = - left * Right.inverse();

    // PCA on S
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_S(S);
    K1(i) = es_S.eigenvalues()(0);
    K2(i) = es_S.eigenvalues()(1);
    D1.row(i) = es_S.eigenvectors()(0, 0) * u + es_S.eigenvectors()(0, 1) * v;
    D2.row(i) = es_S.eigenvectors()(1, 0) * u + es_S.eigenvectors()(0, 1) * v;

  }
}
