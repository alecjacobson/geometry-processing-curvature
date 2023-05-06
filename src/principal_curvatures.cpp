#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <iostream>
#include <math.h>
#include <Eigen/Eigenvalues> 
#include <set>

using namespace std;
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

  // compute adj matrix for computing 2-rings
  Eigen::SparseMatrix<int> A;
  igl::adjacency_matrix(F,A);

  for (int ii = 0; ii < V.rows(); ii++){
    // compute matrix P: k-by-3 matrix
    Eigen::MatrixXd P;
    {
      std::set<int> twoRings; // store two ring indices
      for (Eigen::SparseMatrix<int>::InnerIterator iterOneRing(A,ii); iterOneRing; ++iterOneRing){
        int oneRingVIdx = iterOneRing.row();
        for (Eigen::SparseMatrix<int>::InnerIterator iterTwoRing(A,oneRingVIdx); iterTwoRing; ++iterTwoRing){
          twoRings.insert(iterTwoRing.row());
        }
      }
      twoRings.erase(ii); // erase itself
      P.resize(twoRings.size(), 3);
      int idx = 0;
      for (auto it = twoRings.begin(); it != twoRings.end(); ++it){
        P.row(idx) = V.row(*it) - V.row(ii);
        idx++;
      }
    }

    // Eigen decomposition for P^T*P
    Eigen::VectorXd u, v, w;
    {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(P.transpose() * P);
      u = eig.eigenvectors().col(2); // first PC
      v = eig.eigenvectors().col(1); // second PC
      w = eig.eigenvectors().col(0); // height
    }

    // compute quadratic surface coefficients
    Eigen::VectorXd A;
    {
      Eigen::VectorXd B = P * w; // project P to w coordinates
      
      Eigen::MatrixXd M; // M = [u, v, u^2, uv, v^2]
      M.resize(P.rows(), 5);
      M.col(0) = P * u;
      M.col(1) = P * v;
      M.col(2) = M.col(0).cwiseProduct(M.col(0));
      M.col(3) = M.col(0).cwiseProduct(M.col(1));
      M.col(4) = M.col(1).cwiseProduct(M.col(1));

      Eigen::MatrixXd Minv;
      igl::pinv(M, Minv);

      A = Minv * B; // least square solve
    }

    // compute shape operator
    Eigen::Matrix2d S;
    {
      double E = 1 + A(0)*A(0);
      double F = A(0) * A(1);
      double G = 1 + A(1)*A(1);
      double e = 2*A(2) / sqrt(A(0)*A(0) + 1 + A(1)*A(1));
      double f = A(3) / sqrt(A(0)*A(0) + 1 + A(1)*A(1));
      double g = 2*A(4) / sqrt(A(0)*A(0) + 1 + A(1)*A(1));

      Eigen::Matrix2d leftMat, rightMat;
      leftMat << e, f, f, g;
      rightMat << E, F, F, G;
      S = -leftMat * rightMat.inverse();
    }

    // eigen decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(S);
    double k1 = eig.eigenvalues()(1);
    double k2 = eig.eigenvalues()(0);

    // lift to eucliden space
    Eigen::Vector3d PC1Dir = eig.eigenvectors()(1,1) * u + eig.eigenvectors()(1,0) * v;
    Eigen::Vector3d PC2Dir = eig.eigenvectors()(0,1) * u + eig.eigenvectors()(0,0) * v;

    // Store to output matrix
    K1(ii) = k1;
    K2(ii) = k2;
    D1.row(ii) = PC1Dir;
    D2.row(ii) = PC2Dir;
  }

}
