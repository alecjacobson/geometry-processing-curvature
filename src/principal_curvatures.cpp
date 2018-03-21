#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <vector>
#include <iostream>
#include <Eigen/Eigenvalues>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  // Replace with your code
  Eigen::SparseMatrix<int> A;
  igl::adjacency_matrix(F, A);

  int num_V = F.maxCoeff();
  for(int v_idx = 0; v_idx < num_V; v_idx++){
    // get the two-ring neighbors
    Eigen::MatrixXd P;
    two_ring_P_matrix(V, A, v_idx, P);
    // Eigen decomposition on (P.T P)
    Eigen::EigenSolver<Eigen::MatrixXd> solver_P(P.transpose() * P);

    /// DEBUG
    // Eigen::MatrixXd A(2,3);
    // std::cout << solver_P.eigenvectors().col(0)(0) << std::endl;
    // Eigen::MatrixXd ev = solver_P.eigenvectors().cast<double>();
    // std::cout << ev << std::endl;
    // //A(0,0) = solver_P.eigenvectors().col(0)(0);
    // //std::cout << A << std::endl;
    // break;
    //A.transpose().ldlt().solve(P.transpose());
  }

  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);
}

void two_ring_P_matrix(
  const Eigen::MatrixXd & V,
  Eigen::SparseMatrix<int>& A,
  const int v_idx,
  Eigen::MatrixXd & P)
{
  // get the two-ring neighbors
  std::vector<int> one_ring;
  std::vector<int> two_ring;
  Eigen::SparseMatrix<int> o = A.row(v_idx);
  for (int k=0; k < o.outerSize(); ++k)
  {
      for (Eigen::SparseMatrix<int>::InnerIterator it(o,k); it; ++it)
      {
          one_ring.push_back(it.col());// col index (here it is equal to k)
          two_ring.push_back(it.col());
      }
  }
  for(std::vector<int>::iterator o_itr = one_ring.begin() ; o_itr != one_ring.end(); ++o_itr){
    Eigen::SparseMatrix<int> t = A.row(*o_itr);
    for (int k=0; k < t.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<int>::InnerIterator it(t,k); it; ++it)
        {
            two_ring.push_back(it.col());
        }
    }
  }

  P.resize(two_ring.size(),3);
  int p_idx = 0;
  for(std::vector<int>::iterator t_itr = two_ring.begin() ; t_itr != two_ring.end(); ++t_itr){
    int neighbor_idx = *t_itr;
    P.row(p_idx) = V.row(neighbor_idx) - V.row(v_idx);
    p_idx++;
  }
}
