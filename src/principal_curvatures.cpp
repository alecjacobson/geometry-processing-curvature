#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <Eigen/Eigenvalues>
#include <set>
#include <iostream>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  int num_v = V.rows();
  D1.resize(num_v, 3);
  D2.resize(num_v, 3);
  K1.resize(num_v);
  K2.resize(num_v);
  
  Eigen::SparseMatrix<double> adjacency(num_v, num_v);
  igl::adjacency_matrix(F, adjacency);
  
  Eigen::MatrixXd P;
  Eigen::EigenSolver<Eigen::MatrixXd> es;
  std::set<int> neighbours, neighbours2;
  for (int i = 0; i < num_v; i++) {
    // Find out neighbours
    for (Eigen::SparseMatrix<double>::InnerIterator it(adjacency, i); it; ++it) {
      neighbours.insert(it.index());
    }
    std::cout << "got 1st neighbours" << std::endl;
    
    // Find out neighbours of neighbours
    std::set<int>::iterator s_it;
    for (s_it = neighbours.begin(); s_it != neighbours.end(); ++s_it) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(adjacency, *s_it); it; ++it) {
        neighbours2.insert(it.index());
      }
    }
    std::cout << "got 2nd neighbours" << std::endl;
    
    // Combine neighbours and neighbours-of-neighbours
    neighbours.insert(neighbours2.begin(), neighbours2.end());
    
    // Construct P
    P.resize(neighbours.size(), 3);
    for (int j = 0; j < neighbours.size(); j++) {
      P.row(j) = V.row(j) - V.row(i);
    }
    std::cout << "constructed P" << std::endl;
    
    // Perform PCA on P
    es.compute(P.transpose() * P);
    
    // Now for each neighbour k in P, we find u, v, and w coefficients
    Eigen::MatrixXd new_P(P.rows(), 3), S(P.rows(), 2);
    Eigen::VectorXd B(P.rows());
    Eigen::Matrix3d Q = es.eigenvectors().real();
    new_P = P * Q;
    S.col(0) = new_P.col(0);
    S.col(1) = new_P.col(1);
    B = new_P.col(2);
    std::cout << "set up S and B" << std::endl;
    
    // Now we form the coefficients for a
    Eigen::MatrixXd coeff_a(P.rows(), 5), coeff_a_inv(5, P.rows());
    Eigen::VectorXd a(5);
    coeff_a.col(0) = S.col(0);
    coeff_a.col(1) = S.col(1);
    coeff_a.col(2) = (S.col(0).array().square()).matrix();
    coeff_a.col(3) = (S.col(0).array() * S.col(1).array()).matrix();
    coeff_a.col(4) = (S.col(1).array().square()).matrix();
    igl::pinv(coeff_a, coeff_a_inv);
    a = coeff_a_inv * B;
    std::cout << "solved for a" << std::endl;
    
    // Derive e, f, g, E, F, G
    double E = 1 + pow(a[0], 2);
    double F = a[0] * a[1];
    double G = 1 + pow(a[1], 2);
    double e = (2 * a[2]) / sqrt(pow(a[0], 2) + 1 + pow(a[1], 2));
    double f = a[3] / sqrt(pow(a[0], 2) + 1 + pow(a[1], 2));
    double g = (2 * a[4]) / sqrt(pow(a[0], 2) + 1 + pow(a[1], 2));
    std::cout << "derived efg" << std::endl;
    
    // Construct S
    Eigen::Matrix2d s1, s2;
    s1 << e, f,
          f, g;
    s2 << E, F,
          F, G;
    S = -s1 * s2.inverse();
    std::cout << "got S" << std::endl;
    
    // Eigen Decomposition of S
    es.compute(S);

    // Assign directions and curvatures
    Eigen::Matrix3d eigvecs;
    eigvecs.setZero();
    eigvecs.block(0,0,2,2) = es.eigenvectors().real();
    
    Eigen::Matrix3d final_eigvecs = eigvecs * Q.transpose();
    D1.row(i) = final_eigvecs.col(0);
    D2.row(i) = final_eigvecs.col(1);
    K1[i] = es.eigenvalues().real()[0];
    K2[i] = es.eigenvalues().real()[1];
    std::cout << "finished" << std::endl;
  }
}
