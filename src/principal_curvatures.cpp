#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <vector>
#include <Eigen/Eigenvalues> 
#include <iostream>
#include <cmath>

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

  Eigen::SparseMatrix<int> adjacency_mat;
  igl::adjacency_matrix(F, adjacency_mat);
  // USEFUL PROPERTY of adjacency matrices!!
  // squaring the adjacency matrix gives walks of length 2.
  // since we're working with triangle meshes, all length-1 walks are also accessible via length-2 walks
  // so we only need to look at the squared adjacency matrix
  // https://en.wikipedia.org/wiki/Adjacency_matrix#Matrix_powers
  //adjacency_mat = adjacency_mat * adjacency_mat;

  for (int i = 0; i < V.rows(); i++) {
    std::vector<int> plane;
    Eigen::VectorXd neighbourhood = adjacency_mat.col(i);

    // find indices of all non-zero entries in current row of adjacency matrix
    for (int j = 0; j < neighbourhood.rows(); j++) {
      if (neighbourhood(j) > 0) {
        plane.push_back(j);
      }
    }
    // fill P with corresponding vertices
    int k = plane.size();
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(k, 3);
    for (int j = 0; j < k; j++) {
      P.row(j) = V.row(plane.at(j)) - V.row(i);
    }

    // eigen decomp on P'P
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    solver.compute(P.transpose() * P);
    Eigen::MatrixXd eigenvecs = solver.eigenvectors();

    // extract principal directions
    Eigen::VectorXd u_dir = eigenvecs.row(2);
    Eigen::VectorXd v_dir = eigenvecs.row(1);
    Eigen::VectorXd w_dir = eigenvecs.row(0);

    // get components for each point P
    Eigen::MatrixXd S(k, 2);
    Eigen::VectorXd B(k);
    for (int j = 0; j < k; j++) {
      S(j, 0) = P.row(j).dot(u_dir / u_dir.norm());
      S(j, 1) = P.row(j).dot(v_dir / v_dir.norm());
      B(j) = P.row(j).dot(w_dir / w_dir.norm());
    }

    // solve for 5 unknowns: w = a₁u + a₂v + a₃u² + a₄uv + a₅v²
    // Ax = b

    // compose our matrix of knowns
    Eigen::MatrixXd A(k, 5);
    A.col(0) = S.col(0);
    A.col(1) = S.col(1);
    A.col(2) = S.col(0).cwiseProduct(S.col(0));
    A.col(3) = S.col(0).cwiseProduct(S.col(1));
    A.col(4) = S.col(1).cwiseProduct(S.col(1));

    // compute inverse
    Eigen::MatrixXd A_inv;
    igl::pinv(A, A_inv);
    
    // solve for unknowns
    Eigen::VectorXd x;
    x = A_inv * B;

    // composing shape operator...
    double E = 1.0 + x(0) * x(0);
    double F = x(0) * x(1);
    double G = 1.0 + x(1) * x(1);
    double denom = sqrt( x(0) * x(0) + 1.0 + x(1) * x(1) );
    double e = ( 2.0 * x(2) ) / denom;
    double f = x(3) / denom;
    double g = ( 2.0 * x(4) ) / denom;

    Eigen::MatrixXd lowercase(2, 2);
    lowercase << e, f,
                 f, g;
    Eigen::MatrixXd uppercase(2, 2);
    uppercase << E, F,
                 F, G;
    
    Eigen::MatrixXd shape(2, 2);
    shape = - lowercase * uppercase.inverse(); // inv trivial on such a small matrix
    // eigen decomp on S to get k1 and k2
    // shape operator is self-adjoint
    solver.compute(shape);
    Eigen::MatrixXd tangents = solver.eigenvectors();
    Eigen::VectorXd curvatures = solver.eigenvalues();

    K1(i) = curvatures(0);
    K2(i) = curvatures(1);
    D1.row(i) = tangents(0, 0) * u_dir + tangents(0, 1) * v_dir;
    D2.row(i) = tangents(1, 0) * u_dir + tangents(1, 1) * v_dir;

    // "Lift k1 and k2 back to R3 coords"
  }

}
