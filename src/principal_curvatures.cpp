#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <set>
#include <igl/polar_svd.h>
#include <igl/pinv.h>
#include <cmath>
#include <Eigen/Eigenvalues> 


void compute_shape_operator(double a1, double a2, double a3, double a4, double a5, Eigen::Matrix2d & S);

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

  // adjacency matrix
  Eigen::SparseMatrix<double> A;

  // get the mesh adjacency matrix
  igl::adjacency_matrix(F, A);


  // iterate on the vertices
  for (int i = 0; i < V.rows(); i++) {

    std::set<int> firstneighbors; // will hold only first order neighbors
    std::set<int> neighbors; // will hold 2 rings of neighbors

    // iterate only over the non zero elements of the adjacency matrix
    // use columns to identify the vertex (col-major)

    // iterate over the non-zero elems of column i
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
      if (it.row() == i)
        continue; // dont have self in the list right yet
      neighbors.insert(it.row());
      firstneighbors.insert(it.row());
    }

    // gather the second order neighbors
    // iterate over first order neighbors
    for (int firstorder : firstneighbors) {

      // iterate over 2nd degree neighbors
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, firstorder); it; ++it) {
        if (it.row() == firstorder)
          continue; // dont add again
        neighbors.insert(it.row());
      }
    }

    // now we have all neighbors vertex indices of the current vertex

    Eigen::MatrixXd P(neighbors.size(), 3);

    // make P
    int counter = 0;
    for (int vertex : neighbors) {
      P.row(counter) = V.row(vertex) - V.row(i);
      counter +=1;
    }

  // compute the PCA of P'P
  // dummy matrices for input to polar_svd
  Eigen::MatrixXd R_dummy, T_dummy, V_dummy;

  // singluar values
  Eigen::VectorXd S;

  // eigen vectors
  Eigen::MatrixXd U;

  Eigen::MatrixXd P_corr = P.transpose()*P;

  igl::polar_svd(P_corr, R_dummy, T_dummy, U, S, V_dummy);

  // identify u, v, w based on singular values in S
  Eigen::Vector3d u, v, w;

  // find the max/min
  for (int counter1 = 0; counter1 < 3; counter1++) {
    int idx1 = counter1;
    int idx2 = (counter1 + 1) % 3;
    int idx3 = (counter1 + 2) % 3;

    // max case
    if ((S(idx1) >= S(idx2)) && (S(idx1) >= S(idx3))) 
      u = U.col(idx1);
    else if ((S(idx1) <= S(idx2)) && (S(idx1) <= S(idx3))) // min case
      w = U.col(idx1);
    else
      v = U.col(idx1);
  }

  // now we have our u, v, w!
  // change the x, y, z coordinates of P into u, v, w!
  // if we have a i + b j + c k = p u + q v + r W, then p = (a i + b j + c k).dot(u/|u|)

  Eigen::MatrixXd P_uv(neighbors.size(), 2);
  Eigen::VectorXd S_w(neighbors.size());

  for (int p_idx=0; p_idx < P.rows(); p_idx++) {
    // u
    P_uv(p_idx, 0) = P.row(p_idx).dot(u/u.norm());
    // v
    P_uv(p_idx, 1) = P.row(p_idx).dot(v/v.norm());
    // w
    S_w(p_idx) = P.row(p_idx).dot(w/w.norm());
  }

  // fit the indices using igl::pinv
  // convert the problem to least squares fitting AX = W, where X has the unknowns
  // prepare the matrix A
  Eigen::MatrixXd A_uv(P_uv.rows(), 5);

  // fill in A_uv
  for (int auv_idx=0; auv_idx < P_uv.rows(); auv_idx++) {
    Eigen::VectorXd row_entry(5);

    double u_val = P_uv(auv_idx, 0);
    double v_val = P_uv(auv_idx, 1);

    row_entry(0) = u_val;
    row_entry(1) = v_val;
    row_entry(2) = u_val * u_val;
    row_entry(3) = u_val * v_val;
    row_entry(4) = v_val * v_val;

    A_uv.row(auv_idx) = row_entry;
  }

  // to get the unknowns as X = A.inv * W, get the pseudoinverse of A.
  Eigen::MatrixXd A_uv_inverse;

  // get the pseudoinverse
  igl::pinv(A_uv, A_uv_inverse);

  // X be the solution such that it is (a1, a2, a3, a4, a5)
  Eigen::VectorXd X_solution;
  X_solution = A_uv_inverse * S_w;

  Eigen::Matrix2d Shape_matrix;

  compute_shape_operator(X_solution(0), X_solution(1), X_solution(2), X_solution(3), X_solution(4), Shape_matrix);


  // eigen decompose the Shape matrix!

  // singluar values
  Eigen::Vector2d Shape_eigenvals;

  // eigen vectors
  Eigen::Matrix2d Shape_eigenvectors;

  // Shape matrix would be symmetric so can use self-adjoint solver.
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver1(Shape_matrix);
  Shape_eigenvals = solver1.eigenvalues();
  Shape_eigenvectors = solver1.eigenvectors();

  // get the max and min curvatures and finish!

  if (abs(Shape_eigenvals(0)) >= abs(Shape_eigenvals(1))) {

    K1(i) = Shape_eigenvals(0);
    K2(i) = Shape_eigenvals(1);

    D1.row(i) = Shape_eigenvectors.col(0)(0) * u +
                Shape_eigenvectors.col(0)(1) * v;
    D1.row(i) /= D1.row(i).norm();

    D2.row(i) = Shape_eigenvectors.col(1)(0) * u +
                Shape_eigenvectors.col(1)(1) * v;
    D2.row(i) /= D2.row(i).norm();

  } else {
    K1(i) = Shape_eigenvals(1);
    K2(i) = Shape_eigenvals(0);

    D1.row(i) = Shape_eigenvectors.col(1)(0) * u +
                Shape_eigenvectors.col(1)(1) * v;
    D1.row(i) /= D1.row(i).norm();

    D2.row(i) = Shape_eigenvectors.col(0)(0) * u +
                Shape_eigenvectors.col(0)(1) * v;
    D2.row(i) /= D2.row(i).norm();

  }

  }
}

void compute_shape_operator(double a1, double a2, double a3, double a4, double a5, Eigen::Matrix2d & S) {

  double E, F, G, e, f, g;

  E = 1 + a1*a1;
  F = a1*a2;
  G = 1 + a2*a2;

  e = (2.0 * a3)/sqrt(a1*a1 + 1.0 + a2*a2);
  f = a4/sqrt(a1*a1 + 1.0 + a2*a2);
  g = (2.0 * a5)/sqrt(a1*a1 + 1.0 + a2*a2);

  // form the two matrices as given in the notes - call them A and B
  
  Eigen::Matrix2d A, B;

  A(0,0) = e;
  A(0,1) = f;
  A(1,0) = f;
  A(1,1) = g;

  B(0,0) = E;
  B(0,1) = F;
  B(1,0) = F;
  B(1,1) = G;

  S = -1.0 * A * B.inverse(); // product of two symm matrices or inv of a symm matrix is symm
}
