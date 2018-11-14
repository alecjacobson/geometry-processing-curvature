#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/sort.h>
#include <igl/pinv.h>
#include <igl/per_vertex_normals.h>

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
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);

  // Vertex normals
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  // Construct adjacency matrix
  Eigen::SparseMatrix<int> adj;
  igl::adjacency_matrix(F, adj);

  // Useful containers
  std::vector<int> p_candidates;
  std::vector<int> p_sorted;
  std::vector<size_t> dummy;

  // loop through each vertex
  for (int i = 0; i < V.rows(); ++i) {
    // Construct P *************************************************************

    p_candidates.clear();
    p_sorted.clear();

    // Inner ring of v
    for (Eigen::SparseMatrix<int>::InnerIterator it(adj, i); it; ++it) {
      p_candidates.push_back(it.row());
    } // end loop it

    // Outer ring of v
    int r1_size = p_candidates.size();
    for (int ii = 0; ii < r1_size; ++ii) {
      for (Eigen::SparseMatrix<int>::InnerIterator it(adj, p_candidates[ii]);
          it; ++it) {
        p_candidates.push_back(it.row());
      } // end loop it
    } // end loop ii

    // Sort for purpose of easily identifying duplicates
    igl::sort(p_candidates, true, p_sorted, dummy);

    Eigen::MatrixXd P(p_sorted.size(), 3);
    int k = 0;
    P.row(k) = V.row(p_sorted[0]) - V.row(i);
    ++k;
    for (int ii = 1; ii < p_sorted.size(); ++ii) {
      if (p_sorted[ii] != p_sorted[ii - 1] && p_sorted[ii] != i) {
        P.row(k) = V.row(p_sorted[ii]) - V.row(i);
        ++k;
      }
    } // end loop ii
    P.conservativeResize(k, 3);

    // Fit quadratic ************************************************************
    
    // Construct projection
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(P.transpose()*P);
    Eigen::MatrixXd eig_vec = solver.eigenvectors();
    Eigen::VectorXd eig_val = solver.eigenvalues();

    // Set to be consistent with vertex surface normal
    if (N.row(i).dot(eig_vec.col(0).transpose()) < 0) {
      eig_vec *= -1.0;
    }
    Eigen::MatrixXd proj = P*eig_vec;

    // Fit quadratic
    Eigen::MatrixXd U(k, 5);
    for (int ii = 0; ii < k; ++ii) {
      U(ii,0) = proj(ii,2);
      U(ii,1) = proj(ii,1);
      U(ii,2) = proj(ii,2)*proj(ii,2);
      U(ii,3) = proj(ii,2)*proj(ii,1);
      U(ii,4) = proj(ii,1)*proj(ii,1);
    } // end loop ii
    Eigen::MatrixXd U_pinv;
    igl::pinv(U, U_pinv);
    Eigen::VectorXd a = U_pinv*proj.col(0);

    // Compute principals *******************************************************

    double denom = 1.0/std::sqrt(a(0)*a(0) + 1.0 + a(1)*a(1));
    Eigen::MatrixXd temp1(2,2);
    temp1(0,0) = 2.0*a(2)*denom;
    temp1(0,1) = a(3)*denom;
    temp1(1,0) = temp1(0,1);
    temp1(1,1) = 2.0*a(4)*denom;
    Eigen::MatrixXd temp2(2,2);
    temp2(0,0) = 1.0 + a(0)*a(0);
    temp2(0,1) = a(0)*a(1);
    temp2(1,0) = temp2(0,1);
    temp2(1,1) = 1.0 + a(1)*a(1);
    Eigen::MatrixXd S = -temp1*temp2.inverse();

    // Eigendecomposition of S
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solverS(S);
    Eigen::MatrixXd eig_vecS = solverS.eigenvectors();
    Eigen::VectorXd eig_valS = solverS.eigenvalues();

    // Set principal curvatures
    K1(i) = eig_valS(1);
    K2(i) = eig_valS(0);

    // Lift principal tangent directions back to R3
    Eigen::MatrixXd lift(3,2);
    lift.col(0) = eig_vec.col(2);
    lift.col(1) = eig_vec.col(1);
    D1.row(i) = (lift*eig_vecS.col(1)).transpose();
    D2.row(i) = (lift*eig_vecS.col(0)).transpose();

  } // end loop i
}
