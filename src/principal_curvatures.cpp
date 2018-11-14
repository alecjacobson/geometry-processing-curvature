#include "../include/principal_curvatures.h"
#include <igl/adjacency_list.h>
#include <igl/per_vertex_normals.h>
#include <igl/eigs.h>
#include <igl/pinv.h>
#include <igl/sort.h>
#include <igl/colon.h>
#include <igl/slice.h>

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
  D1.resize(V.rows(),3);
  D2.resize(V.rows(),3);
  
  std::vector<std::vector<int>> adj;
  igl::adjacency_list(F, adj);
  
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);
  
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  for (int i = 0; i < V.rows(); i++) {
  
    // gather the vertices in the two-ring of vi into a matrix P
    std::vector<int> two_ring;
    for (int u : adj[i]) {
      if (find(two_ring.begin(), two_ring.end(), u) == two_ring.end())
        two_ring.push_back(u);
      for (int w : adj[u]) {
        if (find(two_ring.begin(), two_ring.end(), w) == two_ring.end())
          two_ring.push_back(w);
      }
    }
    
    int k = two_ring.size();
    Eigen::MatrixXd P(k, 3);
    for (int j = 0; j < k; j++) {
      P.row(j) = V.row(two_ring[j]) - V.row(i);
    }
    
    // do PCA on Q = P'P
    es.compute(P.transpose() * P);
    Eigen::VectorXd Q_evals, Q_inds;
    igl::sort(es.eigenvalues(), 1, false, Q_evals, Q_inds);
    Eigen::MatrixXd Q_evecs;
    igl::slice(es.eigenvectors(), Q_inds, 2, Q_evecs);
    
    // w should be consistently oriented
    Q_evecs.col(2) *= Q_evecs.col(2).dot(N.row(i)) > 0 ? 1 : -1;
    
    // [Sc|B] = (Q_evecs'*P')' = P*Q_evecs
    Eigen::MatrixXd Q_uv(3, 2);
    Q_uv << Q_evecs.col(0), Q_evecs.col(1);
    Eigen::MatrixXd Sc = P * Q_uv;
    Eigen::VectorXd B = P * Q_evecs.col(2);
    
    // use least squares to fit a quadratic function to the k points (in the above basis)
    Eigen::MatrixXd A(k, 5);
    for (int i = 0; i < Sc.rows(); i++) {
      double u = Sc(i, 0);
      double v = Sc(i, 1);
      A.row(i) << u, v, u*u, u*v, v*v;
    }
    Eigen::VectorXd quad_coeffs = (A.transpose() * A).ldlt().solve(A.transpose() * B);
    double a1 = quad_coeffs[0];
    double a2 = quad_coeffs[1];
    double a3 = quad_coeffs[2];
    double a4 = quad_coeffs[3];
    double a5 = quad_coeffs[4];
    
    // compute shape operator
    Eigen::Matrix2d F1, F2, S;
    F1 << 1 + a1*a1, a1*a2, a1*a2, 1+a2*a2;
    F2 << 2*a3, a4, a4, 2*a5;
    F2 /= sqrt(a1*a1+1+a2*a2);
    
    S = -F2 * F1.inverse();
    
    // do PCA on S
    es.compute(S);
    Eigen::VectorXd S_evals, S_inds;
    igl::sort(es.eigenvalues(), 1, false, S_evals, S_inds);
    Eigen::MatrixXd S_evecs;
    igl::slice(es.eigenvectors(), S_inds, 2, S_evecs);
    
    // finally we have the principal curvatures and directions
    K1[i] = S_evals[0];
    K2[i] = S_evals[1];
    
    D1.row(i) = Q_uv * S_evecs.col(0);
    D2.row(i) = Q_uv * S_evecs.col(1);
  }
}
