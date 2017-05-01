#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/pinv.h>
#include <igl/sort.h>
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
  
  Eigen::MatrixXd N(V.rows(), 3);
  igl::per_vertex_normals(V, F, N);
  
  for (int i = 0; i < num_v; ++i) {
    std::set<int> neighbours, neighbours2;
    
    // Find out neighbours
    for (Eigen::SparseMatrix<double>::InnerIterator it(adjacency, i); it; ++it) {
      neighbours.insert(it.row());
    }
    
    // Find out neighbours of neighbours
    std::set<int>::iterator s_it;
    for (s_it = neighbours.begin(); s_it != neighbours.end(); ++s_it) {
      for (Eigen::SparseMatrix<double>::InnerIterator it2(adjacency, *s_it); it2; ++it2) {
        if (it2.row() != i) {
          neighbours2.insert(it2.row());
        }
      }
    }
    
    // Combine neighbours and neighbours-of-neighbours
    neighbours.insert(neighbours2.begin(), neighbours2.end());
    
    // Construct P
    Eigen::MatrixXd P;
    P.resize(neighbours.size(), 3);
    int j = 0;
    for (s_it = neighbours.begin(); s_it != neighbours.end(); ++s_it) {
      P.row(j) = V.row(*s_it) - V.row(i);
      j++;
    }
    
    // Find basis for tangent plane: PCA was working horribly, so I'm doing this workaround
    Eigen::Matrix3d Q;
    Eigen::Vector3d normal = N.row(i);
    Eigen::Vector3d global_axis, tangent_axis;
    Q.col(2) = normal;
    if (normal[1] == 0 && normal[2] == 0) {
      // this might result in degeneracy so we take y direction
      global_axis << 0,1,0;
      tangent_axis = (Eigen::Matrix3d::Identity() - (normal * normal.transpose()) ) * global_axis;
      Q.col(1) = tangent_axis.normalized();
      Q.col(0) = Q.col(1).cross(Q.col(2));
    } else {
      global_axis << 1,0,0;
      tangent_axis = (Eigen::Matrix3d::Identity() - (normal * normal.transpose()) ) * global_axis;
      Q.col(0) = tangent_axis.normalized();
      Q.col(1) = Q.col(0).cross(Q.col(2));
    }
    
    // Now for each neighbour k in P, we find u, v, and w coefficients
    Eigen::MatrixXd new_P(P.rows(), 3); 
    new_P = (Q.transpose() * P.transpose()).transpose();
    Eigen::MatrixXd S(P.rows(), 2);
    Eigen::VectorXd B(P.rows());
    S.col(0) = new_P.col(0);
    S.col(1) = new_P.col(1);
    B = new_P.col(2);
    
    // Now we form the coefficients for a
    Eigen::MatrixXd coeff_a(P.rows(), 5), coeff_a_inv(5, P.rows());
    Eigen::VectorXd a(5);
    for (int k = 0; k < P.rows(); k++) {
      coeff_a(k,0) = S(k,0);
      coeff_a(k,1) = S(k,1);
      coeff_a(k,2) = S(k,0) * S(k,0);
      coeff_a(k,3) = S(k,0) * S(k,1);
      coeff_a(k,4) = S(k,1) * S(k,1);
    }
    igl::pinv(coeff_a, coeff_a_inv);
    a = coeff_a_inv * B;
    
    // Derive e, f, g, E, F, G
    double E = 1 + pow(a[0], 2);
    double F = a[0] * a[1];
    double G = 1 + pow(a[1], 2);
    double e = (2 * a[2]) / sqrt(pow(a[0], 2) + 1 + pow(a[1], 2));
    double f = a[3] / sqrt(pow(a[0], 2) + 1 + pow(a[1], 2));
    double g = (2 * a[4]) / sqrt(pow(a[0], 2) + 1 + pow(a[1], 2));
    
    // Construct S
    Eigen::Matrix2d s1, s2;
    s1 << e, f,
          f, g;
    s2 << E, F,
          F, G;
    S = (-s1 * s2.inverse()).transpose();
    
    // Eigen Decomposition of S
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(S);
    
    // Assign directions and curvatures
    Eigen::Vector2d curv_eigvals, curv_indices;
    igl::sort(es.eigenvalues(), 1, false, curv_eigvals, curv_indices);
    K1[i] = curv_eigvals[0];
    K2[i] = curv_eigvals[1];
    
    Eigen::VectorXd d1, d2;
    d1 = es.eigenvectors().col(curv_indices[0]);
    d2 = es.eigenvectors().col(curv_indices[1]);
    d1.conservativeResize(3, 1);
    d2.conservativeResize(3, 1);
    d1[2] = 0;
    d2[2] = 0;
    D1.row(i) = Q * d1;
    D2.row(i) = Q * d2;
  }
}
