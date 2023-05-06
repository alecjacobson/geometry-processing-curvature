#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <Eigen/SparseCholesky>	
#include <Eigen/SparseLU>	
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/per_vertex_normals.h>
#include <iostream>
#include <igl/invert_diag.h>

// sign function
double sgn(double x);

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  H = Eigen::VectorXd::Zero(V.rows());

  // get the mass matrix
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

  // get the cotmatrix
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::MatrixXd H_vectors(V.rows(), 3);

  // solve using sparse solver (commented) - Instead invert the diag elems.
  // Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
  // solver.compute(M);
  // H_vectors = solver.solve(L*V);

  // since M is diagonal, do this instead of solving by solvers
  Eigen::SparseMatrix<double> M_inv;
  igl::invert_diag(M, M_inv);
  H_vectors = M_inv * L * V;


  // compute the per vertex normals to get the sign
  Eigen::MatrixXd V_normals(V.rows(), 3);
  igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, V_normals);

  // fill in the signed curvature!
  for (int i = 0; i < V.rows(); i++) {
    H(i) = H_vectors.row(i).norm() * sgn(V_normals.row(i).dot(H_vectors.row(i)));
    //std::cout<<H(i)<<std::endl;
  }

}

// sign function 
double sgn(double x) {
  if (x > 0.0) 
    return +1.0;
  else if (x < 0.0)
    return -1.0;
  else
    return 0;
}
