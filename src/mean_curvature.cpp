#include "../include/mean_curvature.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <Eigen/SparseCholesky>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  //H = Eigen::VectorXd::Zero(V.rows());
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

  Eigen::SimplicialLDLT <Eigen::SparseMatrix<double> > solver;
  solver.compute(M);
  Eigen::SparseMatrix<double> I(M.rows(),M.rows());
  I.setIdentity();
  auto M_inv = solver.solve(I);

  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  H.resize(V.rows());
  Eigen::MatrixXd Hn = (M_inv * L * V);
  for(int i = 0; i < V.rows(); i++){
    double mag = Hn(i,0) / N(i,0);
    H(i) = mag; // TODO: reorient magnitude
  }
}
