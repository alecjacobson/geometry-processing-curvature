#include "../include/mean_curvature.h"
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  Eigen::SparseMatrix<double> L; L.setZero();
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> M; M.setZero();
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, M);

  Eigen::MatrixXd N; N.setZero();
  igl::per_vertex_normals(V, F, N);
   
  Eigen::MatrixXd Hn = M*L*V;
  H.resize(Hn.rows());

  for (int i = 0; i < Hn.rows(); i++)
  {
    double dot = Hn.row(i).dot(N.row(i));
    H[i] = Hn.row(i).norm()*(dot /abs(dot));
  }
    
}
