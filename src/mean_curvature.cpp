#include "../include/mean_curvature.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/invert_diag.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  Eigen::SparseMatrix<double> L,M,M1;
  Eigen::MatrixXd Hn,N;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::per_vertex_normals(V, F, N);
  igl::invert_diag(M,M1);
  Hn=M1*L*V;
  H.resize(V.rows());
  for (int i=0; i<V.rows(); i++){
  	H(i)=Hn.row(i).dot(N.row(i))/(N.row(i).norm());
  }
}
