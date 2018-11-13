#include "../include/mean_curvature.h"

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  int n = F.maxCoeff() + 1;
  H = Eigen::VectorXd::Zero(n);

  Eigen::SparseMatrix<double> L, M, M_inv;
  Eigen::MatrixXd N;

  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,M_inv);
  igl::per_vertex_normals(V,F,N);

  Eigen::MatrixXd Hn = M_inv * L * V;

  for (int i = 0; i < n; i++){
    double dotp = Hn.row(i).dot(N.row(i));
    if (dotp < 0){
      H(i) = - Hn.row(i).norm();
    }
    else{
      H(i) = Hn.row(i).norm();
    }
  }  
}
