#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    Eigen::VectorXd &H)
{
  //H = Eigen::VectorXd::Zero(V.rows());
  Eigen::SparseMatrix<double> L, M, M_in;

  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI, M);
  igl::invert_diag(M, M_in);

  //calculate Hn(mean curvature normal)
  Eigen::MatrixXd Hn, N;
  Hn = M_in * L * V;

  //get per vertex normal
  igl::per_vertex_normals(V, F, N);

  H.resize(V.rows());

  //match normal
  for (int i = 0; i < V.rows(); i++)
  {
    double p = Hn.row(i).dot(N.row(i));
    if (p >= 0)
      H(i) = Hn.row(i).norm();
    //reverse if negative
    else
      H(i) = -1.00 * Hn.row(i).norm();
  }
}
