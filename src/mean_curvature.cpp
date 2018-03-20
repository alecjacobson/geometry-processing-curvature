#include "../include/mean_curvature.h"
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  H.resize(V.rows());

  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

  Eigen::MatrixXd Hn(V.rows(), 3);
  Hn = M.cwiseInverse() * L * V;

  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  for (int i = 0; i < V.rows(); i++)
  {
	  Eigen::RowVector3d c_n = Hn.row(i);
	  Eigen::RowVector3d vertex_n = N.row(i);
	  float dot = c_n.dot(vertex_n);
	  if (signbit(dot))
		  H(i) = c_n.norm();
	  else
		  H(i) = -c_n.norm();
  }

}
