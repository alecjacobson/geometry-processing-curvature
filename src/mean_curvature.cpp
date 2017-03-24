#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>


void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  H = Eigen::VectorXd::Zero(V.rows());

  Eigen::MatrixXd N;
  Eigen::SparseMatrix<double> L, M;

  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::per_vertex_normals(V, F, N);

  // do an element-wise inverse on the mass matrix:
  for(int32_t i = 0; i < M.outerSize(); i++)
  {
    for(Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
    {
      it.valueRef() = 1.0/it.value();
    }
  }

  // 𝑯 n ≈ H = ML⁻¹V
  Eigen::MatrixXd HVec = M*L*V;

  // get the signed mean curvature as the magnitude of each with direction: row*N
  for(int32_t i = 0; i < HVec.rows(); i++)
  {
    H(i) = -HVec.row(i).dot(N.row(i));
  }
}
