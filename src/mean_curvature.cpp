#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
    Eigen::MatrixXd N;
    Eigen::SparseMatrix<double> L,M;
    H = Eigen::VectorXd::Zero(V.rows());
    igl::cotmatrix(V,F,L);
    L = -L;
    igl::per_vertex_normals(V,F,N);

    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT, M);
    for(int i = 0; i < M.rows(); ++i) {
        auto&& v = M.coeffRef(i,i);
        v = 1.0/v;
    }
    H = ((.5 * M * L * V).array() * N.array()).rowwise().sum();
    std::cout << "Mean curvature!" << std::endl;
    std::cout << H.transpose() << std::endl;
}
