#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>

void mean_curvature(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::VectorXd & H)
{
    Eigen::SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    igl::invert_diag(M, Minv);

    Eigen::MatrixXd Hn = Minv * L * V;
    H = Hn.rowwise().norm();
}
