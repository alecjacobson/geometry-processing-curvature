#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/invert_diag.h"
#include "igl/per_vertex_normals.h"

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
    Eigen::SparseMatrix<double> L, M, M_i;
    Eigen::MatrixXd curve_normals, normals;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, M_i);
    curve_normals = (M_i * L) * V;
    igl::per_vertex_normals(V, F, normals);
    
    H.resize(curve_normals.rows());
    for (int i = 0; i < V.rows(); i++) {
        H(i) = 0.0;
        if ((normals.row(i)).dot(curve_normals.row(i)) > 0) {
            H(i) = (curve_normals.row(i)).norm();
        }
        else if ((normals.row(i)).dot(curve_normals.row(i)) < 0) {
            H(i) = - ((curve_normals.row(i)).norm());
        }
    }
}
