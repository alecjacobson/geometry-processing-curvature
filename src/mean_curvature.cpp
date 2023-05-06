#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/per_vertex_normals.h"
#include "igl/invert_diag.h"

void mean_curvature(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::VectorXd &H) {

    H = Eigen::VectorXd::Zero(V.rows());

    //Computing M
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    //Computing L
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    // computing the per vertex normals
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);


    Eigen::SparseMatrix<double> M_inv;
    //this should be fine here.
    igl::invert_diag(M, M_inv);


    Eigen::MatrixXd hn = M_inv * L * V;
    //getting final H
    for (int i = 0; i < V.rows(); i++) {
        double t = N.row(i).dot(hn.row(i));
        H(i) = hn.row(i).norm();
        if (t < 0) {
            H(i) = -H(i);
        }
    }

}
