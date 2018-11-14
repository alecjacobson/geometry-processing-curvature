#include "../include/mean_curvature.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
    // Compute mean curvature normal `Hn = M^{-1}*L*V`
    //      by solving  M*Hn = L*V

    Eigen::SparseMatrix<double> L, M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    Eigen::MatrixXd b;
    b = L*V;
    Eigen::MatrixXd Hn(V.rows(), 3);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(M);
    Hn = solver.solve(b);

    // Determine sign
     H.resize(V.rows());

    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);

    for (int i = 0; i < V.rows(); ++i) {
        H(i) = Hn.row(i).norm();
        if (Hn.row(i).dot(N.row(i)) > 0) {
            H(i) = -H(i);
        }
    }
}

