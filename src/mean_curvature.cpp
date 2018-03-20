#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H) {
    
    // construct the cotangent Laplacian.
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V,F,L);
    
    // construct the mass matrix.
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    
    // construct the mean curvature using a pointwise approximation of the mean curvature
    // normal. Because the mass matrix is diagonal, compute its inverse by simply 
    // inverting each of the entries along the diagonal.
    Eigen::SparseMatrix<double> Minv;
    igl::invert_diag(M, Minv);
    
    Eigen::MatrixXd Hn = Minv*L*V;
    
    // strip the magnitude off each row of Hn to obtain the unsigned mean curvature
    H = Hn.rowwise().norm();
    
    // preserve the sign by checking orientation against the per vertex normals
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);
    
    for (int v = 0; v < V.rows(); v++) {
        if (N.row(v).dot(Hn.row(v)) > 0.0) {
            H(v) *= -1;
        }
    }
}
