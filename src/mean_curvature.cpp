#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/per_vertex_normals.h"
#include "igl/invert_diag.h"

double sign_func(double x) //Copied from internet
{
    if (x > 0)
        return +1.0;
    else if (x == 0)
        return 0.0;
    else
        return -1.0;
}


void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
    Eigen::SparseMatrix<double> L, M, M_inv;
    Eigen::MatrixXd N;
    
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::per_vertex_normals(V,F,N);
    igl::invert_diag(M,M_inv);
    
    Eigen::MatrixXd MinvLV = (M_inv * L) * V;
    Eigen::ArrayXd dot_prod = (MinvLV.array()*N.array()).rowwise().sum();
    
    dot_prod = dot_prod.unaryExpr(std::ptr_fun(sign_func)); //Copied from internet; .sign() apparently doesn't exist in the version we're using
    
    Eigen::ArrayXd norms = (MinvLV.array()*MinvLV.array()).rowwise().sum().sqrt();
    H = -(norms*dot_prod).matrix();
}
