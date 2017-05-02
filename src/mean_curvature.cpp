#include "../include/mean_curvature.h"

#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

// Compute the discrete mean curvature at each vertex of a mesh (V,F) by taking
// the signed magnitude of the mean curvature normal as a pointwise (or
// integral average) quantity.
//
void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
    const int n = V.rows();
  
    // Replace with your code
    H = Eigen::VectorXd::Zero(V.rows());
    Eigen::MatrixXd bigH( n, 3 );
    
    // H = Minv * L * V, and H is of size R^nx3
    Eigen::SparseMatrix<double> M;
    igl::massmatrix( V, F, igl::MASSMATRIX_TYPE_DEFAULT, M );

    // invert by solving: (could also just invert along the diagonal...)
    // Ax = b -> Mx = I
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double > > solver;
    solver.compute( M );
    if( solver.info() != Eigen::Success )
        std::cout<<"SimplicialLDLT failed to prefactor M, might crash..."<<std::endl;

    Eigen::SparseMatrix<double> I( n, n );
    I.setIdentity();
    auto Minv = solver.solve( I );
    if( solver.info() != Eigen::Success )
        std::cout<<"SimplicialLDLT failed to SOLVE M, might crash..."<<std::endl;

    // Lapacian
    Eigen::SparseMatrix<double> L( n, n );
    igl::cotmatrix( V, F, L );

    bigH = Minv * L * V;

    Eigen::MatrixXd N;
    igl::per_vertex_normals( V, F, N );
    H = 0.5 * bigH.rowwise().norm();
    for( int i = 0; i < n; ++i )
    {
        Eigen::RowVector3d c( bigH.row(i) );
        Eigen::RowVector3d n( N.row(i) );
        c.normalize();
        n.normalize();
        double d = c.dot( n );

        if( d < 0.0 )
            H(i) = -1*H(i);            
    }
}
