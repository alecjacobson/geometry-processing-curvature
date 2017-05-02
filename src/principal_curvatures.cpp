#include "../include/principal_curvatures.h"
#include "../include/internal_angles.h"

#include <igl/squared_edge_lengths.h>
#include <igl/principal_curvature.h>
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <igl/per_vertex_normals.h>

#include <vector>
#include <set>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
    const int n = V.rows();

    Eigen::MatrixXd N;
    igl::per_vertex_normals( V, F, N );
    
    K1 = Eigen::VectorXd::Zero(V.rows());
    K2 = Eigen::VectorXd::Zero(V.rows());
    D1 = Eigen::MatrixXd::Zero(V.rows(),3);
    D2 = Eigen::MatrixXd::Zero(V.rows(),3);

    Eigen::SparseMatrix< double > Madj;
    igl::adjacency_matrix( F, Madj );
    
    // for each vertex, 
    // get 2 ring around vertex - results in k points
    // build P^kx3 for that vert
    // define quadratic surface:
    //   - PCA of P -> Eigen Decomp of P^T * P
    //   - S e R^kx2 are the u-v directions for each point in P
    //   - B e R^k (eigenvalues) is the height at each point    
    for( int i=0; i<n; ++i )
    {
        // find the 2-ring verts -- adjacent verts and their shared edge verts.
        // Find the # of non-zero (adjacent) verts in the column (InnerNNZ's)
        // TODO - assumes no degenerate edges or non-unique edges...
        std::set< int > twoRingVerts;
        std::vector< int > oneRingVerts;
        for( Eigen::SparseMatrix<double>::InnerIterator it( Madj, i ); it; ++it )
        {
            oneRingVerts.push_back( it.index() );
            twoRingVerts.insert( it.index() );
        }

        // now add the 2-ring verts...
        const int nOneRing = oneRingVerts.size();
        for( int j = 0; j < nOneRing; ++j )
        {
            for( Eigen::SparseMatrix<double>::InnerIterator it( Madj, oneRingVerts[j] ); it; ++it )
                twoRingVerts.insert( it.index() );
        }
        
        Eigen::MatrixXd P( twoRingVerts.size(), 3 ); // 2-ring
        int count( 0 );
        for( auto &&v : twoRingVerts )
            P.row( count++ ) = V.row( v );

        // get the distance from our vert for each vert in the 2-ring (turing?)
        P.rowwise() -= V.row( i ); 

        // Find the eigen vectors and values of P^T * P
        Eigen::Matrix3d PTP = P.transpose() * P;
        Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > pca( PTP );
        
        Eigen::MatrixXd s( 2, 3 );
        s = pca.eigenvectors().block( 0, 1, 3, 2 ).transpose();
        Eigen::MatrixXd Ps = P * s.transpose();

        // eigenvalues give our height above the surface, B
        Eigen::RowVector3d w = pca.eigenvectors().col( 0 ).transpose();

        // check with normals if above/below surface...
        if( N.row( i ) * w.transpose() < 0.0 )
            w = -w;

        Eigen::MatrixXd B =  P * w.transpose();

        //
        // Least Squares
        //
        // for k sets of u,v,w do least squares and fit for 5 unknowns below
        // w = a_1*u + a_2*v + a_3*u^2 + a_4*u*v + a_5*v^2
        //
        // To avoid confusion with indexes, rewrite as:
        // w = a_0*u + a_1*v + a_2*u^2 + a_3*u*v + a_4*v^2
        //
        // unknowns in columns a_0 through a_4 and use igl::pinv
        // Solve for Ainv given A for our system in P on S (the coefficients).
        //
        Eigen::MatrixXd A( Ps.rows(), 5 );
        Eigen::MatrixXd Ainv; // A^-1
        A.col( 0 ) = Ps.col( 0 );                               // a_0*u
        A.col( 1 ) = Ps.col( 1 );                               // a_1*v
        A.col( 2 ) = Ps.col( 0 ).array() * Ps.col( 0 ).array(); // a_2*u^2
        A.col( 3 ) = Ps.rowwise().prod();                       // a_3*u*v
        A.col( 4 ) = Ps.col( 1 ).array() * Ps.col( 1 ).array(); // a_4*v^2

        igl::pinv( A, Ainv );

        Eigen::VectorXd Ac = Ainv * B;

        // Assemble our Shape operator...osculating jets are interesting....
        double E = 1 + Ac(0) * Ac(0); // 1 + a_0^2
        double F = Ac(0) * Ac(1);     // a_0 * a_1
        double G = 1 + Ac(1) * Ac(1); // 1 + a_1^2

        double denom = std::sqrt( Ac(0)*Ac(0) + 1 + Ac(1)*Ac(1) );
        double e = ( 2 * Ac(2) ) / denom;
        double f = (     Ac(3) ) / denom;
        double g = ( 2 * Ac(4) ) / denom;
        

        // S is the product of the two 2x2 matrices, let's call them lowercase ( e, f, f, g )
        // and uppercase ( E, F, F, G ) 
        Eigen::MatrixXd lowercase( 2, 2 ), uppercase( 2, 2 );
        lowercase << e, f,
                     f, g;

        uppercase << E, F,
                     F, G;

        // Now do the eigen decomposition of shape operator S:
        // get the principal curvatures k1, k2 and principal tangent directions       
        Eigen::Matrix2d S = -lowercase * uppercase.inverse();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> pca_on_S( S );

        Eigen::Vector2d principalCurvatures = pca_on_S.eigenvalues();

        Eigen::Matrix< double, 2, 3 > principalDirs = ( s.transpose() * pca_on_S.eigenvectors() ).transpose();

        // distinguish between min and max components...
        int min( 0 ), max( 1 );
        if( principalCurvatures( min ) > principalCurvatures( max ) )
        {
            max = 0;
            min = 1;
        }

        K1( i ) = principalCurvatures( min );
        K2( i ) = principalCurvatures( max );

        D1.row( i ) = principalDirs.row( min );
        D2.row( i ) = principalDirs.row( max );
    }
}
