#include "../include/angle_defect.h"
#include "../include/internal_angles.h"

#include <igl/doublearea.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
    const int n = V.rows();
    D.resize( n, 1 );
    for( int i = 0; i < n; ++i ) // gotta be a better way...
        D(i) = 2 * M_PI;

    Eigen::MatrixXd l_sqr( F.rows(), F.cols() );
    igl::squared_edge_lengths( V, F, l_sqr );
    
    Eigen::MatrixXd A;
    internal_angles( l_sqr, A );

    
    // Ki = ( 2*pi - Sum_over_faces_of_vert_i( theta_i_face ) ) ) / A_i
    // we can loop over all the faces and sum up their individual contributions
    // for each vertex on the triangle.
    const int nF = F.rows();
    for( int i = 0; i < nF; ++i )
    {
        Eigen::RowVector3i f( F.row(i) );
        Eigen::RowVector3d a( A.row(i) );
        for( int j = 0; j < 3; ++j )
        {
            int idx = f(j);
            D(idx) -= a(j);
        }
    }

    // // now divide each by the Ai (sum of areas for each face)
    // Eigen::VectorXd dblA( n );
    // igl::doublearea( V, F, dblA );

    // //dblA = 0.5 * dblA.inverse().transpose();
    // std::cout << "dblA: " << std::endl <<  dblA << std::endl;
    // // // again, there's probably a nicer Eigen way to do this....find it ;)
    // for( int i = 0; i < n; ++i )
    // {
    //     D(i) = D(i) / ( dblA(i) / 2.0 ) ;
    // }
}

