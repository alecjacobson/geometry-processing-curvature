#include "../include/internal_angles.h"
#include <cmath>
#include <iostream>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    // Given
    // a^2 = b^2 + c^2 - 2*bc*cos(alpha)
    // alpha = acos( ( a^2 - b^2 - c^2 ) ) / -2*bc )
    //
    const int n = l_sqr.rows();
    A.resize( n, 3 );

    for( int i = 0; i < n; ++ i )
    {
        // for each triangle, a row is a set of edges (see igl::sqared_edge_lengths ;)
        Eigen::RowVector3d e( l_sqr.row(i) );
        for( int j = 0; j <3; ++j )
        {
            double a2 = e( j );
            double b2 = e( ( j + 1 ) % 3 );
            double c2 = e( ( j + 2 ) % 3 );
            double alpha = std::acos( ( a2 - b2 -c2 ) / ( -2 * std::sqrt( b2 * c2 ) ));
            
            A( i, j ) = alpha;
        }
    }
}
