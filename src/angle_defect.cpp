#include "../include/angle_defect.h"
#include "../include/internal_angles.h"

#include <igl/adjacency_matrix.h>
#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D) {
    
    D.resize(V.rows());

    // compute the squared edge lengths for each face
    Eigen::MatrixXd l_sqr;
    igl::squared_edge_lengths(V, F, l_sqr);
    
    // use the squared edge lengths to compute the internal angles for each face
    Eigen::MatrixXd in_angle;
    internal_angles(l_sqr, in_angle);
    
    // loop through each face and at each iteration subtract the internal angles from their 
    // corresponding rows in D (where each row corresponds to a vertex).
    D = Eigen::VectorXd::Constant(V.rows(), 2*M_PI);
    
    for (int f = 0; f < F.rows(); f++) {

        // update the appropriate rows
        D(F(f, 0)) -= in_angle(f, 0);
        D(F(f, 1)) -= in_angle(f, 1);
        D(F(f, 2)) -= in_angle(f, 2);
    }
        
}
