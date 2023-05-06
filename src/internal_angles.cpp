#include "../include/internal_angles.h"

void internal_angles(const Eigen::MatrixXd & l_sqr, Eigen::MatrixXd & A) {
    
    A.resize(l_sqr.rows(), 3);
    
    // each row of l_sqr corresponds to a face, so loop through each
    // row and compute the internal angles of the face at each iteration
    for (int i = 0; i < l_sqr.rows(); i++) {
        
        // extract the squared edge lengths for the current face
        double a = l_sqr(i, 0);
        double b = l_sqr(i, 1);
        double c = l_sqr(i, 2);
        
        // compute the cosine of each angle using the law of cosines (i.e.
        // as in cotmatrix.cpp from the smoothing assignment)
        double cosa = (b + c - a)/(2*std::sqrt(b*c));
        double cosb = (a + c - b)/(2*std::sqrt(a*c));
        double cosc = (a + b - c)/(2*std::sqrt(a*b)); 
        
        // compute the arccosine of each angle to obtain the internal angles
        // and use them to populate the ith row of A
        A.row(i) = Eigen::Vector3d(std::acos(cosa), std::acos(cosb), std::acos(cosc));
    }
    
}
