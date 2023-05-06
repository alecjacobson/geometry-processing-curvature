#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
    const Eigen::MatrixXd & l_sqr,
    Eigen::MatrixXd & A)
{
    A.resize(l_sqr.rows(), l_sqr.cols());
    for (int i = 0; i < l_sqr.rows(); i++) 
    {
        Eigen::VectorXd l = l_sqr.row(i);
        for (int j = 0; j < l.size(); j++) {
            
            // cosine law
            double c_sqrd = l[j];
            double b_sqrd = l[(j + 1) % 3];
            double a_sqrd = l[(j + 2) % 3];
            
            A(i,j) = acos((c_sqrd - a_sqrd - b_sqrd) / (-2 * sqrt(a_sqrd * b_sqrd)));
        }
    }
}
