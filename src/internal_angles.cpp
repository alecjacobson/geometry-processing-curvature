#include "../include/internal_angles.h"
#include <igl/cotmatrix.h>
#include <iostream>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    A.resize(l_sqr.rows(), 3);

    for (int t = 0; t < l_sqr.rows(); t++){

        double a = l_sqr(t, 0);
        double b = l_sqr(t, 1);
        double c = l_sqr(t, 2);

        double cos_aij = ((a*a) + (b*b) - (c*c))/(2*a*b);
        double cos_aik = ((c*c) + (b*b) - (a*a))/(2*c*b);
        double cos_ajk = ((c*c) + (a*a) - (b*b))/(2*c*a);

        A(t, 0) = acos(cos_aij);
        A(t, 1) = acos(cos_aik);
        A(t, 2) = acos(cos_ajk);
    }
  
}
