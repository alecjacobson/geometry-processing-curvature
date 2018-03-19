#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    A.resize(l_sqr.rows(),3);
    for(int i = 0; i < 3; i++){
        int c = i;
        int b = (i+1)%3;
        int a = (i+2)%3;
        
        Eigen::ArrayXd numerator = (l_sqr.col(c)-l_sqr.col(a)-l_sqr.col(b)).array();
        Eigen::ArrayXd denominator = -2*(l_sqr.col(a).array()*l_sqr.col(b).array()).sqrt();
        Eigen::ArrayXd cosine = numerator/denominator;
        A.col(i) = cosine.acos().matrix();
    }
}
