#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  A.resize(l_sqr.rows(),3);
  double a,b,c;
  for (int i=0; i<l_sqr.rows(); i++){
  	for (int j=0; j<3; j++){
  	  a=l_sqr(i,j);
  	  b=l_sqr(i,(j+1)%3);
  	  c=l_sqr(i,(j+2)%3);
  	  A(i,j)=acos((b+c-a)/(2*sqrt(b)*sqrt(c)));
	}
  }
}
