#include "../include/internal_angles.h"
#include <iostream>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  int num_faces = l_sqr.rows();
  A.resize(num_faces, 3);
  for(int f_idx = 0; f_idx < num_faces; f_idx++){
    for(int i=0; i < 3; i++){
      double a_squared = l_sqr(f_idx,i);
      double b_squared = l_sqr(f_idx,(i+1)%3);
      double c_squared = l_sqr(f_idx,(i+2)%3);
      // acos is in radians
      A(f_idx, i) = acos(
                    (b_squared + c_squared - a_squared)
                    / (2 * sqrt(c_squared) * sqrt(b_squared)));
    }
  }
}
