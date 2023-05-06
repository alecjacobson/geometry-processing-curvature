#include "../include/internal_angles.h"
#include <cmath>
#include <iostream>
void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  //l_sr: [1, 2], [2, 0], [0, 1]
  // A:  angles that is opposite of each edge
  A.resize(l_sqr.rows(), l_sqr.cols());
  for (int i = 0; i< l_sqr.rows(); i++){
  	for (int j = 0; j < l_sqr.cols(); j++){
  		auto a2 = l_sqr(i, (j + 2) % l_sqr.cols());
  		auto b2 = l_sqr(i, (j + 1) % l_sqr.cols());
  		auto c2 = l_sqr(i, j);
  		//cosine rule
  		auto cosA = (a2 + b2 - c2 )/(2 *sqrt(a2 * b2));
  		A(i, j) = acos(cosA);
  	}
  }
}
