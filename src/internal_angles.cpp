#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  // Add with your code
  Eigen::ArrayXd aSquare, bSquare, cSquare, cosTheta;
	A.resize(l_sqr.rows(), l_sqr.cols());
	for (int ii = 0; ii < l_sqr.cols(); ii++){
		aSquare = l_sqr.col((ii+1)%3).array();
		bSquare = l_sqr.col((ii+2)%3).array();
		cSquare = l_sqr.col((ii  )%3).array();
		cosTheta = (aSquare + bSquare - cSquare) / (2*aSquare.sqrt()*bSquare.sqrt());
		A.col(ii) = cosTheta.acos();
	}
}
