#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <cmath>
#include <iostream>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());
  Eigen::MatrixXd L, A, DA;
  igl::squared_edge_lengths(V, F, L);
  internal_angles(L, A);

  Eigen::MatrixXd angle_sum = Eigen::MatrixXd::Zero(V.rows(), 1);
  for (int i = 0; i < F.rows(); i++){
  	for (int j = 0; j < F.cols(); j++){
  	  int index = F(i, j);
  	  angle_sum(index, 0) += A(i, j);
  	}
  }
  D = (2 * M_PI - angle_sum.array());

}
