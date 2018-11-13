#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <math.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  Eigen::MatrixXd l_sqr,A;
  igl::squared_edge_lengths(V, F, l_sqr);
  internal_angles(l_sqr,A);
  D = Eigen::VectorXd::Ones(V.rows())*2*M_PI;
  for (int i=0; i<F.rows(); i++){
  	for (int j=0; j<3; j++){
  	  D(F(i,j))-=A(i,j);
	}
  }
}
