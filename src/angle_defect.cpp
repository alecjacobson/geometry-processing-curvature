#include "../include/angle_defect.h"

#include <igl/squared_edge_lengths.h>
#include "internal_angles.h"
#include <Eigen/Dense>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());

  // calculate edge length
  Eigen::MatrixXd Flen;
  igl::squared_edge_lengths(V, F, Flen);

  Eigen::MatrixXd VAngles;
  internal_angles(Flen, VAngles);

  D = Eigen::VectorXd::Ones(V.rows()) * 3.1416;
  for (int i = 0; i< F.rows(); i++)
  {
	  for(int j = 0; j< 3; j++)
	  {
		  D(F(i, j)) -= VAngles(F(i, j));
	  }
  }


}
