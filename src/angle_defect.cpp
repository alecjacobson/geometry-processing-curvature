#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <math.h>

using namespace Eigen;

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{

  // get pi
  double pi = acos(-1);

  D = VectorXd::Ones(V.rows()) * pi * 2;

  Eigen::MatrixXd l_sqr, A;
  igl::squared_edge_lengths(V, F, l_sqr);
  internal_angles(l_sqr, A); // get internal angles
  
  for (int i = 0; i < F.rows(); i++) {
    for (int j= 0; j < 3; j++) {
      D(F(i,j)) = D(F(i,j)) - A(i, j);
    }
  }

}
