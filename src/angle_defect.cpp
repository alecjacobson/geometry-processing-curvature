#include "../include/angle_defect.h"
#include <igl/squared_edge_lengths.h>
#include "../include/internal_angles.h"
#include <cmath>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());
  Eigen::MatrixXd l_sqr, A;
  igl::squared_edge_lengths(V,F,l_sqr);
  internal_angles(l_sqr,A);
  for (int i=0; i<F.rows(); i++) {
    for (int j=0; j<3; j++) {
      D(F(i,j)) = 2 * M_PI - A(i,j);
    }
  }
}
