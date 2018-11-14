#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <math.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  Eigen::MatrixXd l_sqr;
  igl::squared_edge_lengths(V, F, l_sqr);
  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);
  
  D = 2*M_PI*Eigen::VectorXd::Ones(V.rows());
  for (int fi = 0; fi < F.rows(); fi++) {
    for (int vi = 0; vi < 3; vi++) {
      D[F(fi, vi)] -= A(fi, vi);
    }
  }
}
