#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Constant(V.rows(), 2 * M_PI);

  Eigen::MatrixXd l_sqr(F.rows(),F.cols());
  igl::squared_edge_lengths(V,F,l_sqr);

  Eigen::MatrixXd A;
  internal_angles(l_sqr,A);
  for(int i = 0; i < F.rows(); ++i) {
      auto&& f = F.row(i);
      auto&& a = A.row(i);
      for(int j = 0; j < 3; ++j) {
          D(f(j)) -= a(j);
      }
  }

}
