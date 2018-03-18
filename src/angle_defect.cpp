#include "../include/angle_defect.h"
#include "../include/internal_angles.h"

#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Ones(V.rows()) * 2.0 * M_PI;

  Eigen::MatrixXd l_sqr;
  igl::squared_edge_lengths(V, F, l_sqr);

  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);

  for(int i = 0; i < F.rows(); i++) {
    for(int j = 0; j < 3; j++) {
      int v_i = F(i, j);
      D[v_i] -= A(i, j);
    }
  }
}
