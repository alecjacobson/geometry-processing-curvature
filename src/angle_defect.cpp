#include "../include/angle_defect.h"
#include "igl/squared_edge_lengths.h"
#include "internal_angles.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Constant(V.rows(), 2 * M_PI);

  Eigen::MatrixXd l_sqr;
  igl::squared_edge_lengths(V, F, l_sqr);

  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);

  int num_F = A.rows();
  for(int f_idx = 0; f_idx < num_F; f_idx++){
    // indices into V
    int v_i = F(f_idx, 0);
    int v_j = F(f_idx, 1);
    int v_k = F(f_idx, 2);

    D(v_i) -= A(f_idx, 0);
    D(v_j) -= A(f_idx, 1);
    D(v_k) -= A(f_idx, 2);
  }
}
