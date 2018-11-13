#include "../include/angle_defect.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  int n = F.maxCoeff() + 1;

  Eigen::MatrixXd L, A;

  igl::squared_edge_lengths(V, F, L);
  internal_angles(L, A);

  D = Eigen::VectorXd::Ones(n);
  D = D * 2 * M_PI;

  for (int i = 0; i < F.rows(); i++){
    for (int j = 0; j < 3; j++){
      int v = F(i, j);
      D(v) = D(v) - A(i,j);
    }
  } 
}
