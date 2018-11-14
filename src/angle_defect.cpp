#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <cmath>

void angle_defect(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    Eigen::VectorXd &D)
{
  //D = Eigen::VectorXd::Zero(V.rows());
  Eigen::MatrixXd L, A;
  igl::squared_edge_lengths(V, F, L);
  internal_angles(L, A);

  //2pi - A
  D.resize(V.rows());
  for (int i = 0; i < F.rows(); i++)
  {
    D(F(i, 0)) = 2.00 * M_PI - A(i, 0);
    D(F(i, 1)) = 2.00 * M_PI - A(i, 1);
    D(F(i, 2)) = 2.00 * M_PI - A(i, 2);
  }
}
