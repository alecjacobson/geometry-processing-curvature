#include "../include/angle_defect.h"
#include <Eigen/Sparse>
#include <igl/squared_edge_lengths.h>
#include "../include/internal_angles.h"
void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  Eigen::MatrixXd l_sqr;
  igl::squared_edge_lengths(V, F, l_sqr);
  
  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);

  D.resize(V.rows());
  D.setConstant(2*M_PI);
  
  for (int i = 0; i < F.rows(); i++)
    for (int j = 0; j < 3; j++)
        D(F(i, j)) = -A(i, j);
}
