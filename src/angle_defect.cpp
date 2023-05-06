#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <igl/PI.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Constant(V.rows(), 2*igl::PI);

  Eigen::MatrixXd L;
  igl::squared_edge_lengths(V,F,L);

  Eigen::MatrixXd A;
  internal_angles(L, A);

  for(int f = 0; f < F.rows(); f++){
    D(F(f,0)-1) -= A(f,0);
    D(F(f,1)-1) -= A(f,1);
    D(F(f,2)-1) -= A(f,2);
  }
}
