#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  int num_v = V.rows();
  D.resize(num_v);
  
  // initialize wit 2pi
  for (int i = 0; i < num_v; i++) {
    D[i] = 2 * M_PI;
  }
  
  // get edge lengths and angles
  Eigen::MatrixXd lengths(F.rows(), 3), angles;
  igl::squared_edge_lengths(V, F, lengths);
  internal_angles(lengths, angles);
  
  // subtract internal angles
  // angle order is [1,2], [2,0], [0,1] due to lengths
  double corner1, corner2, corner3;
  for (int j = 0; j < F.rows(); j++) {
    D[F(j,0)] -= angles(j,0);
    D[F(j,1)] -= angles(j,1);
    D[F(j,2)] -= angles(j,2);
  }
}
