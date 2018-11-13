#include "../include/angle_defect.h"
#include <igl/squared_edge_lengths.h>
#include <igl/adjacency_matrix.h>
#include <Eigen/Sparse>
#include "internal_angles.h"
#include <iostream>

using namespace std;

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());
  Eigen::MatrixXd l_sqr;
  igl::squared_edge_lengths(V, F, l_sqr);
  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);
  D = 2*M_PI * Eigen::VectorXd::Ones(V.rows());

  // Iterate through the faces and subtract internal angles
  // from the corresponding corners
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      int vtx = F(i, j);
      D(vtx) = D(vtx) - A(i, j);
    }
  }
}
