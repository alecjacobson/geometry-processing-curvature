#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/doublearea.h>
#include <igl/squared_edge_lengths.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  // omega_i = (2pi - sum(angles on vertex))
  Eigen::MatrixXd double_areas, edge_lengths, angles;
  igl::doublearea(V, F, double_areas);
  igl::squared_edge_lengths(V, F, edge_lengths);
  internal_angles(edge_lengths, angles);

  // compose sums of incident angles and areas for all vertices
  Eigen::VectorXd incident_angle_sums(V.rows());
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      incident_angle_sums(F(i, j)) += angles(i, j);
    }
  }

  // now fill D!
  D = Eigen::VectorXd::Zero(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    D(i) = (2.0 * M_PI - incident_angle_sums(i));
  }
}
