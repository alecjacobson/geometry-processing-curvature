#include "../include/angle_defect.h"
#include "igl/squared_edge_lengths.h"
#include "internal_angles.h"
#include "igl/PI.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
    Eigen::MatrixXd l_sqr, X;
    igl::squared_edge_lengths(V, F, l_sqr);
    internal_angles(l_sqr, X);
    
    D = 2 * igl::PI * Eigen::VectorXd::Ones(V.rows());
    for (int i = 0; i < F.rows(); i++) {
        D(F(i, 0)) += -X(i, 0);
        D(F(i, 1)) += -X(i, 1);
        D(F(i, 2)) += -X(i, 2);
    }
}
