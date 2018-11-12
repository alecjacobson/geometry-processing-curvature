#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Constant(V.rows(), 2*M_PI);

  // Compute internal angles
  Eigen::MatrixXd l_sqr;
  igl::squared_edge_lengths(V, F, l_sqr);
  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);

  // Loop through all faces
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < F.cols(); ++j) {

        D(F(i, j)) -= A(i, j);

    } // end loop j
  } // end loop i
}
