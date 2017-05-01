#include "../include/angle_defect.h"
#include "internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <igl/massmatrix.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  /*
    the angle defect is measured as the difference between all of the internal angles at a point and 2π
   */

  D = Eigen::VectorXd::Constant(V.rows(), 2*M_PI);
  
  // squared edge lengths (used to get the internal angles)
  Eigen::MatrixXd l_sqr;
  igl::squared_edge_lengths(V, F, l_sqr);

  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);

  // loop over all of the faces and add the contribution of the internal angle to the corresponding point
  for(int32_t i = 0; i < F.rows(); i++)
  {
    Eigen::RowVectorXi tri = F.row(i);
    Eigen::RowVectorXd iAngle = A.row(i);

    for(int32_t j = 0; j < 3; j++)
    {
      D(tri(j)) -= iAngle(j);
    }
  }
}
