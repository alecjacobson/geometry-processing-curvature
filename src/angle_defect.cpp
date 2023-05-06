#include "../include/angle_defect.h"
#include "igl/squared_edge_lengths.h"
#include "internal_angles.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  // Get list of internal angles
  Eigen::MatrixXd L;
  igl::squared_edge_lengths(V, F, L);
  Eigen::MatrixXd A;
  internal_angles(L, A);

  // From angles compute angle defects
  D = Eigen::VectorXd::Ones(V.rows()) * 2 * std::acos(-1); // acos(-1) = pi
  for (int i = 0; i < F.rows(); i++ ){
	for (int j = 0; j <= 2; j++) {
 	D(F(i,j)) = D(F(i,j)) - A(i,j);
	} 
  }
}
