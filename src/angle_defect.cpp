#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>
#include <iostream>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());
  Eigen::MatrixXd K, L;

  igl::squared_edge_lengths(V, F, L);
  internal_angles(L, K);

   for (int i = 0; i< V.rows(); i++){
    D(i) = 2*M_PI;
   }
   
   for (int i = 0; i < F.rows(); i++) {
    for(int j= 0; j < 3; j++){
        D(F(i,j)) -= K(i, j);
      }
    }
    std::cout << D << std::endl;
}
