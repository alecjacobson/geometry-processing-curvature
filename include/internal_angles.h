#ifndef INTERNAL_ANGLES_H
#define INTERNAL_ANGLES_H
#include <Eigen/Core>
// Given (squared) edge-lengths of a triangle mesh `l_sqr` compute the internal
// angles at each corner (a.k.a. wedge) of the mesh.
//
// Inputs:
//    l_sqr  #F by 3 list squared edge-lengths opposite respective corner
// Outputs:
//    A  #F by 3 list of internal angles incident on respective corner
//
void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A);
#endif

