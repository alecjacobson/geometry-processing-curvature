#ifndef MEAN_CURVATURE_H
#define MEAN_CURVATURE_H
#include <Eigen/Core>
// Compute the discrete mean curvature at each vertex of a mesh (`V`,`F`) by
// taking the signed magnitude of the mean curvature normal as a _pointwise_ or
// _integral average_ quantity.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh face indices into V
// Outputs:
//   H  #V list of mean curvature values in units 1/m²
//
void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H);

#endif
