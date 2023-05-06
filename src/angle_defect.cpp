#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#define _USE_MATH_DEFINES
#include <cmath>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());

  // make the length sqared matrix 
  Eigen::MatrixXd l_sqr(F.rows(), F.cols());

  // iterate over faces
  for (int i = 0; i < F.rows(); i ++) {
    l_sqr(i, 0) = (V.row(F(i, 1)) - V.row(F(i, 2))).squaredNorm();
    l_sqr(i, 1) = (V.row(F(i, 0)) - V.row(F(i, 2))).squaredNorm();
    l_sqr(i, 2) = (V.row(F(i, 1)) - V.row(F(i, 0))).squaredNorm();
  }

  // get the internal angle matrix
  Eigen::MatrixXd A;
  internal_angles(l_sqr, A);

  // instead of iterating over all the faces per each vertex i, 
  // let's just iterate over faces and keep adding as we encounter vertices of that face.

  D.resize(V.rows());
  D = Eigen::VectorXd::Ones(V.rows()) * 2 * M_PI;

  // iterate over faces
  for (int i = 0; i < F.rows(); i++) {

    // update D for each vertex of the face
    // decrement from the total 2*pi
    D(F(i, 0)) -= A(i, 0);
    D(F(i, 1)) -= A(i, 1);
    D(F(i, 2)) -= A(i, 2);
  }

}
