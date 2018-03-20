#include "../include/angle_defect.h"
#include "igl/massmatrix.h"
#include "internal_angles.h"
#include <math.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  // Defect lives on vertices
  D = Eigen::VectorXd::Zero(V.rows());

  Eigen::SparseMatrix<double> Mass;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Mass);

  // We will use the edge lengths (squared) to compute internal angles
  Eigen::MatrixXd All_lengths_sq;
  igl::squared_edge_lengths(V, F, All_lengths_sq);
  Eigen::MatrixXd All_angles;
  internal_angles(All_lengths_sq, All_angles);

  // Sum the angles around each vertex
  for(int faceNumber = 0; faceNumber < F.rows(); faceNumber++){
    for(int vertexNumber = 0; vertexNumber < F.cols(); vertexNumber++){
      // Each time we see the vertex, add the defect from its corner of the triangle
      D(F(faceNumber, vertexNumber)) += All_angles(faceNumber, vertexNumber);
    }
  }

  for(int vertexNumber = 0; vertexNumber < V.rows(); vertexNumber++){
    // Seems like the mass is already accounted for
    D(vertexNumber) = (2*M_PI - D(vertexNumber)); // / Mass.coeffRef(vertexNumber, vertexNumber)
  }
}
