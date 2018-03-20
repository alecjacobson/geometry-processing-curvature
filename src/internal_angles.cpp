#include "../include/internal_angles.h"
#include <math.h>

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  int numFaces = l_sqr.rows();
  // Three internal angles per trinagular face
  A.resize(numFaces, 3);

  for(int faceNumber = 0; faceNumber < numFaces; faceNumber++){
    for(int startingVertex = 0; startingVertex < 3; startingVertex++){
      int v1 = startingVertex;
      int v2 = (startingVertex + 1) % 3;
      int v3 = (startingVertex + 2) % 3;

      // Apply the law of cosines
      double c_squared = l_sqr.coeffRef(faceNumber, v1);
      double others_squared = l_sqr.coeffRef(faceNumber, v2) + l_sqr.coeffRef(faceNumber, v3);
      double numerator = c_squared - others_squared;
      double ab = sqrt(l_sqr.coeffRef(faceNumber, v2) * l_sqr.coeffRef(faceNumber, v3));
      double denominator = - 2.0 * ab;

      double ratio = numerator / denominator;

      A.coeffRef(faceNumber, v1) = acos(ratio);

    }
  }
}
