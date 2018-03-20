#include "../include/angle_defect.h"
#include <math.h>
#include "internal_angles.h"
#include <igl/massmatrix.h>
#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

  Eigen::MatrixXd EL;
  igl::squared_edge_lengths(V, F, EL);

  Eigen::MatrixXd A;
  internal_angles(EL, A);

  Eigen::VectorXd angle_sums(V.rows());
  angle_sums = Eigen::VectorXd::Zero(V.rows());

  for (int i = 0; i < F.rows(); i++)
  {
	  for (int j = 0; j < 3; j++)
	  {
		  int vertex = F(i,j);
		  angle_sums(vertex) += A(i,j);
	  }
  }

  for (int i = 0; i < V.rows(); i++)
  {
	  //std::cout << angle_sums(i) << std::endl;
	  float angle_defect = 2 * M_PI - angle_sums(i);
	  D(i) = angle_defect;
  }
  
}
