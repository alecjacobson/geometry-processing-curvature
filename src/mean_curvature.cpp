#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <iostream>

using namespace std;

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V,F,L);

	Eigen::SparseMatrix<double> invM; // inverse mass matrix
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,invM);
  for (int ii = 0; ii < invM.rows(); ii++){
  	invM.coeffRef(ii,ii) = 1.0 / invM.coeffRef(ii,ii);
  }

  Eigen::MatrixXd HMat;
  HMat.resize(V.rows(), V.cols());
  HMat = invM * L * V;

  H = Eigen::VectorXd::Zero(V.rows());
  H = HMat.cwiseProduct(HMat).cwiseSqrt().rowwise().sum(); // double check whether to divide by 2
}
