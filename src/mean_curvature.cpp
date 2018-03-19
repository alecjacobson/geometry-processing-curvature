#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <iostream>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
	std::cout << "test" << std::endl;
	Eigen::SparseMatrix<double> L; 
	igl::cotmatrix(V, F, L);
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	H = (M.cwiseInverse() * L * V).rowwise().squaredNorm();
	std::cout << "test2" << std::endl;
}
