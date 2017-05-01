#include "../include/mean_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/per_vertex_normals.h"

using namespace Eigen;

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{		
	//Compute the information required for the discrete laplacian 
	Eigen::SparseMatrix<double> L, M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> M_inv(M.rows(), M.cols());
	for (int i = 0; i < M.rows();i++) {
		M_inv.coeffRef(i, i) = 1.0f / M.coeff(i, i);		
	}

	//Apply the discrete laplacian to the vertices, should use the 
	//positive definite version of the cotmatrix to keep conventions consistent.
	Eigen::MatrixXd Hn = M_inv * (-L) * V;

	//Calculate the magnitude of each row of Hn, but ensure sign is correct
	//by comparing with the per vertex normal 
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	H.resize(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		if (Hn.row(i).dot(N.row(i)) >0) {
			H(i) = Hn.row(i).norm();
		}
		else {
			H(i) = -Hn.row(i).norm();
		}
	}

}
