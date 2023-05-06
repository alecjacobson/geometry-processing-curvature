#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  // H = Eigen::VectorXd::Zero(V.rows());


	// Get cotangent laplacian
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	// Get mass matrix
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

	// Compute the inverse of the mass matrix. Since we're using the diagonalized version,
	// all we have to do is compute reciprocals along the diagonal. Doing this the same way
	// as in the previous assignment.
	Eigen::SparseMatrix<double> M_inv(M.cols(), M.cols());
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	for (int ii = 0; ii < M.cols(); ii++)
		triplets.push_back(T(ii, ii, 1.0/M.coeff(ii, ii)));

	M_inv.setFromTriplets(triplets.begin(), triplets.end());

	// Now the mean curvature matrix is just the product M_inv*L*V
	Eigen::MatrixXd H_temp = M_inv*L*V;
	H.resize(V.rows());

	// Now compute the norms, and use dot products with each normal vector to fix the signs, as per readme.
	// First, extract normals on each face.

	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	for (int ii = 0; ii < H_temp.rows(); ii++)
	{
		// Dot product of this H vector with this normal
		double this_dot = N(ii, 0)*H_temp(ii, 0) + N(ii, 1)*H_temp(ii, 1) + N(ii, 2)*H_temp(ii, 2);

		// If this is negative, use the norm of this H row with a negative sign out front.
		// Otherwise, use the norm unaltered.
		double this_norm = H_temp.row(ii).norm();
		if (this_dot < 0.0)
			H(ii) = -this_norm;
		else
			H(ii) = this_norm;
	}

	return;

}
