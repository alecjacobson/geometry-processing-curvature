#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <igl/slice.h>
#include <igl/sort.h>
#include <igl/per_vertex_normals.h>
#include <unordered_set>
#include <Eigen/Eigenvalues>


#include <igl/principal_curvature.h>

void principal_curvatures(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd & D1,
	Eigen::MatrixXd & D2,
	Eigen::VectorXd & K1,
	Eigen::VectorXd & K2)
{
  // Replace with your code
	K1 = Eigen::VectorXd::Zero(V.rows());
	K2 = Eigen::VectorXd::Zero(V.rows());
	D1 = Eigen::MatrixXd::Zero(V.rows(), 3);
	D2 = Eigen::MatrixXd::Zero(V.rows(), 3);

	Eigen::SparseMatrix<double> adjacency;
	igl::adjacency_matrix(F, adjacency);

	Eigen::MatrixXd vertex_normals;
	igl::per_vertex_normals(V, F, vertex_normals);

	//For each vertex, need to grab a 2-ring of points and put them in a matrix.
	//Then compute the principal components, do a least squares fit, then another eigen decomposition.
	int num_verts = V.rows();

	for (int v = 0; v < num_verts; v++) {
		//To get the two ring, kind of a pain.
		//First, get the non-zero entries of adjacency.row(v). This is the 1-ring.
		//Then, for each of those rows, grab all entries. But dont want duplicates
		Eigen::VectorXi two_ring_verts;
		{
			std::unordered_set<int> two_ring;
			for (int i = 0; i < adjacency.cols(); i++) {
				if (adjacency.coeff(v, i) != 0) {
					two_ring.emplace(i);
					for (int j = 0; j < adjacency.cols(); j++) {
						if (adjacency.coeff(i, j) != 0) {
							two_ring.emplace(j);
						}
					}
				}
			}
			two_ring_verts.resize(two_ring.size());
			auto it = two_ring.cbegin();
			for (auto i = 0; i < two_ring.size(); i++, it++) {
				two_ring_verts(i) = *it;
			}
		}

		Eigen::MatrixXd P;
		igl::slice(V, two_ring_verts, 1, P);
		P.rowwise() -= V.row(v);

		//Now do principle component analysis on P. 
		Eigen::RowVector3d mean = P.colwise().mean();

		P.rowwise() -= mean;

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(P.transpose()*P);

		Eigen::Vector3d sorted_eigenvalues;
		Eigen::Vector3i sorted_indices;
		igl::sort(eigen.eigenvalues().cwiseAbs(),1, false, sorted_eigenvalues,sorted_indices);

		
		Eigen::Matrix<double, 3, 2> rotation;
		igl::slice(eigen.eigenvectors(), sorted_indices.head(2), 2, rotation);
		Eigen::MatrixXd S(P.rows(), 2);
		Eigen::VectorXd B(P.rows());

		Eigen::Matrix<double, 3, 1> height = eigen.eigenvectors().col(sorted_indices(2));
		Eigen::Vector3d n = vertex_normals.row(v).transpose();
		if (n.dot(height) < 0) {
			height *= -1;
		}

		for (int i = 0; i < P.rows(); i++) {
			S.row(i) = P.row(i) * rotation;
			B.row(i) = P.row(i) * height;
		}
		

		//Do least-squares fitting on these variables.
		Eigen::Matrix<double, 5, 1> A;
		{
			Eigen::MatrixXd quadratic_system(P.rows(), 5);
			for (int i = 0; i < P.rows(); i++) {
				quadratic_system(i, 0) = S(i, 0);
				quadratic_system(i, 1) = S(i, 1);
				quadratic_system(i, 2) = S(i, 0)*S(i, 0);
				quadratic_system(i, 3) = S(i, 0)*S(i, 1);
				quadratic_system(i, 4) = S(i, 1)*S(i, 1);
			}
			Eigen::MatrixXd inv;
			igl::pinv(quadratic_system, inv);
			A = inv*B;
		}

		//Finally, need to assemble the shape operator, then decompose.
		Eigen::Matrix2d shape_operator;
		{
			Eigen::Matrix2d second, first; //Second and first fundamental forms
			double denominator = std::sqrt(A(0)*A(0) + 1 + A(1)*A(1));
			second(0, 0) = 2 * A(2) / denominator;
			second(0, 1) = second(1, 0) = A(3) / denominator;
			second(1, 1) = 2 * A(4) / denominator;

			first(0, 0) = 1 + A(0)*A(0);
			first(0, 1) = first(1, 0) = A(0)*A(1);
			first(1, 1) = 1 + A(1)*A(1);

			shape_operator = -second*first.inverse();
		}

		{
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(shape_operator);
			double k1 = eigensolver.eigenvalues()(0);
			double k2 = eigensolver.eigenvalues()(1);
			auto eigenvecs = eigensolver.eigenvectors();
			//Now need to convert the eigenvectors back to the cartesian space.
			Eigen::Vector3d d1 = rotation*eigenvecs.col(0).real() + mean.transpose();
			Eigen::Vector3d d2 = rotation*eigenvecs.col(1).real() + mean.transpose();
			
			
			if (k1 < k2) {
				std::swap(k1, k2);
				std::swap(d1, d2);
			}

			K1(v) = k1;
			K2(v) = k2;
			D1.row(v) = d1;
			D2.row(v) = d2;
			
		}
	}

}
