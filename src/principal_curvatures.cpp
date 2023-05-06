#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/pinv.h>
#include <math.h>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  // Replace with your code
	K1.resize(V.rows());
	K2.resize(V.rows());
	D1.resize(V.rows(), 3);
	D2.resize(V.rows(), 3);

	// Compute normal for each vertex for deciding height in "w-direction"
	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);
	
	// Compute adjacency matrix for sampling points in the vicinity
	Eigen::SparseMatrix<int> A1;
	igl::adjacency_matrix(F, A1);
	

	// Compute adjacency matrix A^2 which identifies points that are two-hop away from each other
	Eigen::SparseMatrix<int> A2 = A1 * A1;

	// Compute adjacency matrix A for "two-ring" of each vertex
	Eigen::SparseMatrix<int> A = A1 + A2;

	// Loop over each vertex of a given mesh
	for (int i = 0; i < V.rows(); i++) {
		// Calculate matrix P which contains positions of those "two-ring" adjacent vertices for 
		// vertex i relative to vertex i
		int k = A.innerVector(i).nonZeros();
		Eigen::MatrixXd P(k, 3);
		int j = 0;
		for (Eigen::SparseMatrix<int>::InnerIterator it(A, i); it; ++it) {
			P.row(j) = V.row(it.row()) - V.row(i);
			j++;
		}
		
		// Do principal component analysis on P to get S (two most principal directions) and B 
		// (least principal direction) for each corresponding point in P
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> PCA(P.transpose()*P);
		Eigen::Vector3d u, v, w;
		u = PCA.eigenvectors().col(2);
		v = PCA.eigenvectors().col(1);
		w = PCA.eigenvectors().col(0);
		Eigen::Vector3d n = N.row(i);
		if (w.dot(n) < 0)
			w = -w;
		Eigen::MatrixXd S(k, 2);
		S.col(0) = P * u;
		S.col(1) = P * v;
		Eigen::VectorXd B = P * w;
		
		// Solve for 5 unknown coefficients (a1, a2 ...) in the following equation by least-square fitting
		// w = a1u + a2v + a3u^2 + a4uv + a5v^2
		// Construct uvMat in form uvMat.row(i) = [u_i, v_i, u_i*u_i, u_i*v_i, v_i*v_i]
		Eigen::MatrixXd uvMat(k, 5);
		uvMat.col(0) = S.col(0);
		uvMat.col(1) = S.col(1);
		for (int j = 0; j < k; j++) {
			uvMat(j, 2) = S(j, 0)*S(j, 0);
			uvMat(j, 3) = S(j, 0)*S(j, 1);
			uvMat(j, 4) = S(j, 1)*S(j, 1);
		}
		// Solve unknown coefficients by coeff = uvMat^-1*B
		Eigen::MatrixXd uvMat_inv;
		igl::pinv(uvMat, uvMat_inv);
		Eigen::VectorXd coeff = uvMat_inv * B;

		// Construct shape operator SO = −[ef][EF]-1
		//                                [fg][FG]
		double a1sq = coeff(0)*coeff(0);
		double a1a2 = coeff(0)*coeff(1);
		double a2sq = coeff(1)*coeff(1);
		double E = 1 + a1sq;
		double F = a1a2;
		double G = 1 + a2sq;
		double denominator = sqrt(a1sq + 1 + a2sq);
		double e = 2 * coeff(2) / denominator;
		double f = coeff(3) / denominator;
		double g = 2 * coeff(4) / denominator;
		Eigen::Matrix2d small, big;
		small << e, f, f, g;
		big << E, F, F, G;
		Eigen::MatrixXd SO = -small*(big.inverse());
		
		// Decomposition of S to get principal curvatures and principal tangent directions
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> dcp(SO);
		Eigen::VectorXd eigenvalues = dcp.eigenvalues();
		Eigen::MatrixXd eigenvectors = dcp.eigenvectors();
		K1(i) = eigenvalues[0];
		K2(i) = eigenvalues[1];
		D1.row(i) = eigenvectors(0, 0)*u + eigenvectors(1, 0)*v;
		D2.row(i) = eigenvectors(0, 1)*u + eigenvectors(1, 1)*v;
	}
}
