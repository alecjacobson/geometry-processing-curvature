#include "../include/principal_curvatures.h"
#include <iostream>
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/pinv.h>
#include <igl/sort.h>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
	D1 = Eigen::MatrixXd(V.rows(), 3);
	D2 = Eigen::MatrixXd(V.rows(), 3);
	K1 = Eigen::VectorXd(V.rows());
	K2 = Eigen::VectorXd(V.rows());
	
	//Get the one and two-hop adjacent vertices. Adj_two has all the information that Adj_one has
	//because we have a symetric matrix.
	//https://www.coursera.org/learn/advanced-data-structures/lecture/38zHi/support-multiplying-adjacency-matrices	
	Eigen::SparseMatrix<int> Adj_one, Adj_two;
	igl::adjacency_matrix(F, Adj_one);		
	Adj_two = Adj_one * Adj_one; //Very cool property of adjacency matrices, get two hop neighbors

	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	for (int i = 0; i < V.rows(); i++) {

		//Step 1: Find all 2-hop neighbors and append them to the Point matrix P, with relative positions.
		//Essentially the index of non-zero entries along one col/row in the adjacency matrix are the
		//correct points.
		Eigen::MatrixXd P(Adj_two.col(i).nonZeros(), 3);
		int points = 0;
		for (typename Eigen::SparseMatrix<int>::InnerIterator it(Adj_two, i); it; ++it) {
			P.row(points) = V.row(it.index()) - V.row(i);
			points++;
		}		

		///Step 2: Perform PCA to get the 2 principle components - Perform SVD and the V in
		//USV* are the principle components
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinV); 		
		Eigen::MatrixXd M = svd.matrixV();
		Eigen::MatrixXd uv = M.leftCols(2); //First 2 components uv
		Eigen::VectorXd w = M.col(2);		//3rd component w
		//Ensure consistency in direction of normal and the height component
		if (w.dot(N.row(i)) < 0) {
			w *= -1;
		}
		//Compute the S and B matrix
		Eigen::MatrixXd S = P * uv;
		Eigen::VectorXd B = P * w;

		///Step 3: Construct and solve the least squares problem to fit a quadratic to the points
		Eigen::MatrixXd coeff(P.rows(), 5);
		for (int j = 0; j < S.rows(); j++) {
			coeff(j, 0) = S(j, 0);				//u
			coeff(j, 1) = S(j, 1);				//v
			coeff(j, 2) = S(j, 0) * S(j, 0);	//u*u
			coeff(j, 3) = S(j, 0) * S(j, 1);	//u*v
			coeff(j, 4) = S(j, 1) * S(j, 1);	//v*v
		}
		Eigen::MatrixXd invCoeff;
		igl::pinv(coeff, invCoeff);		
		Eigen::VectorXd A = invCoeff * B; //Stores the coefficients of the fit quadratic

		///Step 4: Form and compute the shape operator
		double E, F, G, e, f, g, bottom;
		bottom = sqrt((1 + A(0)*A(0) + A(1)*A(1)));
		E = 1 + A(0)*A(0);
		F = A(0)*A(1);
		G = 1 + A(1)*A(1);
		e = 2 * A(2) / bottom;
		f = A(3) / bottom;
		g = 2 * A(4) / bottom;

		Eigen::Matrix2d shape1, shape2, shape;
		shape1 << e, f, f, g;
		shape2 << E, F, F, G;
		shape = -shape1 * shape2.inverse();

		///Step 5: Eigen Value decomposistion to get the principle curvatures and the tangent directions
		//Eigen::EigenSolver<Eigen::Matrix2d> eig;
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig;
		eig.compute(shape);
		Eigen::Vector2d eigValues = eig.eigenvalues();
		Eigen::Matrix2d eigVectors = eig.eigenvectors();
		
		Eigen::Vector2d vals, ind;
		igl::sort(eigValues.real(), 1, false, vals, ind);
		K1[i] = vals(0);
		K2[i] = vals(1);
		Eigen::VectorXd d1, d2;
		d1 = eigVectors.col(ind(0));
		d2 = eigVectors.col(ind(1));
		D1.row(i) = M.leftCols(2)*d1;
		D2.row(i) = M.leftCols(2)*d2;

	}

}
