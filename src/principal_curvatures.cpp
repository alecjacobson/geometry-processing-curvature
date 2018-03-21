#include "../include/principal_curvatures.h"
#include <vector>
#include "igl/adjacency_matrix.h"
#include "igl/pinv.h"
#include "igl/per_vertex_normals.h"

using namespace Eigen;

void principal_curvatures(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		Eigen::MatrixXd & D1,
		Eigen::MatrixXd & D2,
		Eigen::VectorXd & K1,
		Eigen::VectorXd & K2) {
	int vertexCount = V.rows();

	K1 = Eigen::VectorXd::Zero(vertexCount);
	K2 = Eigen::VectorXd::Zero(vertexCount);
	D1 = Eigen::MatrixXd::Zero(vertexCount, 3);
	D2 = Eigen::MatrixXd::Zero(vertexCount, 3);

	// We can verify the sign of the curvature by comparing against the
	// normals.
	Eigen::MatrixXd normals;
	igl::per_vertex_normals(V, F, normals);

	// The squared ajacency matrix will tell us which nodes are within a
	// distance of two edges.
	Eigen::SparseMatrix<int> adjacency;
	igl::adjacency_matrix(F, adjacency);
	adjacency = adjacency * adjacency;

	// Now we can iterate through each vertex and solve for the quadratic
	// surface.
	for (int i = 0; i < vertexCount; i++) {
		Vector3d vertex = V.row(i);

		// Using the iterator style from the Eigen documentation, we can
		// iterate over all the non-0 entries of each row of the adjacency
		// matrix to locate all neighbours of the vertex (and the vertex
		// itself).
		SparseMatrix<int> neighbourAdjacency = adjacency.row(i);
		std::vector<Vector3d> neighbours;
		for (int k = 0; k < neighbourAdjacency.outerSize(); ++k) {
			for (SparseMatrix<int>::InnerIterator it(neighbourAdjacency, k); 
					it; ++it) {

				int neighbourIndex = it.col();
				Vector3d neighbourVertex = V.row(neighbourIndex);
				neighbours.push_back(neighbourVertex - vertex);
			}
		}

		// Construct the matrix P from the readme.
		int neighbourCount = neighbours.size();
		MatrixXd P(neighbourCount, 3);
		for (int j = 0; j < neighbourCount; j++) {
			P.row(j) = neighbours[j];
		}

		// Perform PCA on P. Eigen-values are in increasing order.
		SelfAdjointEigenSolver<MatrixXd> pcaP(P.transpose() * P);
		const MatrixXd& eigenVectorsP = pcaP.eigenvectors();
		Vector3d u = eigenVectorsP.col(2);
		Vector3d v = eigenVectorsP.col(1);
		Vector3d w = eigenVectorsP.col(0);

		// We want to make sure that w faces along the normal direction.
		if (w.dot(normals.row(i)) < 0) {
			w = -w;
		}

		// Transform the vertices into the space defined by the PCA
		// components.
		VectorXd pU = P * u;
		VectorXd pV = P * v;
		VectorXd pW = P * w;

		// The system ends up being D * a = w where D is a matrix made with
		// the quadratic pU and pV terms, a is the vector of coefficients 
		// and w is the vector of vertex pW values.
		MatrixXd D(P.rows(), 5);
		D.col(0) = pU;
		D.col(1) = pV;
		D.col(2) = pU.cwiseProduct(pU);
		D.col(3) = pU.cwiseProduct(pV);
		D.col(4) = pV.cwiseProduct(pV);

		// We can use the pseudo-inverse to express the solution as matrix
		// multiplication by the pseudo-inverse.
		MatrixXd psuedoInverse;
		igl::pinv(D, psuedoInverse);
		VectorXd a = psuedoInverse * pW;

		// Now its time to construct the shape operator.
		double e = 2 * a(2) / sqrt(a(0) * a(0) + 1 + a(1) * a(1));
		double f = a(3) / sqrt(a(0) * a(0) + 1 + a(1) * a(1));
		double g = 2 * a(4) / sqrt(a(0) * a(0) + 1 + a(1) * a(1));
		double E = 1 + a(0) * a(0);
		double F = a(0) * a(1);
		double G = 1 + a(1) * a(1);

		MatrixXd secondForm(2, 2);
		secondForm(0, 0) = e;
		secondForm(1, 0) = f;
		secondForm(0, 1) = f;
		secondForm(1, 1) = g;

		MatrixXd firstForm(2, 2);
		firstForm(0, 0) = E;
		firstForm(1, 0) = F;
		firstForm(0, 1) = F;
		firstForm(1, 1) = G;
		
		MatrixXd S = -1 * secondForm * firstForm.inverse();

		// Now we have to do PCA again to extract information from the
		// shape operator.
		SelfAdjointEigenSolver<MatrixXd> pcaS(S);

		const VectorXd& eigenValuesS = pcaS.eigenvalues();
		K1(i) = eigenValuesS(0);
		K2(i) = eigenValuesS(1);

		// We have to transform the eigen vectors back into world 
		// coordinates.
		const MatrixXd& eigenVectorsS = pcaS.eigenvectors();
		D1.row(i) = eigenVectorsS(0, 0) * u + eigenVectorsS(1, 0) * v;
		D2.row(i) = eigenVectorsS(0, 1) * u + eigenVectorsS(1, 1) * v;
	}
}
