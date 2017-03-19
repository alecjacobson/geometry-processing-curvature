#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <igl/per_vertex_normals.h>

using namespace Eigen;

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
	SparseMatrix<int> adj;
	igl::adjacency_matrix(F, adj);
	
	MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	int n = V.rows(), f = F.rows();

	D1.resize(n, 3);
	D2.resize(n, 3);
	K1.resize(n);
	K2.resize(n);
	MatrixXd P(n, 3);

	for (int v = 0; v < n; ++v)
	{
		int k = 0;

		VectorXi adj_1 = VectorXi::Zero(n);	//1-ring neighbourhhod adjacency
		adj_1(v) = 1;
		std::vector <int> neib_1;	//1-ring neighbourhood \ {v}
		std::vector <int> neib_2;	//2-ring neighbourhood \ {v' \in neib_1}
		for (SparseMatrix<int>::InnerIterator it(adj, v); it; ++it)
		{
			adj_1(it.row()) = 1;
			neib_1.push_back(it.row());
		}
		
		std::for_each(neib_1.begin(), neib_1.end(), [&] (int v1) {
			for (SparseMatrix<int>::InnerIterator it(adj, v1); it; ++it)
				if (!adj_1(it.row()))
					neib_2.push_back(it.row());

			P.row(k++) = V.row(v1) - V.row(v);
		});

		std::for_each(neib_2.begin(), neib_2.end(), [&](int v2) {
			P.row(k++) = V.row(v2) - V.row(v);
		});

		// PCA
		JacobiSVD <MatrixXd> svd(P.topRows(k), ComputeThinV);

		MatrixXd m = P.topRows(k) * svd.matrixV();
		MatrixXd coeffMat(k, 5);
		VectorXd b = m.col(2);

		// check if b is "height" or "depth", and if latter, reverse it
		if (svd.matrixV().col(2).dot(N.row(v)) < 0)
			b *= -1;

		for (int m_it = 0; m_it < m.rows(); ++m_it)
		{
			coeffMat(m_it, 0) = m(m_it, 0);
			coeffMat(m_it, 1) = m(m_it, 1);
			coeffMat(m_it, 2) = m(m_it, 0) * m(m_it, 0);
			coeffMat(m_it, 3) = m(m_it, 0) * m(m_it, 1);
			coeffMat(m_it, 4) = m(m_it, 1) * m(m_it, 1);
		}

		MatrixXd coeffMatPinv;
		igl::pinv(coeffMat, coeffMatPinv);

		VectorXd A = coeffMatPinv*b;

		auto E = 1 + A(0)*A(0);
		auto F = A(0)*A(1);
		auto G = 1 + A(1)*A(1);
		auto denom = std::sqrt(1 + A(0)*A(0) + A(1)*A(1));
		auto e = 2 * A(2) / denom;
		auto f = A(3) / denom;
		auto g = 2 * A(4) / denom;

		Matrix2d efg;
		efg(0, 0) = e;
		efg(0, 1) = efg(1, 0) = f;
		efg(1, 1) = g;

		Matrix2d EFG;
		EFG(0, 0) = E;
		EFG(0, 1) = EFG(1, 0) = F;
		EFG(1, 1) = G;

		Matrix2d S = -efg * EFG.inverse();

		SelfAdjointEigenSolver<Matrix2d> eig(S);
		Vector2d eigVal = eig.eigenvalues();
		Matrix2d eigVec = eig.eigenvectors();

		//sort eigenpairs (descending)
		if (eigVal(1) > eigVal(0))
		{
			eigVal.reverseInPlace();
			eigVec.rowwise().reverse().eval();
		}

		K1(v) = eigVal(0);
		K2(v) = eigVal(1);

		D1.row(v) = (svd.matrixV().leftCols(2)*eigVec.col(0)).transpose();
		D2.row(v) = (svd.matrixV().leftCols(2)*eigVec.col(1)).transpose();
	}
}
