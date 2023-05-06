#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <set>
#include <Eigen/SVD>
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
	// get adjacency matrix
	Eigen::SparseMatrix<double> Adj;
	igl::adjacency_matrix(F,Adj);

	for (int i = 0; i < V.rows(); i++) {
		std::set<int> s1, s2;
		// first ring of neighbours
		for (int j = 0; j < Adj.row(i).size(); j++) {
			if (Adj.coeff(i, j) != 0) {
				s1.insert(j);
				s2.insert(j);
			}
		}
		// second ring of neighbours
		std::set<int>::iterator it;
		for (it = s1.begin(); it != s1.end(); ++it) {
			for (int j = 0; j < Adj.row(*it).size(); j++) {
				if (Adj.coeff(*it, j) != 0 && i != j) {
					s2.insert(j);
				}
			}
		}
		// construct P
		Eigen::MatrixXd P(s2.size(), 3);
		int idx;
		for (it = s2.begin(), idx = 0; it != s2.end(); ++it, idx++) {
			P.row(idx) = V.row(*it) - V.row(i);
		}

		// principal-component analysis
		Eigen::Vector3d singular, u, v, w;
		Eigen::Matrix3d eigens;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(P.transpose() * P, Eigen::ComputeThinU);
		singular = svd.singularValues();
		eigens = svd.matrixU();
		Eigen::Vector3d::Index max_inx, min_inx;
		singular.maxCoeff(&max_inx);
		singular.minCoeff(&min_inx);
		// maximum, most principal direction
		u = eigens.col(max_inx);
		// minimum, least principal direction
		w = eigens.col(min_inx);
		// the other one
		v = eigens.col(3 - max_inx - min_inx);
		// construct S and B
		// project each direction in P onto uvw plane
		Eigen::MatrixXd S(P.rows(), 2);
		Eigen::VectorXd B(P.rows());
		for (int j = 0; j < P.rows(); j++) {
			S(j, 0) = P.row(j).dot(u/u.norm());
			S(j, 1) = P.row(j).dot(v/v.norm());
			B(j) = P.row(j).dot(w.normalized());
		}

		// compute a's
		Eigen::MatrixXd UV(P.rows(), 5), UV_I;
		for (int j = 0; j < P.rows(); j++) {
			UV(j, 0) = S(j, 0);
			UV(j, 1) = S(j, 1);
			UV(j, 2) = S(j, 0) * S(j, 0);
			UV(j, 3) = S(j, 0) * S(j, 1);
			UV(j, 4) = S(j, 1) * S(j, 1);
		}
		igl::pinv(UV, UV_I);
		Eigen::VectorXd A = UV_I * B;

		// compute shape operator
		Eigen::Matrix2d Shape, m1, m2;
		double E = 1 + A(0) * A(0);
		double F = A(0) * A(1);
		double G = 1 + A(1) * A(1);
		double e = 2 * A(2) / sqrt(1 + A(0) * A(0) + A(1) * A(1));
		double f = A(3) / sqrt(1 + A(0) * A(0) + A(1) * A(1));
		double g = 2 * A(4) / sqrt(1 + A(0) * A(0) + A(1) * A(1));
		m1 << e, f, f, g;
		m2 << E, F, F, G;
		Shape = (-1.0) * m1 * m2.inverse();

		// decompose shape operator
		Eigen::Vector2d eigenVal;
		Eigen::Matrix2d eigenVec;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd2(Shape, Eigen::ComputeThinU);
		eigenVal = svd2.singularValues();
		eigenVec = svd2.matrixU();
		
		// roll back to xyz
		Eigen::MatrixXd uv(2, 3);
		uv.row(0) = u;
		uv.row(1) = v;
		Eigen::Vector2d::Index max_idx, min_idx;
		eigenVal.maxCoeff(&max_idx);
		eigenVal.minCoeff(&min_idx);
		K1.resize(V.rows());
		K2.resize(V.rows());
		D1.resize(V.rows(),3);
		D2.resize(V.rows(),3);
		// maximum, k1
		K1(i) = eigenVal(max_idx);
		D1.row(i) = eigenVec.col(max_idx).transpose() * uv;
		// minimum, k2
		K2(i) = eigenVal(min_idx);
		D2.row(i) = eigenVec.col(min_idx).transpose() * uv;
	}
}
