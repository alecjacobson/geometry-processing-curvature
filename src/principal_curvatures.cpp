#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <iostream>
#include <igl/sort.h>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
	int n_V = V.rows();

	D1 = Eigen::MatrixXd(n_V, 3);
	D2 = Eigen::MatrixXd(n_V, 3);
	K1 = Eigen::VectorXd(n_V);
	K2 = Eigen::VectorXd(n_V);

	Eigen::SparseMatrix<int> Adj;
	igl::adjacency_matrix(F, Adj);

	// Note: This question has made me doubt my habit of declaring all my loop scope variables outside the loop
	// I can't even remember why I started doing it. To save computation cycles in redclaring variables every iteration?
	Eigen::SparseMatrix<int> Adj2 = Adj*Adj; // entries contain number of walks of length 2 from i -> j
	Eigen::MatrixXd P;
	int num_adj;
	int num_adj2;
	int l;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd;
	Eigen::MatrixXd W;
	Eigen::MatrixXd S;
	Eigen::VectorXd B;
	Eigen::VectorXd LS;
	Eigen::MatrixXd T;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double E;
	double J;
	double G;
	double e;
	double f;
	double g;
	Eigen::MatrixXd Sh1(2, 2);
	Eigen::MatrixXd Sh2(2, 2);
	Eigen::MatrixXd Sh;
	Eigen::VectorXd d1;
	Eigen::VectorXd d2;
	double k1;
	double k2;
	Eigen::MatrixXd T_pinv;
	Eigen::EigenSolver<Eigen::MatrixXd> es;
	Eigen::VectorXd e_values;
	Eigen::VectorXd e_sort_index;
	for (int i = 0; i < n_V; ++i) 
	{
		num_adj = Adj.col(i).nonZeros();
		num_adj2 = Adj2.col(i).nonZeros() - 1; // don't count yourself
		P = Eigen::MatrixXd(num_adj + num_adj2, 3);
		
		l = 0;
		for (Eigen::SparseMatrix<int>::InnerIterator it(Adj, i); it; ++it) // loop over rows for column i
		{
			P.row(l) = V.row(it.row()) - V.row(i);
			l++;
		}

		for (Eigen::SparseMatrix<int>::InnerIterator it(Adj2, i); it; ++it) // loop over rows for column i
		{
			if (it.row() == i) continue;

			P.row(l) = V.row(it.row()) - V.row(i);
			l++;
		}

		svd.compute(P, Eigen::ComputeThinU | Eigen::ComputeThinV); 
		svd.computeV();
		W = svd.matrixV();

		S = P*W.leftCols(2);
		B = P*W.col(2);

		T = Eigen::MatrixXd(P.rows(), 5);
		for (int j = 0; j < S.rows(); ++j) 
		{
			T(j, 0) = S(j, 0);
			T(j, 1) = S(j, 1); 
			T(j, 2) = S(j, 0)*S(j,0);
			T(j, 3) = S(j, 0)*S(j, 1);
			T(j, 4) = S(j, 1)*S(j, 1);
		}

		igl::pinv(T, T_pinv);

		LS = T_pinv*B;

		a1 = LS(0);
		a2 = LS(1);
		a3 = LS(2);
		a4 = LS(3);
		a5 = LS(4);

		E = 1 + a1*a1;
		J = a1*a2;
		G = 1 + a2*a2;
		e = 2 * a3 / sqrt(a1*a1 + 1 + a2*a2);
		f = a4 / sqrt(a1*a1 + 1 + a2*a2);
		g = 2 * a5 / sqrt(a1*a1 + 1 + a2*a2);

		Sh1 << e, f, f, g;
		Sh2 << E, J, J, G;
		Sh = -Sh1*Sh2.inverse();

		es.compute(Sh);
		igl::sort(es.eigenvalues().real(), 1, false, e_values, e_sort_index);
		k1 = e_values(0);
		k2 = e_values(1);
		d1 = es.eigenvectors().col(e_sort_index(0)).real();
		d2 = es.eigenvectors().col(e_sort_index(1)).real();

		K1(i) = k1;
		K2(i) = k2; 

		D1.row(i) = W.leftCols(2)*d1;
		D2.row(i) = W.leftCols(2)*d2;
	}
}
