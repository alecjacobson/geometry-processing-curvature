#include "../include/principal_curvatures.h"
#include "igl/adjacency_matrix.h"
#include "igl/pinv.h"
#include "igl/slice.h"
#include "igl/sort.h"
#include "igl/per_vertex_normals.h"
#include <Eigen/Eigenvalues> 
#include "igl/principal_curvature.h"

using namespace Eigen;
using namespace igl;

void get2ringMatrix(const MatrixXd& V, const SparseMatrix<int>& A, const SparseMatrix<int>& A2, MatrixXd& P, int vIndex)
{
	// entry in A(i,j) are non zero if a vertex can be reached with two steps (2 ring).
	// one ring should be automatically included since its an undirected graph

	int ring2Size = A2.col(vIndex).nonZeros();
	P.resize(ring2Size, 3);

	int j = 0;

	for (SparseMatrix<int>::InnerIterator it(A2, vIndex); it; ++it)
	{
		auto curV = V.row(it.row());
		auto refV = V.row(vIndex);
		P.row(j++) = curV - refV;
	}
}
	

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
	D1 = Eigen::MatrixXd::Zero(V.rows(),3);
	D2 = Eigen::MatrixXd::Zero(V.rows(),3);

	//igl::principal_curvature(V, F, D1, D2, K1, K2); return; //For comparing results
	
	SparseMatrix<int> A;
	adjacency_matrix(F, A);
	
	SparseMatrix<int> A2 = A * A; //two ring adjacency. 
	
	for (int i = 0; i < A.outerSize(); ++i)
	{
		MatrixXd P; 
		get2ringMatrix(V, A, A2, P, i);

		JacobiSVD<MatrixXd> svd(P.transpose()*P, ComputeThinU);
		//auto lambdas = svd.singularValues();
		MatrixXd lambdaVec = svd.matrixU();
		MatrixXd uv = lambdaVec.leftCols(2);
		MatrixXd S = P * uv;

		MatrixXd w = lambdaVec.rightCols(1);
		MatrixXd B = P * w;
		
		MatrixXd UV(B.size(), 5);

		for (int i = 0; i < S.rows(); ++i) 
		{
			double u = S(i, 0), v = S(i, 1);
			UV(i, 0) = u;
			UV(i, 1) = v;
			UV(i, 2) = u*u;
			UV(i, 3) = u*v;
			UV(i, 4) = v*v;
		}

		MatrixXd X;
		pinv(UV, X);

		VectorXd coef = X*B;
				
		double p = 1.0/sqrt(coef(0)*coef(0) + 1 + coef(1)*coef(1));
		double e = 2 * coef(2)*p;
		double f = coef(3) * p;
		double g = 2 * coef(4) * p;
		double E = 1 + coef(0)*coef(0);
		double F = coef(0)*coef(1);
		double G = 1 + coef(1)*coef(1);

		MatrixXd M1(2, 2), M2(2,2);
		M1 << e, f, f, g;
		M2 << E, F, F, G;
		
		MatrixXd M = -1 * M1 * M2.inverse();
		
		EigenSolver<MatrixXd> eSolver;
		eSolver.compute(M);
		VectorXd eSorted, eSortedIndex;
		sort(eSolver.eigenvalues().real(), 1, false, eSorted, eSortedIndex);
		
		double k1 = eSorted(0);
		double k2 = eSorted(1);
		VectorXd t1 = eSolver.eigenvectors().col(eSortedIndex(0)).real();
		VectorXd t2 = eSolver.eigenvectors().col(eSortedIndex(1)).real();

		D1.row(i) = lambdaVec.leftCols(2)*t1;
		D2.row(i) = lambdaVec.leftCols(2)*t2;
		K1(i) = k1;	
		K2(i) = k2;
	}
}


