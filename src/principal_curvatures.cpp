#include "../include/principal_curvatures.h"
#include "igl/adjacency_list.h"
#include "igl/pinv.h"
#include <Eigen/Eigenvalues>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);

  std::vector<std::vector<int>> L;
  igl::adjacency_list(F, L); 

  for (int i = 0; i < V.rows(); i++) {
	// BFS on vertex v to get the two ring;
	std::vector<int> V_tr;
        for (int j = 0; j < L[i].size(); j++) {
			V_tr.push_back(L[i][j]);
			for (int k = 0; k < L[j].size(); k++) {
					V_tr.push_back(L[j][k]);
			}
		}
	
        
 	// Compute matrix P of offset positions 
	Eigen::MatrixXd P(V_tr.size(), 3);
	for (int j = 0; j < V_tr.size(); j++) {
		P.row(j) = V.row(V_tr[j]) - V.row(i);
	}	

	// Compute eigenvectors of P * P^T to get principal directions and vectors u, v , w
        Eigen::MatrixXd D = P.transpose() * P;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(D);
	Eigen::Vector3d u = eigensolver.eigenvectors().col(2);
	Eigen::Vector3d v = eigensolver.eigenvectors().col(1);
	Eigen::Vector3d w = eigensolver.eigenvectors().col(0);
  	Eigen::MatrixXd A(V_tr.size(), 5);

	// Compute the coordinates of points in P in (u, v, w) basis
	Eigen::VectorXd b(V_tr.size());
	for (int i = 0; i < V_tr.size(); i++) {
		A(i, 0) = P.row(i).dot(u);
		A(i, 1) = P.row(i).dot(v);
		b(i) = P.row(i).dot(w);
	}

	// Construct quadratic height field coefficients
	A.col(2) = (A.col(0)).cwiseProduct(A.col(0));
	A.col(3) = (A.col(1)).cwiseProduct(A.col(0));
	A.col(4) = (A.col(1)).cwiseProduct(A.col(1));
	Eigen::MatrixXd A_inv;
	igl::pinv(A, -1, A_inv); // use default tolerance
	Eigen::VectorXd a = A_inv *  b;

	// compute shape operator
	Eigen::Matrix2d first_ff;
	first_ff << 1.0 + a(0) * a(0), a(0) * a(1), a(0) * a(1), 1.0 + a(1) * a(1); 
	double div = std::sqrt(1 + a(0) * a(0) + a(1) * a(1));
	Eigen::Matrix2d second_ff;
	second_ff << 2.0 * a(2) / div, a(3) / div, a(3) / div, 2.0 * a(4) /div;
	Eigen::Matrix2d S = - second_ff * first_ff.inverse();

	// eigen decompose the shape operator to obtain curvatures and directions
	eigensolver.compute(S);
	K1(i) = eigensolver.eigenvalues()[1];
	K2(i) = eigensolver.eigenvalues()[0];
	D1.row(i) = eigensolver.eigenvectors()(1,0) * u + eigensolver.eigenvectors()(1,1) * v ;
	D2.row(i) = eigensolver.eigenvectors()(0,0) * u + eigensolver.eigenvectors()(0,1) * v; 
 } 
}
