#include "../include/principal_curvatures.h"
#include <Eigen/Dense>
#include <vector>
#include <igl/per_vertex_normals.h>
#include <igl/slice.h>
#include <igl/pinv.h>

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

  Eigen::MatrixXd P;
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  Eigen::VectorXi C(3);
  C(0) = 0;
  C(1) = 1;
  C(2) = 2;
  
  for (int i = 0; i < V.rows(); i++)
  {
	  std::vector<int> r1;
	  std::vector<int> r2;
	  
	  // one-ring
	  for (int j = 0; j < F.rows(); j++)
	  {
		  auto face = F.row(j);
		  if ((face.array() == i).any())
		  {
			  r1.push_back(face(0));
			  r1.push_back(face(1));
			  r1.push_back(face(2));
		  }
	  }
	  // two-ring
	  for (int j = 0; j < F.rows(); j++)
	  {
		  auto face = F.row(j);
		  for (int k = 0; k < r1.size(); k++)
		  {
			  if ((face.array() == r1[k]).any())
			  {
				  r2.push_back(face(0));
				  r2.push_back(face(1));
				  r2.push_back(face(2));
			  }
		  }
	  }

	  r1.reserve(r1.size() + r2.size());
	  r1.insert(r1.end(), r2.begin(), r2.end());
	  
	  std::sort(r1.begin(), r1.end());
	  r1.erase(std::remove(r1.begin(), r1.end(), i), r1.end());
	  r1.erase(unique(r1.begin(), r1.end()), r1.end());
	  
	  int* arr = &r1[0];
	  Eigen::Map<Eigen::VectorXi> R(arr, r1.size(), 1);
	  
	  igl::slice(V, R, C, P);
	  
	  P = P.rowwise() - V.row(i);
	  
	  Eigen::Matrix3d Pt_P = P.transpose() * P;
	  
	  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(Pt_P);
	  Eigen::MatrixXd p = eigensolver.eigenvectors();
	  //std::cout << eigensolver.eigenvalues() << std::endl;
	  
	  Eigen::VectorXd u = p.col(2);
	  Eigen::VectorXd v = p.col(1);
	  Eigen::VectorXd w = p.col(0);
	  
	  Eigen::MatrixXd H(P.rows(),5);
	  H.col(0) = P * u;
	  H.col(1) = P * v;
	  H.col(2) = H.col(0).cwiseProduct(H.col(0));
	  H.col(3) = H.col(0).cwiseProduct(H.col(1));
	  H.col(4) = H.col(1).cwiseProduct(H.col(1));
	  
	  Eigen::MatrixXd X;
	  igl::pinv(H, X);

	  float dot = w.dot(N.row(i));
	  signbit(dot) ? w = -w : w = w;

	  Eigen::VectorXd a = X * P * w;

	  double e = (2 * a(2)) / (sqrt(pow(a(0), 2) + 1 + pow(a(1), 2)));
	  double f = a(3) / (sqrt(pow(a(0), 2) + 1 + pow(a(1), 2)));
	  double g = (2 * a(4)) / (sqrt(pow(a(0), 2) + 1 + pow(a(1), 2)));
	  double E = 1 + pow(a(0), 2);
	  double F = a(0) * a(1);
	  double G = 1 + pow(a(1), 2);

	  Eigen::Matrix2d s1;
	  s1 << e, f,
			f, g;
	  Eigen::Matrix2d s2;
	  s2 << E, F,
			F, G;
	  Eigen::Matrix2d S = -s1 * s2.inverse();

	  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver_shape(S);
	  Eigen::MatrixXd s = eigensolver_shape.eigenvectors();

	  K1(i) = eigensolver_shape.eigenvalues()(1);
	  K2(i) = eigensolver_shape.eigenvalues()(0);

	  D1.row(i) = s(1, 1) * u + s(1, 0) * v;
	  D1.row(i) = s(0, 1) * u + s(0, 0) * v;
	  
  }

}
