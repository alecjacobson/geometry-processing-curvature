#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <igl/pinv.h>
#include <iostream>
#include <igl/principal_curvature.h>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
#ifdef USE_IGL_PRINCIPAL_CURVATURE
    igl::principal_curvature(V,F,D1,D2,K1,K2);
#else
  // Replace with your code
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);


  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F,A);

  Eigen::SparseMatrix<double> A2 = A * A;


  for(int k=0; k<A2.outerSize(); ++k) {
      int size = A.innerVector(k).nonZeros();
      Eigen::MatrixXd P(size,3);
      int count = 0;
      for(typename Eigen::SparseMatrix<double>::InnerIterator it (A,k); it; ++it) {
          P.row(count++) = V.row(it.index());
      }
      P.rowwise() -= V.row(k);
      Eigen::Matrix3d Q = P.transpose() * P;
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> uv_eig(Q);


      Eigen::Matrix<double,2,3> uv = uv_eig.eigenvectors().block(0,1,3,2).transpose();
      Eigen::RowVector3d w = uv_eig.eigenvectors().col(0).transpose();
      Eigen::VectorXd B = P * w.transpose();
      Eigen::MatrixXd LC = P * uv.transpose();


      Eigen::MatrixXd data(LC.rows(),5);
      data.col(0) = LC.col(0);//u
      data.col(1) = LC.col(1);//v
      data.col(2) = LC.col(0).array().pow(2);//u^2
      data.col(3) = LC.rowwise().prod();//uv
      data.col(4) = LC.col(1).array().pow(2);//v^2
      Eigen::MatrixXd data_inv;
      igl::pinv(data,data_inv);
      Eigen::VectorXd A = data_inv * B;



      Eigen::Matrix2d f1;
      Eigen::Matrix2d f2;
      double a0_2 = std::pow<double>(A(0),2);
      double a1_2 = std::pow<double>(A(1),2);
      f1(0,0) = 1 + a0_2;
      f1(1,0) = f1(0,1) = A(0) * A(1);
      f1(1,1) = 1 + a1_2;

      double denom = std::sqrt(1 + a0_2 + a1_2);
      f2(0,0) = 2 * A(2) / denom;
      f2(1,0) = f2(0,1) = A(3) / denom;
      f2(1,1) = 2 * A(4) / denom;


      Eigen::Matrix2d S = -f2 * f1.inverse();
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> S_eig(S);

      Eigen::Matrix<double,2,3> PD = (uv.transpose() * S_eig.eigenvectors().transpose()).transpose();
      Eigen::Vector2d PC = S_eig.eigenvalues();

     
      K1(k) = PC(0);
      K2(k) = PC(1);

      D1.row(k) = PD.row(0);
      D2.row(k) = PD.row(1);

  }

#endif
}
