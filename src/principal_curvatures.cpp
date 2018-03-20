#include "../include/principal_curvatures.h"
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <vector>
#include <cmath>
#include <iostream>
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

  Eigen::SparseMatrix<int> adjacency;
  igl::adjacency_matrix(F, adjacency);
  for (int i = 0; i < V.rows() ; i ++){
    Eigen::VectorXd v_pos = V.row(i);
    //construct two ring
    std::vector<Eigen::VectorXd> list_of_p;
    list_of_p.push_back(v_pos - v_pos);
    //iteration through all adjcent vertices first-ring
    for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency,i); it; ++it){
      auto first_index = it.row();
      Eigen::VectorXd v1_pos = V.row(first_index);
      list_of_p.push_back(v1_pos - v_pos);
      //second-ring
      for (Eigen::SparseMatrix<int>::InnerIterator it2(adjacency,first_index); it2; ++it2){
        auto second_index = it2.row();
        Eigen::VectorXd v2_pos = V.row(second_index);
          list_of_p.push_back(v2_pos - v_pos);
      }
    }
    //copy over to a matrix
    Eigen::MatrixXd P(list_of_p.size(),3);
    for (int j = 0; j < list_of_p.size(); j++){
      P.row(j) = list_of_p[j];
    }
    //do svd on P
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU|Eigen::ComputeThinV);
    auto score = svd.matrixU() * svd.singularValues().asDiagonal();

    Eigen::MatrixXd u = score.col(0);
    Eigen::MatrixXd v = score.col(1);
    Eigen::MatrixXd w = score.col(2);

    Eigen::MatrixXd to_solve(list_of_p.size(), 5);
    to_solve.col(0) = u;
    to_solve.col(1) = v;
    to_solve.col(2) = u.array() * u.array();
    to_solve.col(3) = u.array() * v.array();
    to_solve.col(4) = v.array() * v.array();

    Eigen::MatrixXd pinv;
    igl::pinv(to_solve, pinv);
    //solve for coeff
    Eigen::VectorXd coeff = pinv * w;

    //Construct shape operator
    double E1 = 1 + coeff(0) * coeff(0);
    double F1 = coeff(0) * coeff(1);
    double G1 = 1 + coeff(1) * coeff(1);
    Eigen::Matrix2d right;
    right << E1 , F1 ,
            F1 , G1;
    double temp = sqrt(1 + coeff(0) * coeff(0) + coeff(1) * coeff(1));
    double e = 2 * coeff(2) / temp;
    double f =     coeff(3) / temp;
    double g = 2 * coeff(4) / temp;
    Eigen::Matrix2d left;
    left << e, f,
            f, g;

    auto S = - left * right.inverse();

    //use self adjoint eigen solver since its allr real number
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(S);
    //suppose k1 = 0, k2 =1
    int k1 = 0;
    int k2 = 1;
    if (eig.eigenvalues()(0) < eig.eigenvalues()(1)){ //pick the largest one to be k1
      k1 = 1;
      k2 = 0;
    }
    K1(i) = eig.eigenvalues()(k1);
    K2(i) = eig.eigenvalues()(k2);

    //transfer back to world
    Eigen::MatrixXd base = svd.matrixV().leftCols(2);
    D1.row(i) = eig.eigenvectors().row(k1) * base.transpose();
    D2.row(i) = eig.eigenvectors().row(k2) * base.transpose();

  }

}
