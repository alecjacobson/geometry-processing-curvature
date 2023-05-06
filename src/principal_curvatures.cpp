#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
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
  Eigen::SparseMatrix<int> A;
  igl::adjacency_matrix(F, A);
  A = A * A;

  for(int i=0; i<V.rows(); i++){
    Eigen::SparseVector<int> Ring = A.col(i);
    Eigen::MatrixXd P(Ring.nonZeros(), 3);

    int k = 0;
    for (Eigen::SparseVector<double>::InnerIterator it(Ring); it; ++it)
        P.row(k++) = V.row(it.index()) - V.row(i);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P.transpose() * P, Eigen::ComputeFullU);
    Eigen::MatrixXd SD = svd.matrixU();
    std::cout << SD.cols() << std::endl;

    Eigen::VectorXd t_1 = P * SD.col(0);
    Eigen::VectorXd t_2 = P * SD.col(1);
    Eigen::VectorXd t_3 = t_1.array().colwise()*t_1.array();
    Eigen::VectorXd t_4 = t_1.array().colwise()*t_2.array();
    Eigen::VectorXd t_5 = t_2.array().colwise()*t_2.array();
    Eigen::MatrixXd t(P.rows(), 5);
    t<< t_1, t_2, t_3, t_4, t_5;

    Eigen::MatrixXd a;
    igl::pinv(t, a);
    a = a * P * SD.col(2);

    double tsqrt = sqrt(a(0) * a(0) + 1 + a(1) * a(1));
    double E = 1 + a(0) * a(0);
    double F1 = a(0) * a(1);
    double G = 1 + a(1) * a(1);
    double e = 2 * a(2) / tsqrt ;
    double f = a(3) / tsqrt;
    double g = 2 * a(4) / tsqrt;

    Eigen::MatrixXd left(2, 2);
    left << e, f,
          f, g;
    Eigen::MatrixXd right(2, 2);
    right << E, F1, 
          F1, G;

    Eigen::MatrixXd S = -1*left*right.inverse();
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;
    ces.compute(S);

    K1(i) = ces.eigenvalues()[0].real();
    K2(i) = ces.eigenvalues()[1].real();
    D1.row(i) = SD.leftCols(2) * ces.eigenvectors().col(0).real();
    D2.row(i) = SD.leftCols(2) * ces.eigenvectors().col(1).real();
  }
}
