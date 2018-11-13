#include "../include/principal_curvatures.h"

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  int n = F.maxCoeff() + 1;

  // Replace with your code
  K1 = Eigen::VectorXd::Zero(n);
  K2 = Eigen::VectorXd::Zero(n);
  D1 = Eigen::MatrixXd::Zero(n,3);
  D2 = Eigen::MatrixXd::Zero(n,3);

  Eigen::SparseMatrix<int> A, AA;
  igl::adjacency_matrix(F, A);

  for (int i = 0; i < n; i++){
    A.coeffRef(i,i) = 1;
    // add diagonal matrix
  }

  AA = A * A; // Contains everything within 2 hops

  for (int i = 0; i < n; i++){
    // Get two ring
    Eigen::MatrixXd P;
    P.resize(AA.col(i).nonZeros(), 3);
    
    int idx = 0;
    for (Eigen::SparseMatrix<int>::InnerIterator it(AA, i); it; ++it)
    {
      P.row(idx) = V.row(it.index()) - V.row(i);
      idx++;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd u = svd.matrixV().col(0);
    Eigen::VectorXd v = svd.matrixV().col(1);
    Eigen::VectorXd w = svd.matrixV().col(2);

    Eigen::VectorXd Pu = P * u, Pv = P * v, Pw = P * w;
    // coordinates along eigenvectors

    Eigen::MatrixXd W;
    W.resize(P.rows(), 5);

    W.col(0) = Pu;
    W.col(1) = Pv;
    W.col(2) = Pu.cwiseProduct(Pu);
    W.col(3) = Pu.cwiseProduct(Pv);
    W.col(4) = Pv.cwiseProduct(Pv);

    Eigen::MatrixXd W_pinv;
    igl::pinv(W, W_pinv);
    Eigen::VectorXd A = W_pinv * Pw;
    // Readme says w, but I think it should be 
    // the height, which is Pw (or B)

    // Nasty equations
    double E = 1 + A(0) * A(0);
    double F = A(0) * A(1);
    double G = 1 + A(1) * A(1);
    double denom = sqrt(E + G - 1);
    double e = (2 * A(2)) / denom;
    double f = A(3) / denom;
    double g = (2 * A(4)) / denom;

    Eigen::Matrix2d Sleft, Sright;
    Sright << E, F, F, G;
    Sleft << e, f, f, g;
    Eigen::Matrix2d S = -Sleft * Sright.inverse();

    Eigen::JacobiSVD<Eigen::MatrixXd> svdS(S, Eigen::ComputeFullU);

    K1(i) = svdS.singularValues()(0);
    K2(i) = svdS.singularValues()(1);
    D1.row(i) = svd.matrixV().leftCols(2) * svdS.matrixU().col(0);
    D2.row(i) = svd.matrixV().leftCols(2) * svdS.matrixU().col(1);
  }
}
