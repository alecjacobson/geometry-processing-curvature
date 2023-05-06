#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>

static Eigen::MatrixXd constructMatrixP (
    const Eigen::MatrixXd & V,
    const Eigen::SparseMatrix<double> A,
    const int v)
{
  std::vector<int> twoRingVertIndices;
  twoRingVertIndices.push_back(v);

  /* Build 1st ring vertices */
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
      if (it.row() == v) {
        twoRingVertIndices.push_back(it.col());
      }
    }
  }

  const int ring1Size = twoRingVertIndices.size();

  /* Build 2nd ring vertices from 1st ring*/
  for (int j = 0; j < ring1Size; j++) {
    for (int k = 0; k < A.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
        if (it.row() == twoRingVertIndices[j]) {
          if (std::find(twoRingVertIndices.begin(), twoRingVertIndices.end(), it.col()) == twoRingVertIndices.end()) {
            twoRingVertIndices.push_back(it.col());
          }
        }
      }
    }
  }

  Eigen::MatrixXd P(twoRingVertIndices.size() - 1, 3);

  for (int k = 1; k < twoRingVertIndices.size(); k++) {
    P.row(k - 1) = V.row(twoRingVertIndices[k]) - V.row(v);
  }

  return P;
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
  K1.resize(V.rows());
  K2.resize(V.rows());
  D1.resize(V.rows(), 3);
  D2.resize(V.rows(), 3);

  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F, A);

  for (int v = 0; v < V.rows(); v++) {
    Eigen::MatrixXd P = constructMatrixP(V, A, v);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P.transpose()*P, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXd S = (P*svd.matrixV()).leftCols(2);
    Eigen::MatrixXd B = (P*svd.matrixV()).rightCols(1);

    Eigen::MatrixXd H(P.rows(),5);
    H.leftCols(2) = S;
    H.col(2) = S.col(0).cwiseProduct(S.col(0));
    H.col(3) = S.col(0).cwiseProduct(S.col(1));
    H.col(4) = S.col(1).cwiseProduct(S.col(1));

    Eigen::MatrixXd Hinv;
    igl::pinv(H,Hinv);

    Eigen::VectorXd a = Hinv*B;

    double e,f,g,E,F,G;
    E = 1 + a(0)*a(0);
    F = a(0)*a(1);
    G = 1 + a(1)*a(1);
    e = 2*a(2)/sqrt(a(0)*a(0)+1+a(1)*a(1));
    f =   a(3)/sqrt(a(0)*a(0)+1+a(1)*a(1));
    g = 2*a(4)/sqrt(a(0)*a(0)+1+a(1)*a(1));

    Eigen::Matrix2d fundamental1;
    fundamental1 << E, F,
                    F, G;
    Eigen::Matrix2d fundamental2;
    fundamental2 << e, f,
                    f, g;

    Eigen::Matrix2d Shape = -fundamental2 * fundamental1.inverse();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(Shape, Eigen::ComputeFullU | Eigen::ComputeFullV);

    K1(v) = svd2.singularValues()(0);
    K2(v) = svd2.singularValues()(1);
    D1.row(v) = svd.matrixU().leftCols(2) * svd2.matrixU().col(0);
    D2.row(v) = svd.matrixU().leftCols(2) * svd2.matrixU().col(1);
  }
}
