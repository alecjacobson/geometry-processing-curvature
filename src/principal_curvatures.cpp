#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/slice.h>
#include <algorithm>
#include <limits>
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
  D1 = Eigen::MatrixXd::Zero(V.rows(), 3);
  D2 = Eigen::MatrixXd::Zero(V.rows(), 3);
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());

  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F, A);
  for (int i = 0; i < A.cols(); i++) {
    std::vector<int> two_rings_list;
    for (int j = 0; j < A.cols(); j++) {
      if (A.coeff(i, j) == 1) {
        two_rings_list.push_back(j);
        for (int k = 0; k < A.cols(); k++) {
          if (A.coeff(j, k) == 1 && k != i) {
            if (std::find(two_rings_list.begin(), two_rings_list.end(), k) != two_rings_list.end()) {
              two_rings_list.push_back(k);
            }
          }
        }
      }
    }
    Eigen::VectorXi two_rings = Eigen::VectorXi::Zero(two_rings_list.size());
    for (int t = 0; t < two_rings.rows(); t++) {
      two_rings(t) = two_rings_list.at(t);
    }
    Eigen::MatrixXd points = Eigen::MatrixXd::Zero(two_rings.rows(), 3);
    igl::slice(V, two_rings, 1, points);
    points.rowwise() -= V.row(i);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(points.transpose() * points);
    Eigen::VectorXd eval = eigensolver.eigenvalues();
    Eigen::MatrixXd evec = eigensolver.eigenvectors();
    int idx1, idx2, idx3;
    double val1 = std::numeric_limits<double>::lowest();
    double val2 = std::numeric_limits<double>::lowest();
    double val3 = std::numeric_limits<double>::max();
    for (int e = 0; e < eval.rows(); e++) {
      if (eval(e) > val1) {
        val2 = val1;
        idx2 = idx1;
        val1 = eval(e);
        idx1 = e;
      }
      else if(eval(e) > val2) {
        val2 = eval(e);
        idx2 = e;
      }
      if (eval(e) < val3) {
        val3 = eval(e);
        idx3 = e;
      }
    }
    Eigen::MatrixXd top2 = Eigen::MatrixXd::Zero(3, 2);
    Eigen::MatrixXd bot = Eigen::MatrixXd::Zero(3, 1);
    top2.col(0) = evec.col(idx1);
    top2.col(1) = evec.col(idx2);
    bot.col(0) = evec.col(idx3);

    Eigen::MatrixXd S = points * top2;
    Eigen::MatrixXd B = points * bot;

    Eigen::MatrixXd coeff(points.rows(), 5);
    coeff.col(0) = S.col(0);
    coeff.col(1) = S.col(1);
    coeff.col(2) = S.col(0).array() * S.col(0).array();
    coeff.col(3) = S.col(0).array() * S.col(1).array();
    coeff.col(4) = S.col(1).array() * S.col(1).array();
    Eigen::MatrixXd pinv;
    igl::pinv(coeff, pinv);
    Eigen::VectorXd a = pinv * B;

    double e, f, g, x, y, z;
    e = 2*a(2)/sqrt(a(0) * a(0) + 1 + a(1) * a(1));
    f = 2*a(3)/sqrt(a(0) * a(0) + 1 + a(1) * a(1));
    g = 2*a(4)/sqrt(a(0) * a(0) + 1 + a(1) * a(1));
    x = 1 + a(0) * a(0);
    y = a(0) * a(1);
    z = 1 + a(1) * a(1);

    Eigen::Matrix2d left, right, Shape;
    left << e, f, f, g;
    right << x, y, y, z;
    Shape = - left * right.inverse();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver2D(Shape);
    evec = eigensolver2D.eigenvectors();
    K1(i) = eigensolver2D.eigenvalues()(0);
    K2(i) = eigensolver2D.eigenvalues()(1);
    D1.row(i) = top2 * evec.col(0);
    D2.row(i) = top2 * evec.col(1);
  }
}
