#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <igl/sort.h>

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

  Eigen::SparseMatrix<int> adj;
  igl::adjacency_matrix(F, adj);
  Eigen::SparseMatrix<int> adj2 = adj*adj;

  for (int i = 0; i < V.rows(); ++i)
  {
    Eigen::MatrixXd P(adj2.col(i).nonZeros(), 3);
    int p = 0;
    for (Eigen::SparseMatrix<int>::InnerIterator it(adj2, i); it; ++it)
    {
      P.row(p) = V.row(it.index()) - V.row(i);
      p++;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd W = svd.matrixV();
    Eigen::MatrixXd UV = P*W.leftCols(2);
    Eigen::MatrixXd B = P*W.rightCols(1);

    Eigen::MatrixXd C(B.rows(), 5);
    for (int j = 0; j < UV.rows(); ++j)
        C.row(j) << UV(j, 0), UV(j, 1), UV(j, 0)*UV(j, 0), UV(j, 0)*UV(j, 1), UV(j, 1)*UV(j, 1);

    igl::pinv(C, C);
    Eigen::VectorXd A = C*B;

    double E = 1 + A(0)*A(0);
    double _F = A(0)*A(1);
    double G = 1 + A(1)*A(1);
    double e = 2 * A(2) / sqrt(A(0)*A(0) + 1 + A(1)*A(1));
    double f = A(3) / sqrt(A(0)*A(0) + 1 + A(1)*A(1));
    double g = 2 * A(4) / sqrt(A(0)*A(0) + 1 + A(1)*A(1));

    Eigen::MatrixXd SL(2, 2);
    Eigen::MatrixXd SR(2, 2);
    SL << e, f, f, g;
    SR << E, _F, _F, G;
    Eigen::MatrixXd S = -SL*SR.inverse();

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(S);
    Eigen::MatrixXd eigen_vecs = eigensolver.eigenvectors().real();
    Eigen::VectorXd eigen_vals = eigensolver.eigenvalues().real().cwiseAbs();
    Eigen::VectorXd eigen_vals_sorted;
    Eigen::VectorXd eigen_index_sorted;
    igl::sort(eigen_vals, 1, false, eigen_vals_sorted, eigen_index_sorted);

    Eigen::VectorXd d1 = eigen_vecs.col(eigen_index_sorted(0));
    Eigen::VectorXd d2 = eigen_vecs.col(eigen_index_sorted(1));

    K1(i) = eigen_vals_sorted(0);
    K2(i) = eigen_vals_sorted(1);

    D1.row(i) = W.leftCols(2)*d1;
    D2.row(i) = W.leftCols(2)*d2;
  }
}
