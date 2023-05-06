#include "../include/principal_curvatures.h"
#include <cmath>
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <Eigen/SVD>
#include <igl/pinv.h>
#include <Eigen/Eigenvalues>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  // // Replace with your code
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);

  Eigen::SparseMatrix<int> A;
  igl::adjacency_matrix(F,A);

  // two rings
  for (int i=0; i<V.rows(); i++) {
    std::vector<Eigen::VectorXd> two_rings;
    two_rings.push_back(V.row(i)-V.row(i));

    for (Eigen::SparseMatrix<int>::InnerIterator iter(A,i);iter;++iter) {
      int pos = iter.row();
      two_rings.push_back(V.row(pos) - V.row(i));
      for (Eigen::SparseMatrix<int>::InnerIterator iters(A,pos);iters;++iters) {
        two_rings.push_back(V.row(iters.row()) - V.row(i));
      }
    }
    
    Eigen::MatrixXd P;
    P.resize(two_rings.size(),3);
    for (int j=0; j<two_rings.size(); j++) {
      P.row(j) = two_rings[j];
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P.transpose()*P, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd uvw = svd.matrixU() * svd.singularValues().asDiagonal();
    Eigen::MatrixXd uv;
    uv.resize(two_rings.size(),5);
    uv.col(0) = P * uvw.col(0);
    uv.col(1) = P * uvw.col(1);
    uv.col(2) = uv.col(0).array() * uv.col(0).array();
    uv.col(3) = uv.col(0).array() * uv.col(1).array();
    uv.col(4) = uv.col(1).array() * uv.col(1).array();

    Eigen::MatrixXd pw = P * uvw.col(2);

    Eigen::MatrixXd uv_1;
    igl::pinv(uv,uv_1);
    Eigen::VectorXd a = uv_1 * pw;
    double E,F,G,e,f,g;
    Eigen::MatrixXd L,U;
    L.resize(2,2);
    U.resize(2,2);
    E = 1 + a[0] * a[0];
    F = a[0] * a[1];
    G = 1 + a[1] * a[1];
    L << E,F,
         F,G;

    double deno = sqrt(1 +  a[0] * a[0] + a[1] * a[1]);
    e = 2 * a[2] / deno;
    f = a[3] / deno;
    g = 2 * a[4] / deno;
    U << e,f,
         f,g;

    Eigen::MatrixXd S = - U * L.inverse();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(S);
    int k1,k2;
    if (eig.eigenvalues()(0) > eig.eigenvalues()(1)) {
      k1 = 0;
      k2 = 1;
      K1(i) = eig.eigenvalues()(0);
      K2(i) = eig.eigenvalues()(1);
    }
    else {
      k2 = 0;
      k1 = 1;
      K2(i) = eig.eigenvalues()(0);
      K1(i) = eig.eigenvalues()(1);
    }

    D1.row(i) = eig.eigenvectors()(1, 1) * uvw.col(0) + eig.eigenvectors()(1, 0) * uvw.col(1);
    D2.row(i) = eig.eigenvectors()(0, 1) * uvw.col(0) + eig.eigenvectors()(0, 0) * uvw.col(1);
  }


}
