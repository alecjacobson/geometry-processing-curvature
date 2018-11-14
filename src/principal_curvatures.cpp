#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <set>
#include <igl/pinv.h>


using namespace Eigen;
using namespace std;

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

  // compute the adjacency matrix: prepare for the two rings
  SparseMatrix<int> A;
  igl::adjacency_matrix(F, A);

  // compute the per vertex normals for the mesh
  MatrixXd N;
  igl::per_vertex_normals(V, F, N);
    
  // iterate over all vertices and solve for the quadratic surface
  for (int i = 0; i < V.rows(); i++) {

    // construct P matrix
    MatrixXd P;
    set<int> rings;
    for (SparseMatrix<int>::InnerIterator iter1(A, i); iter1; ++iter1) {
      int cur1 = iter1.row();
      for (SparseMatrix<int>::InnerIterator iter2(A, cur1); iter2; ++iter2) {
          int cur2 = iter2.row();
          rings.insert(cur2);
      }
    }
    // discard i
    rings.erase(i);
    P.resize(rings.size(), 3);
    int j = 0;
    for (set<int>::iterator iter = rings.begin(); iter != rings.end(); ++iter, j++) {
      int cur = *iter;
      P.row(j) = V.row(cur) - V.row(i);
    }

    // PCA
    SelfAdjointEigenSolver<Matrix3d> PCA(P.transpose() * P);
    MatrixXd evs = PCA.eigenvectors();
		Vector3d u = evs.col(2);
    Vector3d v = evs.col(1);
    Vector3d w = evs.col(0);

    // uv terms
    MatrixXd UV(P.rows(), 5);
    UV.col(0) = P * u;
    UV.col(1) = P * v;
    UV.col(2) = (UV.col(0).array().square()).matrix();
    UV.col(3) = (UV.col(0).array() * UV.col(1).array()).matrix();
    UV.col(4) = (UV.col(1).array().square()).matrix();

    // solve for coefficients
    MatrixXd UV_inv;
		igl::pinv(UV, UV_inv);
		VectorXd a = UV_inv * (P * w);

    // construct shape operators
    double E = 1 + a(0) * a(0);
    double F = a(0) * a(1);
    double G = 1 + a(1) * a(1);
    double denom = sqrt(a(0) * a(0) + 1 + a(1) * a(1));
    double e = 2.0 * a(2) / denom;
    double f = a(3) * 1.0 / denom;
    double g = 2.0 * a(4) / denom;

    MatrixXd MatL(2, 2), MatU(2, 2);
    MatL << e, f, f, g;
    MatU << E, F, F, G;
    MatrixXd S = -MatL * MatU.inverse();

    SelfAdjointEigenSolver<MatrixXd> ei(S);

    // update
    double k1 = ei.eigenvalues()(1);
    double k2 = ei.eigenvalues()(0);
    K1(i) = k1;
    K2(i) = k2;

    Vector3d d1 = ei.eigenvectors()(1, 1) * u + ei.eigenvectors()(1, 0) * v;
    Vector3d d2 = ei.eigenvectors()(0, 1) * u + ei.eigenvectors()(0, 0) * v;
    D1.row(i) = d1;
    D2.row(i) = d2;

  }

}
