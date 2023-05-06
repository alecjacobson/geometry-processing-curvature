#include "../include/principal_curvatures.h"
#include <set>
#include "igl/adjacency_matrix.h"
#include "igl/per_vertex_normals.h"
#include "igl/pinv.h"

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
    Eigen::SparseMatrix<double> Ad;
    Eigen::MatrixXd normals;
    igl::adjacency_matrix(F, Ad);
    igl::per_vertex_normals(V, F, normals);
    
    for (int i = 0; i < V.rows(); i++) {
        Eigen::MatrixXd P, X, X_i;
        Eigen::VectorXd u, v, w, A;
        Eigen::Matrix2d S, S_left, S_right;
        std::set<int> neighbour, neighbour2;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Ad,i); it; ++it) {
            if (it.row() != i) {
                neighbour.insert(it.row());
                neighbour2.insert(it.row());
            }
        }
        for (int n: neighbour) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Ad,n); it; ++it) {
                if (it.row() != n) {
                    neighbour2.insert(it.row());
                }
            }
        }
        P.resize(neighbour2.size(), 3);
        int z = 0;
        for (int y: neighbour2) {
            P.row(z++) = V.row(y) - V.row(i);
        }
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(P.transpose() * P);
        w = es.eigenvectors().col(0);
        v = es.eigenvectors().col(1);
        u = es.eigenvectors().col(2);
        if (w.dot(normals.row(i)) < 0) w *= -1;
        X.resize(P.rows(), 5);
        X.col(0) = P * v;
        X.col(1) = P * u;
        X.col(2) = (X.col(0)).array().square();
        X.col(3) = (X.col(0)).array() * (X.col(1)).array();
        X.col(4) = (X.col(1)).array().square();
        igl::pinv(X, X_i);
        A = X_i * (P * w);
        
        double d = sqrt(A(0) * A(0) + 1 + A(1) * A(1));
        S_right << (1 + A(0) * A(0)), A(0) * A(1), A(0) * A(1), 1 + A(1) * A(1);
        S_left << (2 * A(2) / d), (A(3) / d), (A(3) / d), (2 * A(4) / d);
        S = -S_left * S_right.inverse();
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(S);
        K1(i) = es2.eigenvalues()(1);
        D1.row(i) = es.eigenvectors()(1, 1) * u + es.eigenvectors()(1, 0) * v;
        K2(i) = es2.eigenvalues()(0);
        D2.row(i) = es.eigenvectors()(0, 1) * u + es.eigenvectors()(0, 0) * v;
    }
}
