#include "../include/principal_curvatures.h"
#include "igl/pinv.h"
#include "igl/adjacency_matrix.h"
#include "igl/per_vertex_normals.h"
#include <Eigen/Eigenvalues>
#include <set>

void principal_curvatures(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &D1,
        Eigen::MatrixXd &D2,
        Eigen::VectorXd &K1,
        Eigen::VectorXd &K2) {

    int N = V.rows();
    K1 = Eigen::VectorXd(N);
    K2 = Eigen::VectorXd(N);
    D1 = Eigen::MatrixXd(N, 3);
    D2 = Eigen::MatrixXd(N, 3);

    //we need vertices that share an edge.
    Eigen::SparseMatrix<int> A;
    igl::adjacency_matrix(F, A);



    for (int i = 0; i < N; i++) {
        // Step1 let's get the two_rings of v
        std::set<int> two_rings;
        for (Eigen::SparseMatrix<int>::InnerIterator iterator(A, i); iterator; ++iterator) {
            int current_row = iterator.row();
            two_rings.insert(current_row);
            for (Eigen::SparseMatrix<int>::InnerIterator innerIterator(A,
                                                                          current_row); innerIterator; ++innerIterator) {
                int current_inner_row = innerIterator.row();
                two_rings.insert(current_inner_row);
            }
        }
        int j = 0;
        //Let's get P in kx3 relative to v
        Eigen::MatrixXd P(two_rings.size(), 3);
        for (auto index : two_rings) {
            P.row(j) = V.row(index) - V.row(i);
            ++j;
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(P.transpose() * P);

        //two most principal directions (calling them the u- and v- directions) as in the readme
        Eigen::VectorXd u = eig.eigenvectors().col(2);
        Eigen::VectorXd v = eig.eigenvectors().col(1);

        // w-height= least principal direction
        Eigen::VectorXd w = eig.eigenvectors().col(0);

        // compute quadratic surface coefficients
        Eigen::VectorXd B = P * w; // project P to w coordinates

        // Let's create the quad_term "u, v, u*u, uv, v*v" matrix to compute igl::pinv
        Eigen::MatrixXd quad_term(P.rows(), 5);
        quad_term.col(0) = P * u;
        quad_term.col(1) = P * v;
        quad_term.col(2) = quad_term.col(0).cwiseProduct(quad_term.col(0));
        quad_term.col(3) = quad_term.col(0).cwiseProduct(quad_term.col(1));
        quad_term.col(4) = quad_term.col(1).cwiseProduct(quad_term.col(1));

        //solving
        Eigen::MatrixXd X;
        igl::pinv(quad_term, X);
        Eigen::VectorXd A = X * B;

        // compute shape operator
        // equation 39
        double E = 1.0 + A(0) * A(0);
        double F_ = A(0) * A(1);
        double G = 1.0 + A(1) * A(1);
        Eigen::Matrix2d uppercase;
        uppercase << E, F_, F_, G;

        //---
        double e = 2.0 * A(2) / sqrt(A(0) * A(0) + 1.0 + A(1) * A(1));
        double f = A(3) / sqrt(A(0) * A(0) + 1.0 + A(1) * A(1));
        double g = 2.0 * A(4) / sqrt(A(0) * A(0) + 1.0 + A(1) * A(1));
        Eigen::Matrix2d lowercase;
        lowercase << e, f, f, g;


        Eigen::Matrix2d S = -lowercase * uppercase.inverse();

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigS(S);

        K1(i) = eigS.eigenvalues()(1);
        K2(i) = eigS.eigenvalues()(0);
        D1.row(i) = eigS.eigenvectors()(1, 1) * u + eigS.eigenvectors()(1, 0) * v;
        D2.row(i) = eigS.eigenvectors()(0, 1) * u + eigS.eigenvectors()(0, 0) * v;
    }

}
