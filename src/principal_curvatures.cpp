#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/sort.h>
#include <igl/pinv.h>
#include <Eigen/SVD>
#include <igl/eigs.h>

void principal_curvatures(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    Eigen::MatrixXd &D1,
    Eigen::MatrixXd &D2,
    Eigen::VectorXd &K1,
    Eigen::VectorXd &K2)
{
    // Replace with your code
    K1 = Eigen::VectorXd::Zero(V.rows());
    K2 = Eigen::VectorXd::Zero(V.rows());
    D1 = Eigen::MatrixXd::Zero(V.rows(), 3);
    D2 = Eigen::MatrixXd::Zero(V.rows(), 3);

    Eigen::SparseMatrix<double> A, A2;
    igl::adjacency_matrix(F, A);
    //Compute 2 ring
    A2 = A * A;

    for (int i = 0; i < V.rows(); i++)
    {
        Eigen::MatrixXd P;
        int k = A2.col(i).nonZeros();
        P.resize(k, 3);
        int index = 0;

        //grab all other vertices to P in 2 rings
        for (Eigen::SparseMatrix<double>::InnerIterator it(A2, i); it; ++it)
        {
            P.row(index) = V.row(it.row()) - V.row(i);
            index++;
        }

        //SVD decomposition of P
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(P.transpose() * P, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd U = svd.matrixU();

        //get two most principle directions and one least principle direction
        Eigen::VectorXd u = U.col(0) / U.col(0).norm();
        Eigen::VectorXd v = U.col(1) / U.col(1).norm();
        Eigen::VectorXd w = U.col(2) / U.col(2).norm();

        //maxtrix for w
        Eigen::MatrixXd C(k, 5), Ci;
        Eigen::VectorXd v1, v2, v3, v4, v5;
        v1 = P * u;
        v2 = P * v;
        v3 = (P * u).cwiseProduct(P * u);
        v4 = (P * u).cwiseProduct(P * v);
        v5 = (P * v).cwiseProduct(P * v);
        C << v1, v2, v3, v4, v5;

        //solve for a1-a5
        igl::pinv(C, Ci);
        Eigen::VectorXd a = Ci * P * w;

        //get S coefficients
        double E, F, G, e, f, g;
        E = 1.0 + a(0) * a(0);
        F = a(0) * a(1);
        G = 1 + a(1) * a(1);
        e = (2 * a(2)) / sqrt(a(0) * a(0) + 1 + a(1) * a(1));
        f = (a(3)) / sqrt(a(0) * a(0) + 1 + a(1) * a(1));
        g = (2 * a(4)) / sqrt(a(0) * a(0) + 1 + a(1) * a(1));

        Eigen::Matrix2d first, second;
        first << e, f,
                f, g;

        second << E, F,
                F, G;

        //perform enigen decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(-first * second.inverse());

        K1(i) = solver.eigenvalues()(0);
        K2(i) = solver.eigenvalues()(1);
        D1.row(i) = solver.eigenvectors()(0, 0) * u + solver.eigenvectors()(0, 1) * v;
        D2.row(i) = solver.eigenvectors()(1, 0) * u + solver.eigenvectors()(1, 1) * v;
    }
}
