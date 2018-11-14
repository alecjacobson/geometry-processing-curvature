#include "../include/principal_curvatures.h"
#include <cmath>
#include <Eigen/Eigenvalues>
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/pinv.h>
#include <igl/sort.h>

void minmax(const Eigen::VectorXd eigenvalues,
    double & min,
    double & max,
    int & imin,
    int & imax) 
{
    imin = 0;
    imax = 0;

    max = eigenvalues[0];
    min = eigenvalues[0];

    for (int i = 1; i < eigenvalues.size(); i++) 
    {
        if (eigenvalues[i] > max) 
        {
            max = eigenvalues[i];
            imax = i;
        }
        else if (eigenvalues[i] < min) 
        {
            min = eigenvalues[i];
            imin = i;
        }
    }
}

void shape_operator(
    const Eigen::VectorXd & a,
    Eigen::MatrixXd & S) 
{
    double E = 1 + std::pow(a[0], 2);
    double F = a[0] * a[1];
    double G = 1 + std::pow(a[1], 2);

    double d = E + G - 1;

    double e = 2 * a[2] / sqrt(d);
    double f = a[3] / sqrt(d);
    double g = 2 * a[4] / sqrt(d);

    Eigen::MatrixXd h(2, 2), H(2, 2);
    h << e, f, f, g;
    H << E, F, F, G;
    
    S.resize(2, 2);
    S = -h * H.inverse();
}

void principal_curvatures(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & D1,
    Eigen::MatrixXd & D2,
    Eigen::VectorXd & K1,
    Eigen::VectorXd & K2)
{
    K1.resize(V.rows());
    K2.resize(V.rows());
    D1.resize(V.rows(), 3);
    D2.resize(V.rows(), 3);

    Eigen::SparseMatrix<double> A;
    igl::adjacency_matrix(F, A);

    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);

    for (int k = 0; k < A.outerSize(); k++)
    {
        std::vector<int> vs;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        {
            vs.push_back(it.row());
        }

        Eigen::MatrixXd P(vs.size(), 3);
        for (int v = 0; v < vs.size(); v++) 
        {
            P.row(v) = V.row(vs[v]) - V.row(k);
        }

        // compute eigenvalues and eigenvectors
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        Eigen::MatrixXd PP = P.transpose() * P;
        es.compute(PP);
        Eigen::MatrixXd PP_eigs, IPP;

        // choose the eigenvectors with highest eigenvalues
        double Pmin, Pmax;
        int iPmin, iPmax;
        minmax(es.eigenvalues(), Pmin, Pmax, iPmin, iPmax);
        Eigen::VectorXd Peig1 = es.eigenvectors().col(iPmax);
        Eigen::VectorXd Peig2 = es.eigenvectors().col(3 - iPmin - iPmax);

        // change of basis
        Eigen::MatrixXd principal_dir(3, 2);
        principal_dir << Peig1 / Peig1.norm(), Peig2 / Peig2.norm();
        Eigen::MatrixXd uv = P * principal_dir;
        Eigen::VectorXd w = P * (N.row(k).transpose() / N.row(k).norm());

        // least squares approximation
        Eigen::MatrixXd C(P.rows(), 5), X, S;
        C << uv.col(0),
            uv.col(1),
            uv.col(0).array().square(),
            uv.col(0).array() * uv.col(1).array(),
            uv.col(1).array().square();
        igl::pinv(C, X);
        Eigen::VectorXd a = X * w;

        // eigendecomposition of shape operator
        shape_operator(a, S);
        es.compute(S);
        Eigen::VectorXd S_eigvals;
        Eigen::MatrixXd S_eigvecs;

        double Smin, Smax;
        int iSmin, iSmax;
        minmax(es.eigenvalues(), Smin, Smax, iSmin, iSmax);

        K1[k] = Smax;
        K2[k] = Smin;
        
        D1.row(k) = principal_dir * es.eigenvectors().col(iSmax);
        D2.row(k) = principal_dir * es.eigenvectors().col(iSmin);
        D1.row(k) /= D1.row(k).norm();
        D2.row(k) /= D2.row(k).norm();
    }
}
