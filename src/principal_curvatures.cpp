#include "../include/principal_curvatures.h"
#include <Eigen/Eigenvalues>

#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>

#include <vector>

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


    // Gather  two-ring neighborhood of each vertex

    Eigen::SparseMatrix<int> A;
    igl::adjacency_matrix(F, A);

    // Compute 2-ring by matrix multiplication
    //  vertices in 1-ring of v also in 2-ring of v by nature of triangle mesh
    //  B(i,i) > 0, since every vertex is distance 2 from itself.
    Eigen::SparseMatrix<int> B;
    B = A*A;

    // Construct neighborhood distance matrix
    int k = 0;
    Eigen::MatrixXd P;

    // Eigen Decomposition solver
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver1;
    // Eigen-vectors and projected positions
    Eigen::MatrixXd W, PW;

    // Fitting height surface with quadratic polynomial
    Eigen::MatrixXd U, Uinv;
    Eigen::VectorXd x(5);

    // Shape operator `S`, foundamental forms `ff1`, `ff2`
    Eigen::Matrix2d ff1, ff2, S;
    double e,f,g,E,FF,G;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver2;

    // principle tangent 
    Eigen::Vector2d pt1, pt2;


    for (int i = 0; i < V.rows(); ++i) {

        // Construct `P`, distance from `v` for vertices in 2-ring of `v`

        P.resize(B.innerVector(i).nonZeros(), 3);
        k = 0;
        for (Eigen::SparseMatrix<int>::InnerIterator it(B, i); it; ++it) {
            // Visit nonzero element in i-th column
            P.row(k) = V.row(it.row()) - V.row(i);
            k += 1;
        }

        // PCA on `P`

        solver1.compute(P.transpose() * P);
        W = solver1.eigenvectors();
        PW = P * W;

        // Note `W` corresponds to eigenvalue in increasing order
        //      hence first principle component is W.col(2), i.e. u

        // Fit polynomial to height-field surface
        //  by solving a least squared problem for w=a₁u+a₂v+a₃u²+a₄uv+a₅v²
        //      i.e. find x = [a₁ a₂ a₃ a₄ a₅]^T by solving `Ux = b` where 
        //      `U` consists of value of `u,v,u²,uv,v²` and `b = PW.col(0)`

        U.resize(P.rows(), 5);
        U.col(0) = PW.col(2);
        U.col(1) = PW.col(1);
        U.col(2) = PW.col(2).array().square();
        U.col(3) = PW.col(2).array() * PW.col(1).array();
        U.col(4) = PW.col(1).array().square();

        igl::pinv(U, Uinv);
        x = Uinv * PW.col(0);

        // Compute shape operator `S`

        E = 1 + x(0)*x(0);
        FF = x(0)*x(1);
        G = 1 + x(1)*x(1);
        e = 2*x(2) / sqrt(E + G - 1);
        f =   x(3) / sqrt(E + G - 1);
        g = 2*x(4) / sqrt(E + G - 1);

        ff2 << e, f,
               f, g;
        ff1 << E, FF, 
               FF, G;
        S = -ff2 * ff1.inverse();

        // Eigen-decomposition on S
        //      eigenvalues are principle curvatures
        //      eigenvectors are principle tangent directions in (u,v) basis
        solver2.compute(S.transpose() * S);

        K1(i) = solver2.eigenvalues()(1);
        K2(i) = solver2.eigenvalues()(0);
        
        //  map pt = (a,b) \in span{u,v} to 3D: au + bv
        pt1 = solver2.eigenvectors().col(1);
        pt2 = solver2.eigenvectors().col(0);

        D1.row(i) = pt1(0) * W.col(2) + pt1(1) * W.col(1);
        D2.row(i) = pt2(0) * W.col(2) + pt2(1) * W.col(1);
    }
}
