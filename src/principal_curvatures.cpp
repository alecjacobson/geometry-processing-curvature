#include "../include/principal_curvatures.h"
#include "igl/adjacency_matrix.h"
#include "igl/slice_mask.h"
#include "igl/pinv.h"
#include "igl/sort.h"
#include <Eigen/Eigenvalues>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
    int n = V.rows();
    Eigen::MatrixXd Eye(n,n), DenseA(n,n), DenseA2(n,n);
    Eigen::SparseMatrix<double> A, A2;
    
    igl::adjacency_matrix(F,A);
    Eye = Eigen::MatrixXd::Identity(V.rows(),V.rows());
    A2 = A*A;
    DenseA = Eigen::MatrixXd(A);
    DenseA2 = Eigen::MatrixXd(A2);
    
    //col i of two_ring will have a 1 for every vertex in the
    //two ring of vi (including vi), and 0 otherwise
    Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> two_ring(n,n);
    two_ring = (Eye.array() + DenseA.array() + DenseA2.array()) > 0;

    K1 = Eigen::VectorXd::Zero(V.rows());
    K2 = Eigen::VectorXd::Zero(V.rows());
    D1 = Eigen::MatrixXd::Zero(V.rows(),3);
    D2 = Eigen::MatrixXd::Zero(V.rows(),3);
    
    for(int i = 0; i < V.rows(); i++){
        Eigen::MatrixXd P;
        igl::slice_mask(V,two_ring.col(i),1,P);
        
        P = P.rowwise() - V.row(i);
        
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::MatrixXd scores = svd.matrixU() * svd.singularValues().asDiagonal();
        
        
        Eigen::VectorXd u = scores.col(0);
        Eigen::VectorXd v = scores.col(1);
        Eigen::VectorXd w = scores.col(2);
        
        Eigen::MatrixXd LS_Matrix(P.rows(),5);
        LS_Matrix.col(0) = u;
        LS_Matrix.col(1) = v;
        LS_Matrix.col(2) = (u.array().square()).matrix();
        LS_Matrix.col(3) = (u.array()*v.array()).matrix();
        LS_Matrix.col(4) = (v.array().square()).matrix();
        
        Eigen::MatrixXd MPinv;
        igl::pinv(LS_Matrix,MPinv);
        
        Eigen::VectorXd a = MPinv * w;
        
        double E = 1 + a(0)*a(0);
        double F = a(0)*a(1);
        double G = 1 + a(1)*a(1);
        double denominator = std::sqrt(a(0)*a(0) + 1 + a(1)*a(1));
        double e = 2*a(2)/denominator;
        double f = a(3)/denominator;
        double g = 2*a(4)/denominator;
        
        Eigen::Matrix2d S, S_left, S_right;
        S_left << e, f, f, g;
        S_right << E, F, F, G;
        
        S = - S_left * S_right.inverse();
        std::cout << S << std::endl;
        
        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(S);
        Eigen::Vector2d K;
        Eigen::Vector2i KI;
        Eigen::VectorXcd eivals = eigensolver.eigenvalues();
        
        igl::sort(eivals.real(),1,false,K,KI);
        std::cout << K << std::endl;
        K1(i) = K(0);
        K2(i) = K(1);
        D1.row(i) = svd.matrixV()*eigensolver.eigenvectors().col(KI(0)).real();
        D2.row(i) = svd.matrixV()*eigensolver.eigenvectors().col(KI(1)).real();
        
        
    }
}
