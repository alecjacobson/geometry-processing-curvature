#include "../include/principal_curvatures.h"
#include "igl/adjacency_matrix.h"
#include "igl/slice.h"
#include "igl/pinv.h"
#include "igl/per_vertex_normals.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

#include <iostream>
using namespace std;


void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
    //return;
  // Replace with your code
    int numV = F.maxCoeff() + 1;
    
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(numV);
    
    
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);
    
  Eigen::SparseMatrix<double> Adj, TwoRing, Sliced, Id;
    Adj.resize(numV,numV);
    TwoRing.resize(numV,numV);
    Id.resize(numV,numV);
    Id.setIdentity();
    igl::adjacency_matrix(F,Adj);
    

    Adj = Adj + Id;
    
    //TwoRing.resize(numV,numV);
    TwoRing = Adj  * Adj;
    
    //int numEls = TwoRing.nonZeros();
    //cout << "Non-Zeros: " << numEls << "\n";
    //cout << "Diagonal: " << Adj.diagonal() << endl;
    
    Eigen::MatrixXi rowInd;
    Eigen::MatrixXd P, N;
    Eigen::MatrixXd UVcoeffs, Wcoeffs, AMat, Ainv, sols;
    double Ed,Fd,Gd;
    double ed,fd,gd;
    Eigen::Matrix<double,2,2> LMat, RMat,S, eigsVec;
    Eigen::Vector2d eigsVal;
    Eigen::EigenSolver<Eigen::Matrix<double,2,2>> es_new;
    rowInd.resize(1,1);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;
    
    int numEls;
    int maxInd, minInd;
    int counter = 0;
    double dir;
    /*double max = 0;
    for (int k=0; k<TwoRing.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(TwoRing,k); it; ++it)
        {
            if (it.value() > max) {
                max = it.value();
            }
            
        }
    }
    
    cout << "Max: " << max << endl;
    return;*/
    
    igl::per_vertex_normals(V,F,N);
    
    //cout << "started" << endl;
    for (int vNo = 0; vNo < numV; vNo ++) {

        numEls = TwoRing.innerVector(vNo).nonZeros();
        
        cout << "Non-Zeros: " << numEls << "\n";
        P.resize(numEls,3);
        
        counter = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(TwoRing,vNo); it; ++it) {
            //cout << "Counter " << counter << ": Value: " << it.value() << ", Index: " << it.row() << endl;
            P.row(counter++) = V.row(it.row()) - V.row(vNo);
            
        }
        //cout << endl;
        
        //cout << P << "\n";
        
        es.compute(P.transpose() * P);
        //cout << "The first eigenvector of the 3x3 matrix of ones is:"
        //<< endl << es.eigenvalues() << endl;
        UVcoeffs = P * es.eigenvectors().rightCols(2);
        if (es.eigenvectors().col(0).dot(N.row(vNo)) < 0) {
             Wcoeffs = -P * es.eigenvectors().col(0);
        } else {
             Wcoeffs = P * es.eigenvectors().col(0);
        }
        //Hypothesis, need to keep the right hand rule
        
       
        
        //cout << "W: " << es.eigenvectors().col(0).transpose() << endl;
        AMat.resize(numEls,5);
        AMat.col(0).array() = UVcoeffs.col(0).array();
        AMat.col(1).array() = UVcoeffs.col(1).array();
        AMat.col(2).array() = UVcoeffs.col(0).array().square();
        AMat.col(3).array() =UVcoeffs.col(0).array() * UVcoeffs.col(1).array();
        AMat.col(4).array() = UVcoeffs.col(1).array().square();
        igl::pinv(AMat,Ainv);
        sols = Ainv * Wcoeffs;
        
        Ed = 1.0 + pow(sols(0),2.0);
        Fd = sols(0) * sols(1);
        Gd = 1.0 + pow(sols(1),2.0);
        
        ed = 2.0 * sols(2) / sqrt(1.0 + pow(sols(0),2.0) + pow(sols(1),2.0));
        fd = sols(3) / sqrt(1.0 + pow(sols(0),2.0) + pow(sols(1),2.0));
        gd = 2.0 * sols(4) / sqrt(1.0 + pow(sols(0),2.0) + pow(sols(1),2.0));
        
        LMat<< ed,fd,fd,gd;
        RMat<< Ed,Fd,Fd,Gd;
        //cout << LMat << endl;
        //cout << RMat << endl;
        S = -LMat * RMat.inverse();
        //cout << S << endl;
        es_new.compute(S);
        
        //cout << "The eigenvalues of A are: " << es_new.eigenvalues().transpose() << endl;
        
        eigsVal.array() = es_new.eigenvalues().real().array();
        eigsVec = es_new.eigenvectors().real();

        if ((eigsVal(0)) > (eigsVal(1))) {
            minInd = 1;
            maxInd = 0;
            
        } else {
            minInd = 0;
            maxInd = 1;
        }
        
        K1(vNo) = eigsVal(maxInd);
        K2(vNo) = eigsVal(minInd);
        
        D1.row(vNo) = es.eigenvectors().rightCols(2) * eigsVec.col(maxInd);
        D2.row(vNo) = es.eigenvectors().rightCols(2) * eigsVec.col(minInd);
        
        if (vNo %  1000 == 1) {
        cout << "Num: " << vNo << endl;
        }
        
    }
    //cout << "done" << endl;
}
