#include "../include/principal_curvatures.h"
#include "igl/adjacency_matrix.h"
#include "igl/slice.h"
#include "igl/pinv.h"
#include "igl/per_vertex_normals.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/Geometry>


void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
    
    //Compute number of vertices
    int numV = F.maxCoeff() + 1;
    
    //Reshape principle values and directions
    K1 = Eigen::VectorXd::Zero(V.rows());
    K2 = Eigen::VectorXd::Zero(V.rows());
    D1 = Eigen::MatrixXd::Zero(V.rows(),3);
    D2 = Eigen::MatrixXd::Zero(V.rows(),3);
    
    //Adj is the adjacency matrix
    //TwoRing is matrix storing the two ring of each vertex
    //Id is the identity matrix
    Eigen::SparseMatrix<double> Adj, TwoRing, Id;
    Adj.resize(numV,numV);
    TwoRing.resize(numV,numV);
    Id.resize(numV,numV);
    
    //Set identity
    Id.setIdentity();
    
    //Compute adjacency matrix
    igl::adjacency_matrix(F,Adj);
    
    //Compute the 1-ring (This is somewhat unnecessary, by construction,
    //the adjacency matrix of triangle meshes allow you to travel to any adjacent vertex
    //in 2 steps and travel from and back to the same vertex in 2 steps
    //i.e. 0-ring and 1-ring implicitly contained in the 2-ring
    Adj = Adj + Id;
    
    //Compute the 2-ring
    TwoRing = Adj  * Adj;
    
    //P as defined in README
    //N is the normals
    Eigen::MatrixXd P, N;
    
    //UVcoeffs are the u,v eigenvectors as in README
    //Wcoeffs is the w eigenvector as in README
    
    //AMat used for computing coefficients of the height-field surface
    //Ainv is the pseudo inverse of AMat.
    //sols is the coefficients
    Eigen::MatrixXd UVcoeffs, Wcoeffs, AMat, Ainv, sols;
    
    //Ed,Fd,Gd define the elements of the first fundamental form
    double Ed,Fd,Gd;
    //ed,fd,gd define the elements of the second fundamental form
    double ed,fd,gd;
    //LMat, second fundamental form
    //RMat, first fundamental form
    //S is a product of the first and second fundamental form
    //eigsVec are the eigenvectors of S
    Eigen::Matrix<double,2,2> LMat, RMat,S, eigsVec;
    //eigsVal are the eigenvalues of S
    Eigen::Vector2d eigsVal;
    
    //es_new used as eigensolver for S
    Eigen::EigenSolver<Eigen::Matrix<double,2,2>> es_new;

    //es used as eigensolver for Pt * P
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3>> es;
    
    //numEls stores number of element for resizing P
    int numEls;
    
    //maxInd, minInd used to reference maximum and minimum principle values
    int maxInd, minInd;
    //Used to update P
    int counter;
    

    //Compute the normals
    igl::per_vertex_normals(V,F,N);
    
    //Loop through each vertex and compute K1,K2,D1,D2
    for (int vNo = 0; vNo < numV; vNo ++) {

        //Extract the two-ring neighbourhood for each vertex
        numEls = TwoRing.innerVector(vNo).nonZeros();
        
        //resize P
        P.resize(numEls,3);
        
        //reset counter
        counter = 0;
        
        //Update P with vertex positions
        for (Eigen::SparseMatrix<double>::InnerIterator it(TwoRing,vNo); it; ++it) {
            P.row(counter++) = V.row(it.row()) - V.row(vNo);
            
        }
        
        //Compute eigenvalues and eigenvectors, selfadjoint sorts them.
        //Positive semidefinite so, eigenvalues are all positive
        es.compute(P.transpose() * P);
        //Extract the top 2 principal directions
        UVcoeffs = P * es.eigenvectors().rightCols(2);
        
        //Align the least principal direction with the normal
        //In principle, they should be equal
        if (es.eigenvectors().col(0).dot(N.row(vNo)) < 0) {
             Wcoeffs = -P * es.eigenvectors().col(0);
        } else {
             Wcoeffs = P * es.eigenvectors().col(0);
        }

        //Construct the A matrix
        AMat.resize(numEls,5);
        AMat.col(0).array() = UVcoeffs.col(0).array();
        AMat.col(1).array() = UVcoeffs.col(1).array();
        AMat.col(2).array() = UVcoeffs.col(0).array().square();
        AMat.col(3).array() =UVcoeffs.col(0).array() * UVcoeffs.col(1).array();
        AMat.col(4).array() = UVcoeffs.col(1).array().square();
        
        //Compute its pseudo inverse
        igl::pinv(AMat,Ainv);
        
        //Compute the coefficients
        sols = Ainv * Wcoeffs;
        
        //Compute the elements of the fundamental forms
        Ed = 1.0 + pow(sols(0),2.0);
        Fd = sols(0) * sols(1);
        Gd = 1.0 + pow(sols(1),2.0);
        
        ed = 2.0 * sols(2) / sqrt(1.0 + pow(sols(0),2.0) + pow(sols(1),2.0));
        fd = sols(3) / sqrt(1.0 + pow(sols(0),2.0) + pow(sols(1),2.0));
        gd = 2.0 * sols(4) / sqrt(1.0 + pow(sols(0),2.0) + pow(sols(1),2.0));
        
        //Construct the fundamental forms
        LMat<< ed,fd,fd,gd;
        RMat<< Ed,Fd,Fd,Gd;
        
        //Compute S
        S = -LMat * RMat.inverse();
        
        //Compute eigenvectors and eigenvalues of S.
        es_new.compute(S);
        
        //Extract the real components.
        eigsVal.array() = es_new.eigenvalues().real().array();
        eigsVec = es_new.eigenvectors().real();

        //Compute consistent index to min and max eigenvalues.
        if ((eigsVal(0)) > (eigsVal(1))) {
            minInd = 1;
            maxInd = 0;
            
        } else {
            minInd = 0;
            maxInd = 1;
        }
        
        //Set K1, K2 to be max and min eigenvalues respectively
        K1(vNo) = eigsVal(maxInd);
        K2(vNo) = eigsVal(minInd);
        
        //Lift the eigenvectors to 3D
        D1.row(vNo) = es.eigenvectors().rightCols(2) * eigsVec.col(maxInd);
        D2.row(vNo) = es.eigenvectors().rightCols(2) * eigsVec.col(minInd);
        
    }

}
