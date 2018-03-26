#include "../include/mean_curvature.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/per_vertex_normals.h"
#include "igl/invert_diag.h"
#include <iostream>
using namespace std;
void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  H = Eigen::VectorXd::Zero(V.rows());
    
    int numV = F.maxCoeff() + 1;
    Eigen::SparseMatrix<double> L, M, M_inv;
    Eigen::MatrixXd N,Hn;
    N.resizeLike(V);
    L.resize(numV,numV);
    M.resize(numV,numV);
    M_inv.resize(numV,numV);
    Hn.resizeLike(V);

    
    igl::cotmatrix(V,F,L);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    igl::per_vertex_normals(V,F,N);
    igl::invert_diag(M,M_inv);
    
    Hn = M_inv * L * V;

    double curSign;
    //Check consistency with normals.
    for (int normNo = 0; normNo < numV; normNo ++) {
        curSign =Hn.row(normNo).dot(N.row(normNo));
        H(normNo) = Hn.row(normNo).norm();
        if (curSign > 0) {
            H(normNo) = -1.0 * H(normNo);
        }
        
        
        
    }
    
}
