#include "../include/angle_defect.h"
#include "igl/squared_edge_lengths.h"
#include <math.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());
    Eigen::MatrixXd L;
    L.resizeLike(V);
    
    igl::squared_edge_lengths(V,F, L);
    
    int numF = F.rows();
    int numV = F.maxCoeff() + 1;
    
    Eigen::MatrixXd vals(3,2);
    
    //Makes it easier to reference other vertices
    vals(0,0) = 1;
    vals(0,1) = 2;
    vals(1,0) = 0;
    vals(1,1) = 2;
    vals(2,0) = 0;
    vals(2,1) = 1;
    
    double cosAng;
    for (int curV = 0; curV < numV; curV ++) {
        D(curV) = 2.0 * M_PI;
    }

    for (int curF = 0; curF < numF; curF ++) {
        
        for (int fNo = 0; fNo < 2; fNo ++){
        //Apply cosine law
            cosAng =L(curF,vals(fNo,0)) + L(curF,vals(fNo,1)) - L(curF,fNo);
            
            cosAng = cosAng / (2.0 * sqrt(L(curF,vals(fNo,0)) * L(curF,vals(fNo,1))));
            
            D(F(curF,fNo)) -= acos(cosAng);
        }
        
    }
    
}
